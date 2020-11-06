// A program to calculate the Mutual information expansion (MIE) of macromolecules from trajectories 
// Copyright (C) 2020  Markus Fleck
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3 as 
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.


#define PRECISION double // define the compiled precision 
#define PRECISION4 double4 // used in the 2D histogram kernel for efficient memory fetches


#define MODFITNBINS 100 // used for finding the longest non-zero stretch for the circular dihedrals
#define N_STREAMS 32 // how many cuda streams to use
#define HISTOGRAM_THREAD_WORK_MULTIPLE 8 // increase to do more work per thread durig building the 2D histograms
#define HISTOGRAM_THREADS 512 // the number of threads per block during building the histograms
#define PLNP_THREADS 32
#define USE_SHARED_MEM_HISTOGRAMS // use the histo2D_shared_block kernel

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <vector>

#include "../util/util.h"
#include "../util/io/io.h"
#include "../util/types.h"

#include "../util/classes/Bat.h"
#include "../util/classes/Entropy_Matrix.h"
// TODO: proper inclusion
#include "PARENT_GPU_kernels.cu"



using namespace std;

/* Variable names of the dofs:
a _g suffix means that the variable is global. If there are 1000 atoms,
there are 3 * 1000 - 6 = 2994 dofs. Thus, a _g variable would range from
0 to 2993. A _l suffix means that the variable is local to the block.
E.g., if the block starts at dof 300 and dof_g = 400, then dof_l = 100.
A _t suffix (type) means that the variable either refers to only bonds, only 
angles or only dihedrals. In our example, there are 1000 - 1 = 999 bonds
and 1000 - 2 = 998 angles. So, if dof_g = 1050, then dof_gt = 1050 - 999 = 51.
The order is bonds before angles, angles before dihedrals.
 */
 
 
// this class represents a chunk of degrees of freedom (dofs) in CPU RAM which can be loaded to the GPU
class GPU_RAM_Block {
public:
  unsigned char type; // the type of the degrees of freedom the block contains, i. e. either bonds(0), angles(1) or dihedrals(2)
  unsigned int dof_id_start_g; // the global id of the first dof in the block
  unsigned int dof_id_end_g; // the global id of the last dof in the block (inclusive)
  unsigned int n_dofs; // the total number of dofs in the block
  size_t n_bytes; // the amount of data the block occupies in bytes
  PRECISION *cpu_ram_start; // a pointer to the data on the CPU (set upon object instantiation). Trajectories are stored contigeously for each dof
  PRECISION *gpu_ram_start; // a pointer to the data on the GPU (set only upon calling the deploy function). Trajectories are stored contigeously for each dof

  GPU_RAM_Block(PRECISION *cpu_ram_start, unsigned int dof_id_start_g,
                unsigned int dof_id_end_g, unsigned int n_frames,
                unsigned int n_dihedrals) {
    // initialize the variables
    this->cpu_ram_start = cpu_ram_start;
    this->type = get_dof_type_from_id(dof_id_start_g, n_dihedrals);
    this->dof_id_start_g = dof_id_start_g;
    this->dof_id_end_g = dof_id_end_g;
    n_dofs = dof_id_end_g - dof_id_start_g + 1; // note that dof_id_end_g is inclusive 
    n_bytes = n_dofs * n_frames * sizeof(PRECISION); // each dof features a trajectory with n_frames, each frame of each dof represented by a single number of PRECISION
  }

  // This function loads the data of the GPU_RAM_Block from the CPU to the GPU. The pointer to the data on the CPU is set upon object instantiation, the pointer for the GPU is handed over as a function parameter
  void deploy(PRECISION *gpu_ram_start) {
    gpuErrchk(cudaMemcpy(gpu_ram_start, cpu_ram_start, n_bytes,
                         cudaMemcpyHostToDevice)); // copy the data to one of the two GPU RAM banks, which have already been alloctated
    this->gpu_ram_start = gpu_ram_start; // save the location of the data on the GPU, which will be one of the two GPU RAM banks
  }
};


// this class represents a chunk of degrees of freedom (dofs) which can be deployed to a CPU RAM bank. Extrema as well as 1D entropies are calculated upon deployement.
class CPU_RAM_Block {
public:
  unsigned int dof_id_start_g; // the global id of the first dof in the block
  unsigned int dof_id_end_g; // the global id of the last dof in the block (inclusive)
  unsigned int n_dofs; // the total number of dofs in the block
  int type_id_start_g[3] = {-1, -1, -1}; // the global ids of the first dof in the block for each dof type (i. e. bond, angle, dihedral)
  int type_id_end_g[3] = {-1, -1, -1}; // the global ids of the last dof (inclusive) in the block for each dof type (i. e. bond, angle, dihedral)
  unsigned int type_n_dofs[3] = {0, 0, 0}; // the number of dofs in the block for each type
  unsigned int gpu_ram_blocks_per_type[3]; // the number of GPU_RAM_Blocks needed for each dof type
  unsigned int n_dihedrals; // the number of total number of dihedrals in the molecule. From this number, the number of bonds(==n_dihedrals + 2), angles(==n_dihedrals + 1) and atoms(==n_dihedrals + 3) can be calculated  
  unsigned int gpu_dofs_per_block; // the number of dofs which fits in a GPU_RAM_Block
  unsigned int n_frames; // the number of frames of the trajectory 
  unsigned int n_frames_padded; // the padded number of frames. The padding is done to a multiple of 4 and used in the histo2D_shared_block kernel for efficient data loading
  vector<GPU_RAM_Block> blocks; // the GPU_RAM_Blocks associated to this CPU_RAM_Block
  PRECISION *block_start; // the location of the trajectories. Trajectories are stored contigeously for each dof
  PRECISION *minima; // the location of the minima for each dof
  PRECISION *maxima; // the location of the maxima for each dof
  bool extrema_calculated = false; // a flag indicating if minima, maxima and 1D entropy values have already been calculated
  unsigned int n_bins; // the number of bins used for building the histograms
  PRECISION *result_entropy1D; // the location to store the 1D entropy results
  PRECISION *type_addr[3]; // the start of the trajectories for each dof type 
  PRECISION *bonds; // == type_addr[0]
  PRECISION *angles; // == type_addr[1]
  PRECISION *dihedrals; // == type_addr[2]
  Bat *bat; // the BAT trajectory class to load the data from hard disk

  CPU_RAM_Block(unsigned int dof_id_start_g, unsigned int dof_id_end_g,
                unsigned int gpu_dofs_per_block, Bat *bat, PRECISION *minima,
                PRECISION *maxima, unsigned int n_bins,
                PRECISION *result_entropy1D) {

    // initilize the variables            
    this->dof_id_start_g = dof_id_start_g;
    this->dof_id_end_g = dof_id_end_g;
    this->n_dofs = dof_id_end_g - dof_id_start_g + 1; // note that dof_id_end_g is inclusive 
    this->n_dihedrals = bat->get_n_dihedrals();
    this->gpu_dofs_per_block = gpu_dofs_per_block;
    this->n_frames = bat->get_n_frames();
    this->n_frames_padded = bat->get_n_frames_padded(4);
    this->minima = minima;
    this->maxima = maxima;
    this->n_bins = n_bins;
    this->result_entropy1D = result_entropy1D;
    this->bat = bat;

    for (unsigned short type = get_dof_type_from_id(dof_id_start_g, n_dihedrals);
         type <= get_dof_type_from_id(dof_id_end_g, n_dihedrals); type++) { // for all dof types contained in the CPU_RAM_Block
         
      if (dof_id_start_g < get_min_id_for_type(type, n_dihedrals)) { // if the first dof id of the CPU_RAM_Block is smaller than the first id of the type 
        type_id_start_g[type] = get_min_id_for_type(type, n_dihedrals); // the block contains the first dof of the type
      } else {
        type_id_start_g[type] = dof_id_start_g; // otherwise the first dof of this type in the block is also the very first dof in the block
      }

      if (dof_id_end_g > get_max_id_for_type(type, n_dihedrals)) { // if the last dof id of the CPU_RAM_Block is larger than the last id of the type
        type_id_end_g[type] = get_max_id_for_type(type, n_dihedrals); // the block contains the last dof of the type
      } else {
        type_id_end_g[type] = dof_id_end_g; // otherwise the last dof of this type in the block is also the very last dof in the block
      }
      type_n_dofs[type] = type_id_end_g[type] - type_id_start_g[type] + 1; // set the number of dofs for each type in the block

      gpu_ram_blocks_per_type[type] = (type_n_dofs[type] + (gpu_dofs_per_block - 1) ) / gpu_dofs_per_block; // calculate how many GPU_RAM_Blocks are needed for each dof type. There is an additional GPU_RAM_Block used if the devision has a remainder
    }
  }

  // This function loads the dofs from the hard disk into this CPU RAM block and calculates minima, maxima and 1D entropy values
  void deploy(PRECISION *block_start) { // load the data to the RAM bank with the adrress block_start

    this->block_start = block_start; // save the RAM bank where the data is loaded to
    type_addr[TYPE_B] = block_start; // save the start of the bonds
    type_addr[TYPE_A] = block_start + type_n_dofs[TYPE_B] * n_frames_padded; // save the start of the angles
    type_addr[TYPE_D] =
        block_start + (type_n_dofs[TYPE_B] + type_n_dofs[TYPE_A]) * n_frames_padded; // save the start of the dihedrals
    bonds = type_addr[TYPE_B]; // create an alias for the bonds
    angles = type_addr[TYPE_A]; // create an alias for the angles
    dihedrals = type_addr[TYPE_D]; // create an alias for the dihedrals

    bat->load_dofs(type_addr, type_id_start_g, type_id_end_g, 4); // load the dofs contigeously (order bonds, angles, dihedrals) and with a padding of 4

    modfit_dihedrals(); // apply modiftting to the dihedrals, see the function below
    if (!extrema_calculated) { // if the extrema have not already been calculated 
      calculate_extrema(); // calculate minima and maxima 
      calculate_entropy1D(); // also calculate 1D entropy values
      extrema_calculated = true; // set the flag which indicates that minima, maxima as well as 1D entropy values have already been calculated
    }

    blocks.clear(); // clear the GPU_RAM_Blocks to reinitialize them
    for (unsigned char type = get_dof_type_from_id(dof_id_start_g, n_dihedrals); // for all dof types in the CPU_RAM_Block
         type <= get_dof_type_from_id(dof_id_end_g, n_dihedrals); type++) {
      for (unsigned int i = 0; i < gpu_ram_blocks_per_type[type]; i++) { // and for all GPU_RAM_Blocks of this type
        unsigned int gpu_block_id_start_g =
            type_id_start_g[type] + i * gpu_dofs_per_block; // calculate the id of the first dof in the GPU_RAM_Block to be created
        unsigned int gpu_block_id_end_g =
            type_id_start_g[type] + (i + 1) * gpu_dofs_per_block - 1; // calculate the id of the last dof in the GPU_RAM_Block to be created
        if (int(gpu_block_id_end_g) > type_id_end_g[type])
          gpu_block_id_end_g = type_id_end_g[type]; // if gpu_block_id_end_g refers to a different dof type than gpu_block_id_start_g, set gpu_block_id_end_g to the last id of the current dof type 
        PRECISION *cpu_ram_start =
            block_start + (gpu_block_id_start_g - dof_id_start_g) * n_frames_padded; // caluclate the address of the first dof of the current GPU_RAM_Block in the CPU RAM
        blocks.push_back(*new GPU_RAM_Block(cpu_ram_start, gpu_block_id_start_g,
                                            gpu_block_id_end_g, n_frames_padded,
                                            n_dihedrals)); // create the GPU_RAM_Block and add it to the vector blocks, containing all GPU_RAM_Blocks 
      }
    }
  }

    //this function finds the longest non-zero stretch of the trajectory for each dihedral and shifts the circular dihedral coordinates
    //so that the longest non-zero stretch has coordinate value 0. This leads to efficiently filled bins and improves numerical accuracy 
  void modfit_dihedrals() {
    if (type_n_dofs[TYPE_D] == 0)
      return; // return if there are no dihedrals in the CPU_RAM_Block
    const PRECISION pi = acos(-1); // get the value of pi from the library 

#pragma omp parallel
    {
      PRECISION modFit, binsize;
      int longestZeroStretch, longestZeroStretchPos, currentZeroStretch,
          currentZeroStretchPos;
      bool zeroExists;
      long long int histo[MODFITNBINS];

#pragma omp for
      for (unsigned int j = 0; j < type_n_dofs[TYPE_D]; j++) { // for all dihedrals
        // first build a histogram of the dihedral values over the trajectory
        for (int k = 0; k < MODFITNBINS; k++)
          histo[k] = 0;

        binsize = (2 * pi +
                   5e-9 * (sizeof(PRECISION) == sizeof(float) ? 100000 : 1)) /
                  MODFITNBINS; // divide 2 pi into MODFITNBINS bins, making sure that we cover a little more than 2 pi (depending on the precision used for the calculations) 
        for (unsigned int i = 0; i < n_frames; i++) // bin the dihedrals
          histo[int((dihedrals[j * n_frames_padded + i]) / binsize)] += 1;

        zeroExists = false;
        for (int k = 0; k < MODFITNBINS; k++) // check if there are empty bins
          zeroExists = zeroExists || (histo[k] == 0);

        if (zeroExists) { // if any of the bins of the histogram is empty find
                          // the longest consecutive stretch of empty bins
          longestZeroStretch = 0;
          currentZeroStretch = 0;
          longestZeroStretchPos = -1;
          for (int k = 0; k < 2 * MODFITNBINS;
               k++) {                // for all bins of the histogram
            int l = k % MODFITNBINS; // taking car of zero stretches which span
                                     // the circular boundaries
            if ((currentZeroStretch == 0) &&
                (histo[l] == 0)) { // find and save a beginning zero stretch
              currentZeroStretch = 1;
              currentZeroStretchPos = k;
            }
            else if ((currentZeroStretch > 0) && (histo[l] == 0)) {
              currentZeroStretch += 1;
            }
            else if ((currentZeroStretch > 0) &&
                (histo[l] != 0)) { // and the end of it. If it is currently the
                                   // longest zero stretch, save it
              if (currentZeroStretch > longestZeroStretch) {
                longestZeroStretch = currentZeroStretch;
                longestZeroStretchPos = currentZeroStretchPos;
              }
              currentZeroStretch = 0;
            }
          }
        } else { // if none of the bins is empty
          longestZeroStretchPos =
              0; // misuse the zeroStretch variables for determining the minimally filled bin
          longestZeroStretch = histo[0];
          for (int k = 0; k < MODFITNBINS; k++) {
            if (histo[k] < longestZeroStretch) {
              longestZeroStretch = histo[k];
              longestZeroStretchPos = k;
            }
          }
        }
        modFit = 2 * pi - (longestZeroStretchPos + 0.5) *
                              binsize; // calculate the shift to put the zero
                                       // stretch to the 2pi end
        for (unsigned int k = 0; k < n_frames; k++) {
          dihedrals[j * n_frames_padded + k] =
              dihedrals[j * n_frames_padded + k] + modFit -
              2 * pi *
                  int((dihedrals[j * n_frames_padded + k] + modFit) /
                      (2 * pi)); // and apply it taking care of circularity
        }
      }
    }
  }

  // calculate the extrema of the dofs for binning
  void calculate_extrema() {
#pragma omp parallel
    {
      PRECISION tmpMin, tmpMax;

#pragma omp for
      for (unsigned int j = 0; j < n_dofs; j++) { // for all dofs
        tmpMax = block_start[j * n_frames_padded]; // set the temporary minimum/maximum to the coordinate value of the first frame
        tmpMin = block_start[j * n_frames_padded];
        for (unsigned int i = 1; i < n_frames; i++) { // the for all frames
          if (block_start[j * n_frames_padded + i] > tmpMax) { // set the temporary maximum/minimum the a new value if the coordinates value is larger/smaller in the current frame  
            tmpMax = block_start[j * n_frames_padded + i];
          }
          if (block_start[j * n_frames_padded + i] < tmpMin) {
            tmpMin = block_start[j * n_frames_padded +
                                 i]; 
          }
        }
        if ((tmpMin < 0.0) || (tmpMax < 0.0)) { //check for errors
          cerr << "ERROR: Degree of freedom " << dof_id_start_g + j
               << " is smaller than 0.0" << endl;
          exit(EXIT_FAILURE);
        }
        tmpMin -= 5e-9 * (sizeof(PRECISION) == sizeof(float)
                              ? 1e5
                              : 1); // and increase the boundaries a tiny bit depending on the precision used
        tmpMax += 5e-9 * (sizeof(PRECISION) == sizeof(float) ? 1e5 : 1);
        if (tmpMin < 0.0) {
          tmpMin = 0.0;
        }
        if ((tmpMax - tmpMin) < 1.0e-4) {
          tmpMax += 0.05;
        }
        minima[dof_id_start_g + j] = tmpMin; // and store the calculated values
        maxima[dof_id_start_g + j] = tmpMax;
      }
    }
  }

  void calculate_entropy1D() {
#pragma omp parallel
    {
      PRECISION binsize, probDens, binval, plnpsum, Jac;
      long long int histo[n_bins];
      int occupbins;
#pragma omp for
      for (unsigned int j = dof_id_start_g; j <= dof_id_end_g;
           j++) { // for all dofs (using omp threads)
        for (unsigned int k = 0; k < n_bins; k++) {
          histo[k] = 0; // initialize a histogram with zeros
        }
        binsize = (maxima[j] - minima[j]) /
                  n_bins; // calculate the size of the bins
        for (unsigned int i = 0; i < n_frames;
             i++) { // and fill the histogram using all frames of the trajectory
          histo[int(
              (block_start[(j - dof_id_start_g) * n_frames_padded + i] - minima[j]) /
              binsize)] += 1;
        }
        occupbins = 0; // then use the histogram to calculate the (discretized)
                       // entropy, taking care of the Jacobian
        plnpsum = 0;
        binval = minima[j] + (binsize / 2.0);
        Jac = 1;
        for (unsigned int k = 0; k < n_bins; k++) {
          switch (get_dof_type_from_id(j, n_dihedrals)) {
          case TYPE_B:
            Jac = binval * binval;
            break;
          case TYPE_A:
            Jac = sin(binval);
            break;
          case TYPE_D:
            break;
          }
          probDens = histo[k] / (n_frames * binsize * Jac);
          if (probDens > 0) {
            plnpsum = plnpsum + Jac * probDens * log(probDens);
            occupbins = occupbins + 1;
          }
          binval += binsize;
        }
        plnpsum = -plnpsum * binsize;
        result_entropy1D[j] =
            plnpsum +
            (occupbins - 1.0) /
                (2.0 * n_frames); // and apply Herzel entropy unbiasing
      }
    }
  }
};

// this class provides the memory organization on the GPU via two RAM banks
class GPU_RAM_Layout {
public:
  char *gpu_ram_start; // the address of the data on the GPU
  size_t gpu_n_bytes; // the size of the data on the GPU
  unsigned int dofs_per_block; // how many dofs fit into one ram bank, considering the additional storage for the (temporary) results
  PRECISION *dof_block_1; // the address of the first ram bank
  PRECISION *dof_block_2; // the address of the second ram bank
  PRECISION *result; // the address of the 2D entropy entropy results
  unsigned int *occupied_bins; // the address of the result which stores the number of occupied bins in the histograms for Herzel entropy unbiasing
  unsigned int *histograms; // the address of the histograms to be build

  GPU_RAM_Layout(unsigned int n_frames, unsigned int n_bins,
                 size_t gpu_n_bytes, unsigned int n_dofs_total) {
    // calculate the maximum number of dofs (for one of the two dof_blocks) so
    // that everything still fits into GPU RAM. The formula is:  gpu_n_bytes == 2 * dofs_per_block * n_frames * sizeof(PRECISION) + dofs_per_block^2 * ( (n_bins^2 + 1) * sizeof(unsigned int) + sizeof(PRECISION) )
    double a = 2 * n_frames * sizeof(PRECISION);
    double b = sizeof(PRECISION) + sizeof(unsigned int) * (n_bins * n_bins + 1);
    dofs_per_block =
        (unsigned int)((-a / 2 + sqrt(a * a / 4 + gpu_n_bytes * b)) / b);
                 
    
                 
    // if all dofs fit into RAM, still set up two blocks to be consistent with
    // the algorithm
    if (2 * dofs_per_block > n_dofs_total) {
      dofs_per_block = (n_dofs_total + 1 )/ 2;
    }
    
    // calculate the number of bytes to be allocated on the GPU 
    gpu_n_bytes = (static_cast<size_t>(2) * dofs_per_block * n_frames + dofs_per_block * dofs_per_block) * sizeof(PRECISION) + ( dofs_per_block * dofs_per_block * (n_bins * n_bins + 1) ) * sizeof(unsigned int); 
    gpuErrchk(cudaMalloc((void **)&gpu_ram_start, gpu_n_bytes)); // and allocate it
    this->gpu_n_bytes = gpu_n_bytes; // store the number of bytes in the class

    // set the pointers for partitioning the GPU RAM according to the calculated
    // dofs_per_block
    dof_block_1 = (PRECISION *)gpu_ram_start; // set the starting address for GPU RAM bank 1 
    dof_block_2 = dof_block_1 + dofs_per_block * n_frames; // set the starting address for GPU RAM bank 1
    result = dof_block_2 + dofs_per_block * n_frames; // set the starting address for the 2D entropy results
    occupied_bins = (unsigned int *)(result + dofs_per_block * dofs_per_block); // set the starting address for the occupied bins result
    histograms = occupied_bins + dofs_per_block * dofs_per_block; // set the starting address for the histograms to be build
  }
};


// this class provides the memory organization on the CPU via two RAM banks
class CPU_RAM_Layout {
public:
  char *cpu_ram_start; // The address of the whole data block
  size_t cpu_n_bytes; // The number of bytes to use
  unsigned int dofs_per_block; // The number of dofs which each of the two blocks holds
  PRECISION *dof_block_1; // The address of the first dof block (==cpu_ram_start)
  PRECISION *dof_block_2; // The address of the second dof block
  PRECISION *result_entropy; // The address of the entropy results
  PRECISION *result_entropy1D; // The address of the 1D entropy results (==result_entropy)
  PRECISION *result_entropy1D_b; // The address of the 1D bonds entropy results (==result_entropy)
  PRECISION *result_entropy1D_a; // The address of the 1D angles entropy results
  PRECISION *result_entropy1D_d; // The address of the 1D dihedrals entropy results
  PRECISION *result_entropy2D; // The address of the 2D entropy results
  PRECISION *result_entropy2D_bb; // The address of the bond-bond 2D entropy results (==result_entropy2D)
  PRECISION *result_entropy2D_ba; // The address of the bond-angle 2D entropy results
  PRECISION *result_entropy2D_bd; // The address of the bond-dihedral 2D entropy results
  PRECISION *result_entropy2D_aa; // The address of the angle-angle 2D entropy results
  PRECISION *result_entropy2D_ad; // The address of the angle-dihedral 2D entropy results
  PRECISION *result_entropy2D_dd; // The address of the dihedral-dihedral 2D entropy results
  PRECISION *extrema; // The address of the extrema
  PRECISION *minima; // The address of the minima (==extrema)
  PRECISION *maxima; // The address of the maxima
  PRECISION *tmp_result_entropy; // The address for copying the entropy results from GPU to CPU RAM
  unsigned int *tmp_result_occupied_bins; // The address for copying the occupied_bins results from GPU to CPU RAM

  CPU_RAM_Layout(unsigned int n_frames, size_t cpu_n_bytes, unsigned int gpu_dofs_per_block,
                 unsigned int n_dihedrals) {
    unsigned int n_dofs_total = 3 * (n_dihedrals + 1); // == n_bonds + n_angles + n_dihedrals
    unsigned int n_bonds = n_dihedrals + 2;
    unsigned int n_angles = n_dihedrals + 1;
    
    // calculate the maximum number of dofs (for one of the two dof_blocks) so that everything still fits into CPU RAM. From the storage provided (cpu_n_bytes), subtract the storage used for the temporary results(tmp_result_entropy and tmp_result_occupied_bins),
    // the storage for the minima, maxima and 1D entropy values as well as the 2D entroy values. What remains can be used for the two dof blocks, where each dof block needs ( dofs_per_block * n_frames * sizeof(PRECISION) ) bytes.     
    dofs_per_block = ( cpu_n_bytes - gpu_dofs_per_block * gpu_dofs_per_block * ( sizeof(PRECISION) + sizeof(unsigned int) )  - (3 * n_dofs_total + n_dofs_total * ( n_dofs_total - 1 ) / 2) * sizeof(PRECISION) )
                   / double( 2 * n_frames * sizeof(PRECISION) );

    
    if (dofs_per_block < gpu_dofs_per_block) {
      cerr << "WARNING: You probably have a GPU with a lot of RAM but your CPU "
              "RAM is rather small. ";
      cerr << "I recommend to get more CPU RAM, as this should significantly "
              "enhance performance."
           << endl;
    }
    
    // if all dofs fit into RAM, still set up two blocks to be consistent with
    // the algorithm
    if (2 * dofs_per_block > n_dofs_total) dofs_per_block = (n_dofs_total + 1) / 2; //TODO: In this case, only load the trajectory once
    
    //calculate the actual number of CPU RAM bytes used
    cpu_n_bytes = (static_cast<size_t>(2) * dofs_per_block * n_frames + 3 * n_dofs_total + n_dofs_total * (n_dofs_total - 1) / 2) * sizeof(PRECISION) + gpu_dofs_per_block * gpu_dofs_per_block * ( sizeof(PRECISION) + sizeof(unsigned int) );
    this->cpu_n_bytes = cpu_n_bytes; // store this number
    cpu_ram_start = new char[cpu_n_bytes]; // allocate those bytes


    // set the pointers for partitioning the CPU RAM according to the calculated
    // dofs_per_block
    dof_block_1 = (PRECISION *) cpu_ram_start;
    dof_block_2 = dof_block_1 + dofs_per_block * n_frames;

    result_entropy = dof_block_2 + dofs_per_block * n_frames;
    result_entropy1D = result_entropy;
    result_entropy1D_b = result_entropy1D;
    result_entropy1D_a = result_entropy1D_b + n_bonds;
    result_entropy1D_d = result_entropy1D_a + n_angles;

    result_entropy2D = result_entropy1D + n_dofs_total;
    result_entropy2D_bb = result_entropy2D;
    result_entropy2D_ba = result_entropy2D_bb + n_bonds * (n_bonds - 1) / 2; // binomial formula for pairs 
    result_entropy2D_bd = result_entropy2D_ba + n_bonds * n_angles; // all combinations
    result_entropy2D_aa = result_entropy2D_bd + n_bonds * n_dihedrals; // all combinations
    result_entropy2D_ad = result_entropy2D_aa + n_angles * (n_angles - 1) / 2; // binomial formula for pairs
    result_entropy2D_dd = result_entropy2D_ad + n_angles * n_dihedrals; // all combinations

    extrema = result_entropy2D + n_dofs_total * (n_dofs_total - 1) / 2; // binomial formula for pairs
    minima = extrema;
    maxima = extrema + n_dofs_total;
    tmp_result_entropy = maxima + n_dofs_total;
    tmp_result_occupied_bins =
        (unsigned int *)(tmp_result_entropy + gpu_dofs_per_block * gpu_dofs_per_block);

  }
};

// This class manages the RAM on the CPU as well as on the GPU
class RAM {
public:
  GPU_RAM_Layout *gpu_ram_layout; // the RAM layout on the GPU
  CPU_RAM_Layout *cpu_ram_layout; // the RAM layout on the CPU
  unsigned int n_dihedrals; // the number of dihedrals of the molecule(s)
  unsigned int n_dofs_total; // the total number of degrees of freedom
  vector<CPU_RAM_Block> blocks; // a vector holding all the CPU_RAM_Blocks to be deployed to the two RAM banks

  RAM(size_t cpu_n_bytes, size_t gpu_n_bytes,
      Bat *bat, unsigned int n_bins) {
      
    this->n_dihedrals = bat->get_n_dihedrals(); // get the number of dihedrals from the .(g)bat file
    this->n_dofs_total = 3 * (n_dihedrals + 1); // calculate the total number number of degrees of freedom from the number of dihedrals 

    gpu_ram_layout = new GPU_RAM_Layout(bat->get_n_frames_padded(4), n_bins,
                                        gpu_n_bytes, n_dofs_total); // create the RAM layout on the GPU 
    cpu_ram_layout =
        new CPU_RAM_Layout(bat->get_n_frames_padded(4), cpu_n_bytes, gpu_ram_layout->dofs_per_block, n_dihedrals); // create the RAM layout on the CPU using gpu_ram_layout->dofs_per_block
    for (unsigned int i = 0; i < n_dofs_total;
         i += cpu_ram_layout->dofs_per_block) { // create CPU_RAM_Blocks according to cpu_ram_layout->dofs_per_block
      unsigned int end_g = i + cpu_ram_layout->dofs_per_block - 1;
      if (end_g > n_dofs_total - 1) // the last block can contain less dofs than the others
        end_g = n_dofs_total - 1;
      blocks.push_back(*new CPU_RAM_Block( // save the created blocks in the blocks vector
          i, end_g, gpu_ram_layout->dofs_per_block, bat, cpu_ram_layout->minima,
          cpu_ram_layout->maxima, n_bins, cpu_ram_layout->result_entropy1D));
    }
  }
};

// This class orchestrates the whole entropy calculation
class PARENT_GPU {
public:
  RAM *ram; // an object containing the CPU as well as the GPU RAM layout
  unsigned int n_bins;
  unsigned int n_frames;
  unsigned int n_frames_padded;
  Entropy_Matrix *ent_mat;
  unsigned int n_dihedrals;
  Bat *bat;
    cudaStream_t streams[N_STREAMS];

  PARENT_GPU(size_t cpu_n_bytes,
             size_t gpu_n_bytes, char const *bat_str,
             unsigned int n_bins) {
    this->n_bins = n_bins;

    bat = new Bat(bat_str);

    this->n_frames = bat->get_n_frames();
    this->n_frames_padded = bat->get_n_frames_padded(4);
    this->n_dihedrals = bat->get_n_dihedrals();

    ram = new RAM(cpu_n_bytes, gpu_n_bytes, bat, n_bins);
    ent_mat = new Entropy_Matrix(bat_str, ram->cpu_ram_layout->result_entropy,
                                 n_bins);
  }

  void calculate_entropy() {
    unsigned int skip = 0;
    unsigned int skip_next = 0;
  
    for (unsigned int i = 0; i < N_STREAMS; i++) {
        cudaStreamCreate(streams + i);
    }
    gpuErrchk(cudaPeekAtLastError());

    for (unsigned int i = 0; i < ram->blocks.size() - 1; i++) { //using this more complicated scheme, every but the first hard-dsik read is converted into a pair calculation
        skip = skip_next;
        if((skip!=0) && (i==skip)){
            cout << endl
           << "Deploying Block " << i + 2 << " (dofs "
           << ram->blocks[i].dof_id_start_g << " to " << ram->blocks[i].dof_id_end_g
           << ") to RAM bank 1."
           << endl;
            ram->blocks[i+1].deploy(ram->cpu_ram_layout->dof_block_1);
            
        cout<<"Calculating Block "<<skip + 1<<"."<<endl;
            calculate_block_cpu(&(ram->blocks[skip]));
        
        cout<<"Calculating Block pair "<<skip + 1<<" and "<<i+2<<"."<<endl;
        calculate_block_pair_cpu(&(ram->blocks[skip]), &(ram->blocks[i+1]));
            break;
        }
        else{
            cout << endl
           << "Deploying Block " << i + 1 << " (dofs "
           << ram->blocks[i].dof_id_start_g << " to " << ram->blocks[i].dof_id_end_g
           << ") to RAM bank 1."
           << endl;
            ram->blocks[i].deploy(ram->cpu_ram_layout->dof_block_1);
            cout<<"Calculating Block "<<i+1<<"."<<endl;
            calculate_block_cpu(&(ram->blocks[i]));
        }
        
        if(skip != 0){
            calculate_block_pair_cpu(&(ram->blocks[i]), &(ram->blocks[skip]));
             cout<<"Calculating Block pair "<<i+1<<" and "<<skip+1<<"."<<endl;
        }
      for (unsigned int j = i + 1; j < ram->blocks.size(); j++) {
        if(skip != j){
            skip_next = j;
            cout << endl
                 << "Deploying Block " << j + 1 << " (dofs "
                 << ram->blocks[j].dof_id_start_g << " to "
                 << ram->blocks[j].dof_id_end_g << ") to RAM bank 2." << endl;
            ram->blocks[j].deploy(ram->cpu_ram_layout->dof_block_2);
            cout<<"Calculating Block pair "<<i+1<<" and "<<j+1<<"."<<endl;
            calculate_block_pair_cpu(&(ram->blocks[i]), &(ram->blocks[j]));
            
        }

      }
    }
    cout<<"Calculating Block "<<ram->blocks.size()<<"."<<endl;
    calculate_block_cpu(&(ram->blocks[ram->blocks.size() - 1]));
    
    for (unsigned int i = 0; i < N_STREAMS; i++) {
        cudaStreamDestroy(streams[i]);
    }
    gpuErrchk(cudaPeekAtLastError());
  }

  void calculate_block_pair_cpu(CPU_RAM_Block *block1, CPU_RAM_Block *block2) {
    for (unsigned int k = 0; k < block1->blocks.size(); k++) {
        block1->blocks[k].deploy(ram->gpu_ram_layout->dof_block_1);
      for (unsigned int l = 0; l < block2->blocks.size(); l++) {
        block2->blocks[l].deploy(ram->gpu_ram_layout->dof_block_2);
        calculate_block_pair_gpu(&(block1->blocks[k]), &(block2->blocks[l]));
      }
    }
  }

  void calculate_block_cpu(CPU_RAM_Block *block) { //play the same game as for the hard-disk loads in calculate_entropy() 
  //TODO:does not yield significant performance gains, maybe revert for clearer code
    unsigned int skip = 0;
    unsigned int skip_next = 0;
    
    for (unsigned int i = 0; i < block->blocks.size() - 1; i++) {
        skip = skip_next;
        if((skip!=0) && (i==skip)){
            block->blocks[i+1].deploy(ram->gpu_ram_layout->dof_block_1);
            calculate_block_gpu(&(block->blocks[skip]));
            calculate_block_pair_gpu(&(block->blocks[skip]), &(block->blocks[i+1]));
            break;
        }
        else{
            block->blocks[i].deploy(ram->gpu_ram_layout->dof_block_1);
            calculate_block_gpu(&(block->blocks[i]));
        }
        if(skip != 0){
            block->blocks[i].deploy(ram->gpu_ram_layout->dof_block_1);
            calculate_block_pair_gpu(&(block->blocks[i]), &(block->blocks[skip]));
        }
        for (unsigned int j = i + 1; j < block->blocks.size(); j++) {
            if(skip != j){
                skip_next = j;
                block->blocks[j].deploy(ram->gpu_ram_layout->dof_block_2);
                calculate_block_pair_gpu(&(block->blocks[i]), &(block->blocks[j]));
            }
        }
    }
    calculate_block_gpu(&(block->blocks[block->blocks.size() - 1]));
  }

  unsigned char get_pair_type(unsigned char type1, unsigned char type2) {
    if ((type1 == TYPE_B) || (type2 == TYPE_B)) {
      return type1 + type2;
    } else {
      return type1 + type2 + 1;
    }
  }

  void calculate_block_pair_gpu(GPU_RAM_Block *block1, GPU_RAM_Block *block2) {
        size_t bytes_to_zero = ram->gpu_ram_layout->dofs_per_block *ram->gpu_ram_layout->dofs_per_block*
        (sizeof(PRECISION) +
         sizeof(unsigned int) *
             (n_bins * n_bins +
              1)); 
    gpuErrchk(cudaMemset(ram->gpu_ram_layout->result, 0, bytes_to_zero));


    for (unsigned int i = 0; i < block1->n_dofs; i++) {
      for (unsigned int j = 0; j < block2->n_dofs; j++) {
        
        unsigned int dof1_g = i + block1->dof_id_start_g;
        unsigned int dof2_g = j + block2->dof_id_start_g;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1_g];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2_g];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1_g] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2_g] - min2) / n_bins;
        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block2->n_dofs + j) * n_bins * n_bins;
        
        int blocks = (n_frames + (HISTOGRAM_THREADS * HISTOGRAM_THREAD_WORK_MULTIPLE - 1) ) / (HISTOGRAM_THREADS * HISTOGRAM_THREAD_WORK_MULTIPLE);
        #ifdef USE_SHARED_MEM_HISTOGRAMS
        histo2D_shared_block<<<blocks, HISTOGRAM_THREADS, n_bins * n_bins * sizeof(unsigned int), streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(
        #else
        histo2D<<<blocks, HISTOGRAM_THREADS, 0, streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(
        #endif
            block1->gpu_ram_start + i * n_frames_padded,
            block2->gpu_ram_start + j * n_frames_padded, n_frames, histogram, n_bins,
            bin_size1, bin_size2, min1, min2);
      }
    }
    //gpuErrchk(cudaPeekAtLastError());
    //gpuErrchk(cudaDeviceSynchronize());
    //for (unsigned int i = 0; i < N_STREAMS; i++) {
    //    cudaStreamDestroy(streams[i]);
    //}
    //gpuErrchk(cudaPeekAtLastError());
    

    for (unsigned int i = 0; i < block1->n_dofs; i++) {
      for (unsigned int j = 0; j < block2->n_dofs; j++) {
        unsigned int dof1_g = i + block1->dof_id_start_g;
        unsigned int dof2_g = j + block2->dof_id_start_g;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1_g];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2_g];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1_g] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2_g] - min2) / n_bins;

        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block2->n_dofs + j) * n_bins * n_bins;
        PRECISION *plnpsum =
            ram->gpu_ram_layout->result + i * block2->n_dofs + j;
        unsigned int *occupbins =
            ram->gpu_ram_layout->occupied_bins + i * block2->n_dofs + j;
        int blocks = (n_bins * n_bins + (PLNP_THREADS - 1) ) / PLNP_THREADS;

        switch (get_pair_type(block1->type, block2->type)) {
        case (TYPE_BB):
          cu_bbEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_BA):
          cu_baEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(
              histogram, n_frames, n_bins, n_bins, bin_size1, bin_size2, min1,
              min2, plnpsum, occupbins);
          break;
        case (TYPE_BD):
          cu_bdEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_AA):
          cu_aaEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_AD):
          cu_adEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_DD):
          cu_ddEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block2->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, plnpsum,
                                                  occupbins);
          break;
        }
      }
    }
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
    
    gpuErrchk(cudaMemcpy(ram->cpu_ram_layout->tmp_result_entropy,
                         ram->gpu_ram_layout->result,
                         block1->n_dofs * block2->n_dofs * sizeof(PRECISION),
                         cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(ram->cpu_ram_layout->tmp_result_occupied_bins,
                         ram->gpu_ram_layout->occupied_bins,
                         block1->n_dofs * block2->n_dofs * sizeof(unsigned int),
                         cudaMemcpyDeviceToHost));
    

    for (unsigned int i = 0; i < block1->n_dofs; i++) {
      for (unsigned int j = 0; j < block2->n_dofs; j++) {
        unsigned int dof1_gt = i + block1->dof_id_start_g -
                            get_min_id_for_type(block1->type, n_dihedrals);
        unsigned int dof2_gt = j + block2->dof_id_start_g -
                            get_min_id_for_type(block2->type, n_dihedrals);
        double entropy =
            -ram->cpu_ram_layout->tmp_result_entropy[i * block2->n_dofs + j] +
            (ram->cpu_ram_layout
                 ->tmp_result_occupied_bins[i * block2->n_dofs + j] -
             1.0) /
                (2.0 * n_frames);
        ent_mat->set2DEntropy(block1->type, block2->type, dof1_gt + 1, dof2_gt + 1,
                              entropy);
      }
    }
  }

  void calculate_block_gpu(GPU_RAM_Block *block) {
    
    size_t bytes_to_zero = ram->gpu_ram_layout->dofs_per_block;
    bytes_to_zero *=
        ram->gpu_ram_layout->dofs_per_block *
        (sizeof(PRECISION) +
         sizeof(unsigned int) *
             (n_bins * n_bins +
              1));
    gpuErrchk(cudaMemset(ram->gpu_ram_layout->result, 0, bytes_to_zero));
  

    for (unsigned int i = 0; i < block->n_dofs - 1; i++) {
      for (unsigned int j = i + 1; j < block->n_dofs; j++) {
        unsigned int dof1_g = i + block->dof_id_start_g;
        unsigned int dof2_g = j + block->dof_id_start_g;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1_g];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2_g];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1_g] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2_g] - min2) / n_bins;
        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block->n_dofs + j) * n_bins * n_bins;
      
        int blocks = (n_frames + (HISTOGRAM_THREADS * HISTOGRAM_THREAD_WORK_MULTIPLE - 1) ) / (HISTOGRAM_THREADS * HISTOGRAM_THREAD_WORK_MULTIPLE);
        #ifdef USE_SHARED_MEM_HISTOGRAMS
        histo2D_shared_block<<<blocks, HISTOGRAM_THREADS, n_bins * n_bins * sizeof(unsigned int), streams[(i * block->n_dofs + j) % N_STREAMS]>>>(
        #else
        histo2D<<<blocks, HISTOGRAM_THREADS, 0, streams[(i * block->n_dofs + j) % N_STREAMS]>>>(
        #endif
            block->gpu_ram_start + i * n_frames_padded,
            block->gpu_ram_start + j * n_frames_padded, n_frames, histogram, n_bins,
            bin_size1, bin_size2, min1, min2);
      }
    }
    //gpuErrchk(cudaPeekAtLastError());
    //gpuErrchk(cudaDeviceSynchronize());
    //for (unsigned int i = 0; i < N_STREAMS; i++) {
    //    cudaStreamDestroy(streams[i]);
    //}
    //gpuErrchk(cudaPeekAtLastError());


    for (unsigned int i = 0; i < block->n_dofs - 1; i++) {
      for (unsigned int j = i + 1; j < block->n_dofs; j++) {
        unsigned int dof1_g = i + block->dof_id_start_g;
        unsigned int dof2_g = j + block->dof_id_start_g;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1_g];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2_g];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1_g] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2_g] - min2) / n_bins;

        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block->n_dofs + j) * n_bins * n_bins;
        PRECISION *plnpsum =
            ram->gpu_ram_layout->result + i * block->n_dofs + j;
        unsigned int *occupbins =
            ram->gpu_ram_layout->occupied_bins + i * block->n_dofs + j;
        int blocks = (n_bins * n_bins + (PLNP_THREADS - 1) ) / PLNP_THREADS;

        switch (get_pair_type(block->type, block->type)) {
        case (TYPE_BB):
          cu_bbEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_BA):
          cu_baEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block->n_dofs + j) % N_STREAMS]>>>(
              histogram, n_frames, n_bins, n_bins, bin_size1, bin_size2, min1,
              min2, plnpsum, occupbins);
          break;
        case (TYPE_BD):
          cu_bdEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_AA):
          cu_aaEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_AD):
          cu_adEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_DD):
          cu_ddEnt<<<blocks, PLNP_THREADS, 0, streams[(i * block->n_dofs + j) % N_STREAMS]>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, plnpsum,
                                                  occupbins);
          break;
        }
      }
    }

    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaMemcpy(ram->cpu_ram_layout->tmp_result_entropy,
                         ram->gpu_ram_layout->result,
                         block->n_dofs * block->n_dofs * sizeof(PRECISION),
                         cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(ram->cpu_ram_layout->tmp_result_occupied_bins,
                         ram->gpu_ram_layout->occupied_bins,
                         block->n_dofs * block->n_dofs * sizeof(unsigned int),
                         cudaMemcpyDeviceToHost));
    

    for (unsigned int i = 0; i < block->n_dofs - 1; i++) {
      for (unsigned int j = i + 1; j < block->n_dofs; j++) {
        unsigned int dof1_gt = i + block->dof_id_start_g -
                            get_min_id_for_type(block->type, n_dihedrals);
        unsigned int dof2_gt = j + block->dof_id_start_g -
                            get_min_id_for_type(block->type, n_dihedrals);
        double entropy =
            -ram->cpu_ram_layout->tmp_result_entropy[i * block->n_dofs + j] +
            (ram->cpu_ram_layout
                 ->tmp_result_occupied_bins[i * block->n_dofs + j] -
             1.0) /
                (2.0 * n_frames); // includes Herzel entropy unbiasing
        ent_mat->set2DEntropy(block->type, block->type, dof1_gt + 1, dof2_gt + 1,
                              entropy);
      }
    }
  }
};

int main(int argc, char *argv[]) {

  // start the stopwatch for the execution time
  timeval tv_start, tv_end;
  gettimeofday(&tv_start, NULL);

  int deviceCount;
  gpuErrchk(cudaGetDeviceCount(&deviceCount));

  unsigned int device = 0; // TODO: implement choices for graphics card
  cout << "Found " << deviceCount
       << " CUDA device(s). Chose CUDA device number " << device << "." << endl;
  struct cudaDeviceProp prop;
  gpuErrchk(cudaGetDeviceProperties(&prop, device));
  cout << "Device name: " << prop.name << endl;
  cout << "CUDA capability: " << prop.major << "." << prop.minor << endl;
  cout << "Global memory: " << prop.totalGlobalMem / 1024 / 1024 << " MiB"
       << endl;
  cout << "Shared memory per block: " << prop.sharedMemPerBlock / 1024 << " kiB"
       << endl;
  cout << "Maximum threads per block dimension: " << prop.maxThreadsDim[0]
       << " " << prop.maxThreadsDim[1] << " " << prop.maxThreadsDim[2] << endl;
  cout << "Maximum blocks per grid dimension: " << prop.maxGridSize[0] << " "
       << prop.maxGridSize[1] << " " << prop.maxGridSize[2] << endl;
  cout << "Warp size: " << prop.warpSize << endl << endl;

  unsigned int n_bins;
  vector< vector<int> > dihedrals_top;
  vector<float> masses;
  vector<string> residues;
  vector<int> residueNumbers;
  vector<string> atomNames;
  vector<string> belongsToMolecule;

  if (argc != 11) {
    cerr << "USAGE: " << argv[0] << " -f input.[g]bat -o entropy.par -b #bins --cpu_ram #GiB --gpu_ram #GiB\n";
    exit(EXIT_FAILURE);
  }

  Arg_Parser arg_parser(argc, argv);
  if (!arg_parser.exists("-f") ||
      !arg_parser.exists("-o") ||
      !arg_parser.exists("-b")) {
    // check for correct command line options
    cerr << "USAGE: " << argv[0] << " -f input.[g]bat -o entropy.par -b #bins --cpu_ram #GiB --gpu_ram #GiB\n";
    exit(EXIT_FAILURE);
  }

  if ((strcmp(arg_parser.get_ext(arg_parser.get("-f")),
             "bat")
        && strcmp(arg_parser.get_ext(arg_parser.get("-f")),
             "gbat"))||
      strcmp(arg_parser.get_ext(arg_parser.get("-o")),
             "par")) {
    // check for the extensions of the input and output file
    cerr << "USAGE: " << argv[0] << " -f input.[g]bat -o entropy.par -b #bins --cpu_ram #GiB --gpu_ram #GiB\n";
    exit(EXIT_FAILURE);
  }
  if (sscanf(arg_parser.get("-b"), "%ud", &n_bins) != 1) {
    // read the number of bins and check for correctness
    cerr << "ERROR: Could not read number of bins from command line! Aborting"
         << endl;
    exit(EXIT_FAILURE);
  }

  stringstream cpu_ram_str(arg_parser.get("--cpu_ram"));
  double cpu_ram_provided;
  cpu_ram_str >> cpu_ram_provided;
  size_t cpu_ram_available = static_cast<size_t>(1024) * 1024 * 1024 * cpu_ram_provided;
  
  stringstream gpu_ram_str(arg_parser.get("--gpu_ram"));
  double gpu_ram_provided;
  gpu_ram_str >> gpu_ram_provided;
  size_t gpu_ram_available = static_cast<size_t>(1024) * 1024 * 1024 * gpu_ram_provided;


  PARENT_GPU parent_gpu(cpu_ram_available, gpu_ram_available,
                        arg_parser.get("-f"), n_bins);
  parent_gpu.calculate_entropy();
  cout << "Writing .par file." << endl;
  parent_gpu.ent_mat->write(arg_parser.get("-o"));
  cout << ".par file written." << endl;

  gettimeofday(&tv_end, NULL);
  cout << endl << endl;
  cout << "Total execution time: "
       << tv_end.tv_sec + 1e-6 * tv_end.tv_usec - tv_start.tv_sec -
              1e-6 * tv_start.tv_usec
       << endl;
  cout << "PROGRAM FINISHED SUCCESSFULLY." << endl << endl << endl;

  cudaDeviceReset();
  return 0;
}
