#define PRECISION double
#define MODFITNBINS 100
#define WARPMULTIPLES 1
#define MEMORY_USAGE 0.8 // GPU
#define RAM_USAGE 0.8    // CPU
#define DEBUG false

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <vector>

#include "../util/io/io.h"
#include "../util/types.h"

#include "../util/classes/Bat.h"
#include "../util/classes/Entropy_Matrix.h"
// TODO: proper inclusion
#include "PARENT_GPU_kernels.cu"

#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t return_code, const char *file, int line) {
  if (return_code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(return_code),
            file, line);
    exit(return_code);
  }
}

using namespace std;

class GPU_RAM_Block {
public:
  unsigned char type;
  unsigned int first_dof;
  unsigned int last_dof;
  unsigned int n_dofs;
  unsigned long long int n_bytes;
  PRECISION *cpu_ram_start;
  PRECISION *gpu_ram_start;

  GPU_RAM_Block(PRECISION *cpu_ram_start, unsigned int first_dof,
                unsigned int last_dof, unsigned int n_frames,
                unsigned int n_dihedrals) {
    this->cpu_ram_start = cpu_ram_start;
    this->type = get_dof_type_from_id(first_dof, n_dihedrals);
    this->first_dof = first_dof;
    this->last_dof = last_dof;
    n_dofs = last_dof - first_dof + 1;
    n_bytes = n_dofs * n_frames * sizeof(PRECISION);
  }

  void deploy(PRECISION *gpu_ram_start) {
    gpuErrchk(cudaMemcpy(gpu_ram_start, cpu_ram_start, n_bytes,
                         cudaMemcpyHostToDevice));
    this->gpu_ram_start = gpu_ram_start;
  }
};

class CPU_RAM_Block {
public:
  unsigned int dof_id_start;
  unsigned int dof_id_end; // inclusive
  unsigned int n_dofs;
  int type_id_start[3] = {-1, -1, -1};
  int type_id_end[3] = {-1, -1, -1}; // inclusive
  unsigned int type_n_dofs[3] = {0, 0, 0};
  unsigned int gpu_ram_blocks_per_type[3];
  unsigned int n_dihedrals;
  unsigned int gpu_dofs_per_block;
  unsigned int n_frames;
  unsigned char precision_traj;
  vector<GPU_RAM_Block> blocks;
  PRECISION *block_start;
  PRECISION *minima;
  PRECISION *maxima;
  bool extrema_calculated = false;
  unsigned int n_bins;
  PRECISION *result_entropy1D;
  PRECISION *type_addr[3];
  PRECISION *bonds;
  PRECISION *angles;
  PRECISION *dihedrals;
  Bat *bat;

  CPU_RAM_Block(unsigned int dof_id_start, unsigned int dof_id_end,
                unsigned int gpu_dofs_per_block, Bat *bat, PRECISION *minima,
                PRECISION *maxima, unsigned int n_bins,
                PRECISION *result_entropy1D) {

    this->dof_id_start = dof_id_start;
    this->dof_id_end = dof_id_end;
    this->n_dofs = dof_id_end - dof_id_start + 1;
    this->n_dihedrals = bat->get_n_dihedrals();
    this->gpu_dofs_per_block = gpu_dofs_per_block;
    this->n_frames = bat->get_n_frames();
    this->precision_traj = bat->get_precision();
    this->minima = minima;
    this->maxima = maxima;
    this->n_bins = n_bins;
    this->result_entropy1D = result_entropy1D;
    this->bat = bat;

    for (unsigned short type = get_dof_type_from_id(dof_id_start, n_dihedrals);
         type <= get_dof_type_from_id(dof_id_end, n_dihedrals); type++) {
      if (dof_id_start < get_min_id_for_type(type, n_dihedrals)) {
        type_id_start[type] = get_min_id_for_type(type, n_dihedrals);
      } else {
        type_id_start[type] = dof_id_start;
      }

      if (dof_id_end > get_max_id_for_type(type, n_dihedrals)) {
        type_id_end[type] = get_max_id_for_type(type, n_dihedrals);
      } else {
        type_id_end[type] = dof_id_end;
      }
      type_n_dofs[type] = type_id_end[type] - type_id_start[type] + 1;

      gpu_ram_blocks_per_type[type] = type_n_dofs[type] / gpu_dofs_per_block;
      if (type_n_dofs[type] % gpu_dofs_per_block > 0)
        gpu_ram_blocks_per_type[type] += 1;
    }
  }

  void deploy(PRECISION *block_start) {

    this->block_start = block_start;
    type_addr[TYPE_B] = block_start;
    type_addr[TYPE_A] = block_start + type_n_dofs[TYPE_B] * n_frames;
    type_addr[TYPE_D] =
        block_start + (type_n_dofs[TYPE_B] + type_n_dofs[TYPE_A]) * n_frames;
    bonds = type_addr[TYPE_B];
    angles = type_addr[TYPE_A];
    dihedrals = type_addr[TYPE_D];

    load_dofs(block_start);
    modfit_dihedrals();
    if (!extrema_calculated) {
      calculate_extrema();
      calculate_entropy1D();
      extrema_calculated = true;
    }

    blocks.clear();
    for (unsigned char type = get_dof_type_from_id(dof_id_start, n_dihedrals);
         type <= get_dof_type_from_id(dof_id_end, n_dihedrals); type++) {
      for (unsigned int i = 0; i < gpu_ram_blocks_per_type[type]; i++) {
        unsigned int block_id_start =
            type_id_start[type] + i * gpu_dofs_per_block;
        unsigned int block_id_end =
            type_id_start[type] + (i + 1) * gpu_dofs_per_block - 1;
        if (block_id_end > type_id_end[type])
          block_id_end = type_id_end[type];
        PRECISION *cpu_ram_start =
            block_start + (block_id_start - dof_id_start) * n_frames;
        blocks.push_back(*new GPU_RAM_Block(cpu_ram_start, block_id_start,
                                            block_id_end, n_frames,
                                            n_dihedrals));
      }
    }
  }

  unsigned char load_dofs(PRECISION *dof_bank) {
    unsigned char fail = 0;
    ifstream *file = bat->get_file();

    file->seekg(bat->get_dofs_begin());
    for (unsigned int frame = 0; frame < n_frames; frame++) {
      // to read a frame from the .bat trajectory
      double ddummy[6];
      float fdummy[11];
      unsigned int a_start = get_min_id_for_type(TYPE_A, n_dihedrals);
      unsigned int d_start = get_min_id_for_type(TYPE_D, n_dihedrals);
      unsigned int b_counter = 0;
      unsigned int a_counter = a_start;
      unsigned int d_counter = d_start;
      unsigned int b_counter_local = 0;
      unsigned int a_counter_local = 0;
      unsigned int d_counter_local = 0;

      if (frame % 100000 == 0) {
        cout << "Reading frame " << frame
             << " and the following.\n"; // every 10000 frames issue an
                                         // information to stdout
        cout.flush();
      }

      file->read(
          (char *)fdummy,
          11 *
              sizeof(float)); // read time, precision and box vectors to dummies
      fail = fail | (file->rdstate() & std::ifstream::failbit);
      unsigned char inc;
      if (precision_traj == 1) {
        inc = sizeof(double); // if trajectory is in double precision
      } else {
        inc = sizeof(float);
      }
      file->read((char *)ddummy, 6 * inc); // external coordinates to dummies
      fail = fail | (file->rdstate() & std::ifstream::failbit);

      file->read((char *)ddummy,
                 inc); // read the lengths of the two bonds connecting the root
                       // atoms (internal coordinates)
      fail = fail | (file->rdstate() & std::ifstream::failbit);
      if ((b_counter >= type_id_start[TYPE_B]) &&
          (b_counter <= type_id_end[TYPE_B]))
        bonds[b_counter_local++ * n_frames + frame] = ddummy[0];
      b_counter++;

      file->read((char *)ddummy,
                 inc); // read the lengths of the two bonds connecting the root
                       // atoms (internal coordinates)
      fail = fail | (file->rdstate() & std::ifstream::failbit);
      if ((b_counter >= type_id_start[TYPE_B]) &&
          (b_counter <= type_id_end[TYPE_B]))
        bonds[b_counter_local++ * n_frames + frame] = ddummy[0];
      b_counter++;

      file->read((char *)ddummy, inc); // and the angle between the two
                                       // rootbonds (internal coordinates)
      fail = fail | (file->rdstate() & std::ifstream::failbit);
      if ((a_counter >= type_id_start[TYPE_A]) &&
          (a_counter <= type_id_end[TYPE_A]))
        angles[a_counter_local++ * n_frames + frame] = ddummy[0];
      a_counter++;

      for (int i = 0; i < n_dihedrals;
           i++) {                        // then for all dihedrals in the system
        file->read((char *)ddummy, inc); // read the bondlength between the last
                                         // two atoms in the dihedral
        fail = fail | (file->rdstate() & std::ifstream::failbit);
        if ((b_counter >= type_id_start[TYPE_B]) &&
            (b_counter <= type_id_end[TYPE_B]))
          bonds[b_counter_local++ * n_frames + frame] = ddummy[0];
        b_counter++;

        file->read((char *)ddummy, inc); // read the angle between the last
                                         // threee atoms of the dihedral#
        fail = fail | (file->rdstate() & std::ifstream::failbit);
        if ((a_counter >= type_id_start[TYPE_A]) &&
            (a_counter <= type_id_end[TYPE_A]))
          angles[a_counter_local++ * n_frames + frame] = ddummy[0];
        a_counter++;

        file->read((char *)ddummy, inc); // and the value of the dihedral itself
        fail = fail | (file->rdstate() & std::ifstream::failbit);
        if ((d_counter >= type_id_start[TYPE_D]) &&
            (d_counter <= type_id_end[TYPE_D]))
          dihedrals[d_counter_local++ * n_frames + frame] = ddummy[0];
        d_counter++;
      }
    }

    return fail; // if anything failed return a 1, otherwise a 0
  }

  void modfit_dihedrals() {
    if (type_n_dofs[TYPE_D] == 0)
      return;
    const PRECISION pi = acos(-1);

#pragma omp parallel
    {
      PRECISION modFit, binsize;
      int longestZeroStretch, longestZeroStretchPos, currentZeroStretch,
          currentZeroStretchPos;
      bool zeroExists;
      long long int histo[MODFITNBINS];

#pragma omp for
      for (int j = 0; j < type_n_dofs[TYPE_D]; j++) { // for all dihedrals
        // first build a histogram of the dihedral values over the trajectory
        for (int k = 0; k < MODFITNBINS; k++)
          histo[k] = 0;

        binsize = (2 * pi +
                   5e-9 * (sizeof(PRECISION) == sizeof(float) ? 100000 : 1)) /
                  MODFITNBINS;
        for (int i = 0; i < n_frames; i++)
          histo[int((dihedrals[j * n_frames + i]) / binsize)] += 1;

        zeroExists = false;
        for (int k = 0; k < MODFITNBINS; k++)
          zeroExists = zeroExists || (histo[k] == 0);

        if (zeroExists) { // if any of the bins of the histogram is empty find
                          // the longest consecutive stretch of  empty bins
          longestZeroStretch = 0;
          currentZeroStretch = 0;
          longestZeroStretchPos = -1;
          for (int k = 0; k < 2 * MODFITNBINS;
               k++) {                // for all bins of the histogram
            int l = k % MODFITNBINS; // taking car of zero stretches which span
                                     // the boundaries
            if ((currentZeroStretch == 0) &&
                (histo[l] == 0)) { // find and save a beginning zero stretch
              currentZeroStretch = 1;
              currentZeroStretchPos = k;
            }
            if ((currentZeroStretch > 0) && (histo[l] == 0)) {
              currentZeroStretch += 1;
            }
            if ((currentZeroStretch > 0) &&
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
              0; // misuse the zeroStretch variables for determining the minimum
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
        for (int k = 0; k < n_frames; k++) {
          dihedrals[j * n_frames + k] =
              dihedrals[j * n_frames + k] + modFit -
              2 * pi *
                  int((dihedrals[j * n_frames + k] + modFit) /
                      (2 * pi)); // and apply it taking care of circularity
        }
      }
    }
  }

  void calculate_extrema() {
#pragma omp parallel
    {
      PRECISION tmpMin, tmpMax;

#pragma omp for
      for (int j = 0; j < n_dofs; j++) { // for all dofs
        tmpMax = block_start[j * n_frames];
        tmpMin = block_start[j * n_frames];
        for (int i = 1; i < n_frames; i++) { // and all frames
          if (block_start[j * n_frames + i] > tmpMax) {
            tmpMax = block_start[j * n_frames + i];
          }
          if (block_start[j * n_frames + i] < tmpMin) {
            tmpMin = block_start[j * n_frames +
                                 i]; // find the maximum and minmum values
          }
        }
        if ((tmpMin < 0.0) || (tmpMax < 0.0)) {
          cerr << "ERROR: Degree of freedom " << dof_id_start + j
               << " is smaller than 0.0" << endl;
          exit(EXIT_FAILURE);
        }
        tmpMin -= 5e-9 * (sizeof(PRECISION) == sizeof(float)
                              ? 1e5
                              : 1); // and increase the boundaries a tiny bit
        tmpMax += 5e-9 * (sizeof(PRECISION) == sizeof(float) ? 1e5 : 1);
        if (tmpMin < 0.0) {
          tmpMin = 0.0;
        }
        if ((tmpMax - tmpMin) < 1.0e-4) {
          tmpMax += 0.05;
        }
        minima[dof_id_start + j] = tmpMin; // and store the calculated values
        maxima[dof_id_start + j] = tmpMax;
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
      for (int j = dof_id_start; j <= dof_id_end;
           j++) { // for all dofs (using threads)
        for (int k = 0; k < n_bins; k++) {
          histo[k] = 0; // initialize a histogram with zeros
        }
        binsize = (maxima[j] - minima[j]) /
                  n_bins; // and calculate the size of the bins
        for (int i = 0; i < n_frames;
             i++) { // and fill the histogram using all frames of the trajectory
          histo[int(
              (block_start[(j - dof_id_start) * n_frames + i] - minima[j]) /
              binsize)] += 1;
        }
        occupbins = 0; // then use the histogram to calculate the (discretized)
                       // entropy, taking care of the Jacobian
        plnpsum = 0;
        binval = minima[j] + (binsize / 2.0);
        for (int k = 0; k < n_bins; k++) {
          switch (get_dof_type_from_id(j, n_dihedrals)) {
          case TYPE_B:
            Jac = binval * binval;
            break;
          case TYPE_A:
            Jac = sin(binval);
            break;
          case TYPE_D:
            Jac = 1;
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

class GPU_RAM_Layout {
public:
  unsigned int dofs_per_block;
  PRECISION *dof_block_1;
  PRECISION *dof_block_2;
  PRECISION *result;
  unsigned int *occupied_bins;
  unsigned int *histograms;

  GPU_RAM_Layout(unsigned int n_frames, unsigned int n_bins,
                 unsigned long long int gpu_n_bytes, char *gpu_ram_start) {
    // calculate the maximum number of dofs (for one of the two dof_blocks) so
    // that everything still fits into GPU RAM
    double a = 2 * n_frames * sizeof(PRECISION);
    double b = sizeof(PRECISION) + sizeof(unsigned int) * (n_bins * n_bins + 1);
    this->dofs_per_block =
        (unsigned int)((-a / 2 + sqrt(a * a / 4 + gpu_n_bytes * b)) / b);

    // set the pointers for partitioning the GPU RAM according to the calculated
    // dofs_per_block
    dof_block_1 = (PRECISION *)gpu_ram_start;
    dof_block_2 = dof_block_1 + dofs_per_block * n_frames;
    result = dof_block_2 + dofs_per_block * n_frames;
    occupied_bins = (unsigned int *)(result + dofs_per_block * dofs_per_block);
    histograms = occupied_bins + dofs_per_block * dofs_per_block;
  }
};

class CPU_RAM_Layout {
public:
  unsigned int dofs_per_block;
  PRECISION *dof_block_1;
  PRECISION *dof_block_2;
  PRECISION *result_entropy;
  PRECISION *result_entropy1D;
  PRECISION *result_entropy1D_b;
  PRECISION *result_entropy1D_a;
  PRECISION *result_entropy1D_d;
  PRECISION *result_entropy2D;
  PRECISION *result_entropy2D_bb;
  PRECISION *result_entropy2D_ba;
  PRECISION *result_entropy2D_bd;
  PRECISION *result_entropy2D_aa;
  PRECISION *result_entropy2D_ad;
  PRECISION *result_entropy2D_dd;
  PRECISION *extrema;
  PRECISION *minima;
  PRECISION *maxima;
  PRECISION *tmp_result_entropy;
  unsigned int *tmp_result_occupied_bins;
  double *tmp_read;

  CPU_RAM_Layout(unsigned int n_frames, unsigned long long int cpu_n_bytes,
                 char *cpu_ram_start, unsigned int gpu_dofs_per_block,
                 unsigned int n_dihedrals) {
    unsigned int n_dofs_total = 3 * (n_dihedrals + 1);
    unsigned int n_bonds = n_dihedrals + 2;
    unsigned int n_angles = n_dihedrals + 1;
    // calculate the maximum number of dofs (for one of the two dof_blocks) so
    // that everything still fits into CPU RAM
    dofs_per_block =
        (unsigned int)((cpu_n_bytes -
                        n_dofs_total * ((n_dofs_total + 3) * sizeof(PRECISION) +
                                        sizeof(double)) +
                        (2 * gpu_dofs_per_block - 1) * gpu_dofs_per_block *
                            (sizeof(PRECISION) + sizeof(unsigned int))) /
                       (2 * n_frames * sizeof(PRECISION)));
    if (dofs_per_block < gpu_dofs_per_block) {
      cerr << "WARNING: You probably have a GPU with a lot of RAM but your CPU "
              "RAM is rather small. ";
      cerr << "I recommend to get more CPU RAM, as this should significantly "
              "enhance performance."
           << endl;
    }
    // if all dofs fit into RAM, still set up two blocks to be consistent with
    // the algorithm
    if (2 * dofs_per_block >= n_dofs_total) {
      dofs_per_block = n_dofs_total / 2;
      dofs_per_block += n_dofs_total % 2;
    }

    // set the pointers for partitioning the CPU RAM according to the calculated
    // dofs_per_block
    dof_block_1 = (PRECISION *)cpu_ram_start;
    dof_block_2 = dof_block_1 + dofs_per_block * n_frames;

    result_entropy = dof_block_2 + dofs_per_block * n_frames;
    result_entropy1D = result_entropy;
    result_entropy1D_b = result_entropy1D;
    result_entropy1D_a = result_entropy1D_b + n_bonds;
    result_entropy1D_d = result_entropy1D_a + n_angles;

    result_entropy2D = result_entropy1D + n_dofs_total;
    result_entropy2D_bb = result_entropy2D;
    result_entropy2D_ba = result_entropy2D_bb + n_bonds * (n_bonds - 1) / 2;
    result_entropy2D_bd = result_entropy2D_ba + n_bonds * n_angles;
    result_entropy2D_aa = result_entropy2D_bd + n_bonds * n_dihedrals;
    result_entropy2D_ad = result_entropy2D_aa + n_angles * (n_angles - 1) / 2;
    result_entropy2D_dd = result_entropy2D_ad + n_angles * n_dihedrals;

    extrema = result_entropy2D + n_dofs_total * (n_dofs_total - 1) / 2;
    minima = extrema;
    maxima = extrema + n_dofs_total;
    tmp_result_entropy = maxima + n_dofs_total;
    tmp_result_occupied_bins =
        (unsigned int *)(tmp_result_entropy +
                         (2 * gpu_dofs_per_block - 1) * gpu_dofs_per_block);
    tmp_read = (double *)(tmp_result_occupied_bins +
                          (2 * gpu_dofs_per_block - 1) * gpu_dofs_per_block);
  }
};

class RAM {
public:
  char *cpu_ram_start;
  char *cpu_ram_end;
  unsigned long long int cpu_n_bytes;
  char *gpu_ram_start;
  char *gpu_ram_end;
  unsigned long long int gpu_n_bytes;
  GPU_RAM_Layout *gpu_ram_layout;
  CPU_RAM_Layout *cpu_ram_layout;
  unsigned int n_dihedrals;
  unsigned int n_dofs_total;
  vector<CPU_RAM_Block> blocks;

  RAM(unsigned long long int cpu_n_bytes, unsigned long long int gpu_n_bytes,
      Bat *bat, unsigned int n_bins) {
    cpu_ram_start = new char[cpu_n_bytes];
    cpu_ram_end = cpu_ram_start + cpu_n_bytes - 1;
    this->cpu_n_bytes = cpu_n_bytes;
    gpuErrchk(cudaMalloc((void **)&gpu_ram_start, gpu_n_bytes));
    gpu_ram_end = gpu_ram_start + gpu_n_bytes - 1;
    this->gpu_n_bytes = gpu_n_bytes;
    gpu_ram_layout = new GPU_RAM_Layout(bat->get_n_frames(), n_bins,
                                        gpu_n_bytes, gpu_ram_start);
    this->n_dihedrals = bat->get_n_dihedrals();
    this->n_dofs_total = 3 * (n_dihedrals + 1);
    cpu_ram_layout =
        new CPU_RAM_Layout(bat->get_n_frames(), cpu_n_bytes, cpu_ram_start,
                           gpu_ram_layout->dofs_per_block, n_dihedrals);
    for (unsigned int i = 0; i < n_dofs_total;
         i += cpu_ram_layout->dofs_per_block) {
      unsigned int end = i + cpu_ram_layout->dofs_per_block - 1;
      if (end > n_dofs_total - 1)
        end = n_dofs_total - 1;
      blocks.push_back(*new CPU_RAM_Block(
          i, end, gpu_ram_layout->dofs_per_block, bat, cpu_ram_layout->minima,
          cpu_ram_layout->maxima, n_bins, cpu_ram_layout->result_entropy1D));
    }
  }
};

class PARENT_GPU {
public:
  RAM *ram;
  int threads_per_block;
  unsigned int n_bins;
  unsigned int n_frames;
  Entropy_Matrix *ent_mat;
  unsigned int n_dihedrals;
  Bat *bat;

  PARENT_GPU(unsigned long long int cpu_n_bytes,
             unsigned long long int gpu_n_bytes, char const *bat_str,
             unsigned int n_bins, int threads_per_block) {
    this->threads_per_block = threads_per_block;
    this->n_bins = n_bins;

    bat = new Bat(bat_str);

    this->n_frames = bat->get_n_frames();
    this->n_dihedrals = bat->get_n_dihedrals();

    ram = new RAM(cpu_n_bytes, gpu_n_bytes, bat, n_bins);
    ent_mat = new Entropy_Matrix(bat_str, ram->cpu_ram_layout->result_entropy,
                                 n_bins);
  }

  void calculate_entropy() {

    for (unsigned int i = 0; i < ram->blocks.size() - 1; i++) {
      cout << endl
           << "Deploying Block " << i + 1 << " (dofs "
           << ram->blocks[i].dof_id_start << " to " << ram->blocks[i].dof_id_end
           << ") to RAM bank 1."
           << endl; // TODO: think of more efficient deploying
      ram->blocks[i].deploy(ram->cpu_ram_layout->dof_block_1);
      calculate_block_cpu(&(ram->blocks[i]));
      for (unsigned int j = i + 1; j < ram->blocks.size(); j++) {

        cout << endl
             << "Deploying Block " << j + 1 << " (dofs "
             << ram->blocks[j].dof_id_start << " to "
             << ram->blocks[j].dof_id_end << ") to RAM bank 2." << endl;
        ram->blocks[j].deploy(ram->cpu_ram_layout->dof_block_2);

        calculate_block_pair_cpu(&(ram->blocks[i]), &(ram->blocks[j]));
      }
    }

    calculate_block_cpu(&(ram->blocks[ram->blocks.size() - 1]));
  }

  void calculate_block_pair_cpu(CPU_RAM_Block *block1, CPU_RAM_Block *block2) {
    for (unsigned int k = 0; k < block1->blocks.size(); k++) {
      for (unsigned int l = 0; l < block2->blocks.size(); l++) {
        calculate_block_pair_gpu(&(block1->blocks[k]), &(block2->blocks[l]));
      }
    }
  }

  void calculate_block_cpu(CPU_RAM_Block *block) {
    for (unsigned int k = 0; k < block->blocks.size() - 1; k++) {
      calculate_block_gpu(
          &(block->blocks[k])); // TODO: think of more efficient deploying
      for (unsigned int l = k + 1; l < block->blocks.size(); l++) {
        calculate_block_pair_gpu(&(block->blocks[k]), &(block->blocks[l]));
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
    // TODO: deploy both blocks at once
    block1->deploy(ram->gpu_ram_layout->dof_block_1);
    block2->deploy(ram->gpu_ram_layout->dof_block_2);

    size_t bytes_to_zero = ram->gpu_ram_layout->dofs_per_block;
    bytes_to_zero *=
        ram->gpu_ram_layout->dofs_per_block *
        (sizeof(PRECISION) +
         sizeof(unsigned int) *
             (n_bins * n_bins +
              1)); // TODO: zero only what is needed(block1->n_dofs instead of
                   // ram->gpu_ram_layout->dofs_per_block)
    gpuErrchk(cudaMemset(ram->gpu_ram_layout->result, 0, bytes_to_zero));

    for (unsigned int i = 0; i < block1->n_dofs; i++) {
      for (unsigned int j = 0; j < block2->n_dofs; j++) {
        unsigned int dof1 = i + block1->first_dof;
        unsigned int dof2 = j + block2->first_dof;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2] - min2) / n_bins;
        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block2->n_dofs + j) * n_bins * n_bins;
        int blocks = n_frames / threads_per_block;
        if (n_frames % threads_per_block > 0)
          blocks++;
        histo2D<<<blocks, threads_per_block>>>(
            block1->gpu_ram_start + i * n_frames,
            block2->gpu_ram_start + j * n_frames, n_frames, histogram, n_bins,
            bin_size1, bin_size2, min1, min2);
      }
    }
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    for (unsigned int i = 0; i < block1->n_dofs; i++) {
      for (unsigned int j = 0; j < block2->n_dofs; j++) {
        unsigned int dof1 = i + block1->first_dof;
        unsigned int dof2 = j + block2->first_dof;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2] - min2) / n_bins;

        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block2->n_dofs + j) * n_bins * n_bins;
        PRECISION *plnpsum =
            ram->gpu_ram_layout->result + i * block2->n_dofs + j;
        unsigned int *occupbins =
            ram->gpu_ram_layout->occupied_bins + i * block2->n_dofs + j;
        int blocks = n_bins * n_bins / threads_per_block;
        if (n_bins * n_bins % threads_per_block > 0)
          blocks++;

        switch (get_pair_type(block1->type, block2->type)) {
        case (TYPE_BB):
          cu_bbEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_BA):
          cu_baEnt<<<blocks, threads_per_block>>>(
              histogram, n_frames, n_bins, n_bins, bin_size1, bin_size2, min1,
              min2, plnpsum, occupbins);
          break;
        case (TYPE_BD):
          cu_bdEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_AA):
          cu_aaEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_AD):
          cu_adEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_DD):
          cu_ddEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
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
                         block1->n_dofs * block2->n_dofs * sizeof(PRECISION),
                         cudaMemcpyDeviceToHost));

    for (unsigned int i = 0; i < block1->n_dofs; i++) {
      for (unsigned int j = 0; j < block2->n_dofs; j++) {
        unsigned int dof1 = i + block1->first_dof -
                            get_min_id_for_type(block1->type, n_dihedrals);
        unsigned int dof2 = j + block2->first_dof -
                            get_min_id_for_type(block2->type, n_dihedrals);
        double entropy =
            -ram->cpu_ram_layout->tmp_result_entropy[i * block2->n_dofs + j] +
            (ram->cpu_ram_layout
                 ->tmp_result_occupied_bins[i * block2->n_dofs + j] -
             1.0) /
                (2.0 * n_frames);
        ent_mat->set2DEntropy(block1->type, block2->type, dof1 + 1, dof2 + 1,
                              entropy);
      }
    }
  }

  void calculate_block_gpu(GPU_RAM_Block *block) {
    // TODO: Use already deployed blocks
    block->deploy(ram->gpu_ram_layout->dof_block_1);

    size_t bytes_to_zero = ram->gpu_ram_layout->dofs_per_block;
    bytes_to_zero *=
        ram->gpu_ram_layout->dofs_per_block *
        (sizeof(PRECISION) +
         sizeof(unsigned int) *
             (n_bins * n_bins +
              1)); // TODO: zero only what is needed(block1->n_dofs instead of
                   // ram->gpu_ram_layout->dofs_per_block)
    gpuErrchk(cudaMemset(ram->gpu_ram_layout->result, 0, bytes_to_zero));

    for (unsigned int i = 0; i < block->n_dofs - 1; i++) {
      for (unsigned int j = i + 1; j < block->n_dofs; j++) {
        unsigned int dof1 = i + block->first_dof;
        unsigned int dof2 = j + block->first_dof;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2] - min2) / n_bins;
        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block->n_dofs + j) * n_bins * n_bins;
        int blocks = n_frames / threads_per_block;
        if (n_frames % threads_per_block > 0)
          blocks++;
        histo2D<<<blocks, threads_per_block>>>(
            block->gpu_ram_start + i * n_frames,
            block->gpu_ram_start + j * n_frames, n_frames, histogram, n_bins,
            bin_size1, bin_size2, min1, min2);
      }
    }
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    for (unsigned int i = 0; i < block->n_dofs - 1; i++) {
      for (unsigned int j = i + 1; j < block->n_dofs; j++) {
        unsigned int dof1 = i + block->first_dof;
        unsigned int dof2 = j + block->first_dof;
        PRECISION min1 = ram->cpu_ram_layout->minima[dof1];
        PRECISION min2 = ram->cpu_ram_layout->minima[dof2];

        PRECISION bin_size1 =
            (ram->cpu_ram_layout->maxima[dof1] - min1) / n_bins;
        PRECISION bin_size2 =
            (ram->cpu_ram_layout->maxima[dof2] - min2) / n_bins;

        unsigned int *histogram = ram->gpu_ram_layout->histograms +
                                  (i * block->n_dofs + j) * n_bins * n_bins;
        PRECISION *plnpsum =
            ram->gpu_ram_layout->result + i * block->n_dofs + j;
        unsigned int *occupbins =
            ram->gpu_ram_layout->occupied_bins + i * block->n_dofs + j;
        int blocks = n_bins * n_bins / threads_per_block;
        if (n_bins * n_bins % threads_per_block > 0)
          blocks++;

        switch (get_pair_type(block->type, block->type)) {
        case (TYPE_BB):
          cu_bbEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_BA):
          cu_baEnt<<<blocks, threads_per_block>>>(
              histogram, n_frames, n_bins, n_bins, bin_size1, bin_size2, min1,
              min2, plnpsum, occupbins);
          break;
        case (TYPE_BD):
          cu_bdEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_AA):
          cu_aaEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  bin_size1, bin_size2, min1,
                                                  min2, plnpsum, occupbins);
          break;
        case (TYPE_AD):
          cu_adEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
                                                  n_bins, bin_size1, bin_size2,
                                                  min1, plnpsum, occupbins);
          break;
        case (TYPE_DD):
          cu_ddEnt<<<blocks, threads_per_block>>>(histogram, n_frames, n_bins,
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
                         block->n_dofs * block->n_dofs * sizeof(PRECISION),
                         cudaMemcpyDeviceToHost)); // TODO: use GPU fwd block
    gpuErrchk(cudaMemcpy(ram->cpu_ram_layout->tmp_result_occupied_bins,
                         ram->gpu_ram_layout->occupied_bins,
                         block->n_dofs * block->n_dofs * sizeof(PRECISION),
                         cudaMemcpyDeviceToHost));

    for (unsigned int i = 0; i < block->n_dofs - 1; i++) {
      for (unsigned int j = i + 1; j < block->n_dofs; j++) {
        unsigned int dof1 = i + block->first_dof -
                            get_min_id_for_type(block->type, n_dihedrals);
        unsigned int dof2 = j + block->first_dof -
                            get_min_id_for_type(block->type, n_dihedrals);
        double entropy =
            -ram->cpu_ram_layout->tmp_result_entropy[i * block->n_dofs + j] +
            (ram->cpu_ram_layout
                 ->tmp_result_occupied_bins[i * block->n_dofs + j] -
             1.0) /
                (2.0 * n_frames); // includes Herzel entropy unbiasing
        ent_mat->set2DEntropy(block->type, block->type, dof1 + 1, dof2 + 1,
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

  int threads_per_block = prop.warpSize * WARPMULTIPLES;

  unsigned int n_bins;
  vector<vector<int>> dihedrals_top;
  vector<float> masses;
  vector<string> residues;
  vector<int> residueNumbers;
  vector<string> atomNames;
  vector<string> belongsToMolecule;

  if (argc != 7) {
    cerr << "USAGE: " << argv[0] << " -f input.bat -o entropy.par -b #bins\n";
    exit(EXIT_FAILURE);
  }
  bool fail = false;
  fail |= !cmdOptionExists(argv, argv + argc, "-f");
  fail |= !cmdOptionExists(argv, argv + argc, "-o");
  fail |= !cmdOptionExists(argv, argv + argc, "-b");

  if (fail) {
    // check for correct command line options
    cerr << "USAGE: " << argv[0] << " -f input.bat -o entropy.par -b #bins\n";
    exit(EXIT_FAILURE);
  }

  string tmp1(getCmdOption(argv, argv + argc,
                           "-f")); // first argument is the .bat trajectory file
  string tmp2(getCmdOption(argv, argv + argc,
                           "-o")); // second argument is the .par output file
  char *ptr, *type1, *type2;
  char delimiter[] = ".";

  ptr = strtok((char *)tmp1.c_str(), delimiter);
  while (ptr != NULL) {
    type1 = ptr;
    ptr = strtok(NULL, delimiter);
  }

  ptr = strtok((char *)tmp2.c_str(), delimiter);
  while (ptr != NULL) {
    type2 = ptr;
    ptr = strtok(NULL, delimiter);
  }
  if ((strcmp(type1, "bat")) || (strcmp(type2, "par"))) {
    // check for the extensions of the input and output file
    cerr << "USAGE: " << argv[0] << " -f input.bat -o entropy.par -b #bins\n";
    exit(EXIT_FAILURE);
  }
  if (sscanf(getCmdOption(argv, argv + argc, "-b"), "%ud", &n_bins) != 1) {
    // read the number of bins and check for correctness
    cerr << "ERROR: Could not read number of bins from command line! Aborting"
         << endl;
    exit(EXIT_FAILURE);
  }

  unsigned long long int cpu_ram_available =
      static_cast<unsigned long long int>(1024) * 1024 * 1024 * 4;
  unsigned long long int gpu_ram_available =
      static_cast<unsigned long long int>(1024) * 1024 * 1024 * 1;

  PARENT_GPU parent_gpu(cpu_ram_available, gpu_ram_available,
                        getCmdOption(argv, argv + argc, "-f"), n_bins,
                        threads_per_block);
  parent_gpu.calculate_entropy();
  cout << "Writing .par file." << endl;
  parent_gpu.ent_mat->write(getCmdOption(argv, argv + argc, "-o"));
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
