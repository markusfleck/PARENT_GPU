# PARENT_GPU
CUDA-enabled computation of biomolecular configurational entropy from molecular dynamics trajectories.
<br />  
<br />  
<br />  
# 0) QUICK AND DIRTY <br />  
Install the requirements. The program needs at least CUDA 9.0. Additionally, libgromacs-dev 
needs to be installed to read the trajectories. On Debian/Ubuntu/Linux Mint issue

    sudo apt install libgromacs-dev

inside a shell. Then unzip and navigate to the top folder. In the file run.sh change the parameter
IN_NAME to the path and name of your .top and .xtc files ignoring the file extensions.
(if those two files have different base names, consider changing them or using symbolic links).
Change OUT_NAME to "output/{name_of_your_protein}". Review the commented configuration in the 
first ten lines of run.sh. Then run

    bash run.sh
    
inside your shell. The configurational entropy output using the MIST approximation 
is contained in "output/{name_of_your_protein}.txt". 

If you run the code for testing purposes without any modifications,
the calculation should take something like 2-5 minutes. The resulting values at the end of the
files should match those in "test_system/sample_output".

If you get permission errors that means that your program has resided (or still resides) on a
filesystem which does not support Linux permissions. Move to a supported filesystem and 
unzip again or "chmod +x" the files which throw errors.

You are strongly encouraged to read the rest of this document, but at least the next section.
 


# 1) INSTALLATION AND TESTING

  This code uses NVIDIA CUDA, so you need to make
  sure that your system supports this. The minium requirements are CUDA 9.0
  as well as a graphics card supporting CUDA compute capability 6.1.
  
  The code was developed and tested on a GTX 1060 as well as a RTX 2060 Super,
  with CUDA 11.0 installed, Linux Mint 19.3 Tricia as the GNU/Linux operating system,
  libgromacs-dev version 2018.1-1, g++ version 7.5.0. If you run into issues, please
  consider sending me a bug report containing the according information as just given.

  Sample trajectory and topology files are shipped with this package.
  To check out if your system correctly compiles and runs all of the provided
  programs, in a bash shell type:
  
    make clean
    make
    make checks
    
  Your system executes the code correctly if the following lines are the last output of ```make checks```:
     
    PARENT_GPU: pass
    MIST_GPU: pass
    MIST_openMP: pass
      
  
  After executing the "run.sh" script, the compiled executables are located in the folder
  "bin". For compilation, the "Makefile" in the top directory is processed by "run.sh", so this is where you might 
  want to start troubleshooting (or maybe fine tuning). Note that the Makefile supports compiling for a CUDA capability higher than
  6.1 by issuing ```make CUDA_ARCH={cuda_capability_without_dot}```.
  
  
#  2) RUNNING YOUR OWN TRAJECTORIES
  
  The easiest way to do this is by just modifying the file "run.sh" in the top directory as described 
  chapter 0 as well as in run.sh by the comments:
  
  The input files (GROMACS .top and .xtc) as well as the output should be specified either 
  relative to the top directory of the package or as an absolute path.
  
  The "OUT_NAME" parameter specifies the prefix of all output files. If you want to do multiple 
  calculations using the same output folder, you should change this parameter every time. Otherwise 
  your previous results will be overwritten without warning.
  
  The BINS parameter controls how many bins are used for building the histograms which sample the 
  probability densities of the degrees of freedom. For the 2D values (on which the mutual information terms are based), 
  the square of this value will be used for building the histograms. This means if you set NBINS=50, 2500 bins will be used for every 2D histogram.
  This parameter is used exclusively during the Mutual Information Expansion (MIE) calculations, performed 
  by the program PARENT_GPU, which can be considered the core of this suite. Using 50 of them seems 
  to be a good starting point for ~10^6 frames of the MD trajectory. If you have considerably more
  (less) frames to process, you might want to tweak these parameters slightly up (down).
  In doubt, just leave this parameter as it is.
  
  The "BACKBONE_ATOMS" parameter is designed to find rigid torsion angles in your topology, 
  e. g. at a protein backbone. The names of the atoms here should match the names from 
  the GROMACS .top file. Other dihedrals are declared relative to the backbone dihedrals if they 
  share 3 atoms with them (termed phaseangles). This improves the accuracy of the
  entropy calculations, as the relative phase angles explore a much smaller range of values.
  
  In the case you want to run a considerable amount of trajectories so that compilation/check time 
  is an issue for you, you might want to uncomment the lines
  
    # make clean
    # make CUDA_ARCH=${CUDA_ARCH}
    # make checks
    
especially the line
    
    make checks
    
  in run.sh (after you compiled successfully for the first time), but keep in mind to recompile with
  	
    make clean; make CUDA_ARCH={cuda_capability_without_dot}
  
  when you change computer architecture (e. g. heterogeneous cluster).
  
#  3) EXPLANATION OF THE PROGRAMS
  
  
##  3.1) BAT_builder
    
This program converts your GROMACS topology (.top) and trajectory (.xtc) files to binary .bat
trajectory files, which are then processed by PARENT_GPU, the core program of
this suite. The main purpose of this program is to convert every frame in the .xtc file,
which is stored in Cartesian coordinates, to internal bond-angle-torsion (BAT, Z-matrix) coordinates.
Furthermore additional information is attached to the header of the resulting .bat file, namely 

	-a version number
	
	-the precision of the file (single or double precision)
	
	-the number of non-redundant torsion angles(dihedrals) in the system 
	(which relates to the number of atoms by #atoms = #torsions + 3)
	
	-the number of frames in the trajectory
	
	-a list of all non-redundant torsion angles
	(specifying their constituent atoms using the atom numbers from the .top file)
	
	-atom weights of all atoms in the system (not used yet)
	
	-atom names, residue names, residue numbers and molecule names for every atom in the system



The program is used in the following manner:

	bin/BAT_builder -t input.top -x input.xtc -o output.bat -bb "BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ..." [--single_precision]

input.top, input.xtc and output.bat are self-explanatory.

"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ..." lists the names of atoms belonging to a rigid backbone
as stated in the .top file, e. g.  "CA C N H1 O1" for a protein. Phaseangles are defined relative to a rigid dihedral. Also see section 2 for further information.

\[single_precision\] (the square brackets indicate optional), if set, writes the .bat trajectory in single precision instead of double precision,
which is discouraged, since all calculation is done in double precision anyway. Only use single precision if you are short of harddisk storage.

Additionally, the program can perform a back-conversion from .bat to .xtc, which is done by issuing the following command:

	./BAT_builder.x -b input.bat -o output.xtc


When the trajectory of a complex consisting of more than one molecule is converted, non-physical bonds (termed pseudo-bonds)
are added in order for the complex to have a connected topology. This is done automatically. Also the algorithm guarantees
that the chosen topology for every molecule in the complex is consistent with the topology which would be chosen if the molecules would be treated separately, in 
isolation.

## 3.2) convert_BAT_to_GBAT
This program was developed especially for the present GPU version of PARENT, i. e. PARENT_GPU. While the original version of PARENT was targeted for a CPU 
cluster,
PARENT_GPU is designed for workstations (and outperforms medium-sized CPU clusters with a modern consumer grade graphics card). This means that PARENT_GPU will 
in many cases reload parts of the trajectory from the harddisk if the whole trajectory does not fit into CPU RAM. GROMACS trajectories are stored in a per-frame 
fashion: for 
every time step, the coordinates of the whole molecule are stored together, followed by the next frame. This means that the trajectory of a specific coordinate
is scattered across the .xtc file. The .bat file format follows this convention. However, this format is inefficient for expansion-based entropy calculation:
To compute 2D entropies, the whole timeline of both coordinates needs to be processed, i. e. loaded from harddisk if the CPU RAM is too small to hold the whole
trajectory. The .gbat format now stores the trajectory of coordinates in a contiguous fashion, significantly improving harddisk reading times. For large projects,
this speeds up the computations significantly, given that every pair of coordinate trajectories needs to reside in CPU RAM once, requiring many harddisk fetches.
Considerably improving this I/O performance is the purpose of convert_BAT_to_GBAT. As a rule of thumb, consider using it for molecules larger than 2500 atoms
with at least a million frames (for small systems, the time this conversion takes does hardly pay off) It is used in the following manner: 

    bin/convert_BAT_to_GBAT -f input.bat -o ouput.gbat --ram #GiB
    
input.bat and output.bat are self-explanatory. #GiB is the amount of CPU RAM you can provide in GiB. Be a little conservative here: If you have 16 GB RAM 
installed, set 
10 GiB here. 


## 3.3) PARENT_GPU

This program can be considered the core of this suite. It calculates the configurational entropy according to the pairwise 
Mutual Information Expansion (MIE). A dedicated publication on this GPU version of PARENT is in preparation. Meanwhile, please read our article in the Journal of 
Chemical Theory and Computation  

"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"  
DOI: 10.1021/acs.jctc.5b01217 

and additionally

B. J. Killian, J. Y. Kravitz, and M. K. Gilson, J. Chem. Phys. 127: 024107 (2007).  
J. Numata, and E.-W. Knapp, J. Chem. Theory Comput. 8: 1235 (2012).  

PARENT_GPU takes either a .bat or a .gbat file as an input. 

The program is used in the following manner:

	bin/PARENT_GPU -f input.[g]bat -o entropy.par -b #bins --cpu_ram #GiB --gpu_ram #GiB

input.[g]bat is the result from the conversion to internal BAT coordinates done with BAT_builder or (subsequently) convert_BAT_to_GBAT.

entropy.par is the binary output file, containing all the 1D and 2D entropy terms (from which the mutual information terms can be easiliy calculated).
The header of the file includes the same information as for the .bat file as well as the numbers of bins which were used for the entropy calculation. See section 
2 for further information. The cpu_ram and gpu_ram parameters are the respective RAM amounts you want to provide. Assuming your system is not performing any other
calculations (which I would discourage), be somewhat conservative with the cpu_ram and a little conservative with the gpu_ram. E. g., if your machine has 32 GB 
RAM and your GPU has 6 GB RAM, set 25 GiB for cpu_ram and 5.25 for cpu_ram. Being to aggressive here may result in either the Linux kernel or the NVIDIA driver 
killing your program, which might be a painful experience if you run a large system where the calculation takes a long time. From my experience by now, PARENT_GPU 
is surprisingly fast even with small CPU/GPU RAM provided.

## 3.4) MIST_GPU

Although based on a different mathematical framework than MIE, the Maximum Information Spanning Tree (MIST) approximation relies on the same 
terms to be computed as for MIE. Empirically it seems to demonstrate far superior convergence properties, so from a computational perspective
one is tempted to consider MIST a refinement of MIE. We highly recommend applying MIST_GPU to your output .par file of PARENT.x (if you are 
interested in total configurational entropy values, you should consider this mandatory). Its computation time is negligible compared to PARENT_GPU. 
In addition to our article in the Journal of Chemical Theory and Computation   

"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"  
DOI: 10.1021/acs.jctc.5b01217 

please read

B. M. King, N. W. Silver, and B. Tidor. J. Phys. Chem. B 116: 2891 (2012).

The program is used in the following manner:

	bin/MIST_GPU -f input.par -o output.par

The only difference bewtween input.par and output.par is that only the significant mutual information terms are non-zero in output.par.


## 3.5) get_values_from_PAR
In order to read all entropy/mutual information terms, you need to decode the binary .par files, wether they come from the MIE (bin/PARENT_GPU) or the MIST 
(bin/MIST_GPU) calculation. For this purpose, get_values_from_PAR is provided. It lists the 1D entropies of all degrees of freedom as well as all 2D entropies and 
mutual information values of all pairs of degrees of freedom. The program is used in the following manner:

bin/get_values_from_PAR -p input.par [--short]

Its output is generally large. If you are only interested in a summary, specify the --short paramter. If you plan to e. g. use these values for machine learning,
I highly recommend to have a look at the src/util/classes/Entropy_Matrix.cpp class (in fact, I plan to wrap this class into a Python3 library).


# 4) CONTACT INFORMATION

If you have any questions, feel free to contact me: 

markus.fleck@univie.ac.at
or
fleck.markus@gmx.at

