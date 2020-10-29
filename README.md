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
