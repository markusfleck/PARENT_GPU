# PARENT_GPU
CUDA-enabled computation of biomolecular configurational entropy from molecular dynamics trajectories
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
Change OUTNAME to "output/{name_of_your_protein}". Review the commented configuration in the 
first ten lines of run.sh. Then run

    bash run.sh
    
inside your shell. The configurational entropy output using the MIST approximation 
is contained in " output/{name_of_your_protein}.txt". 

If you run the code for testing purposes without any modifications,
the calculation should take something like 2-5 minutes. The resulting values at the end of the
files should match those in "test_system/sample_output".

If you get permission errors that means that your program has resided (or still resides) on a
filesystem which does not support Linux permissions. Move to a supported filesystem and 
unzip again or "chmod +x" the files which throw errors.

You are strongly encouraged to read the rest of this document, but at least the next section.
 


1) INSTALLATION AND TESTING

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
    make bility
    make checks
    
  Your system executes the code correctly if the following lines are the last output of ```make checks```:
     
    PARENT_GPU: pass
    MIST_GPU: pass
    MIST_openMP: pass
      
  
  After executing the "run.sh" script, the compiled executables are located in the folder
  "bin". For compilation, the "Makefile" in the top directory is processed by "run.sh", so this is where you might 
  want to start troubleshooting (or maybe fine tuning). Note that the Makefile supports compiling for a CUDA capability higher than
  6.1 by issuing ```make CUDA_ARCH={cuda_capability_without_dot}```.
