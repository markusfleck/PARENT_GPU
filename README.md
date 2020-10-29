# PARENT_GPU
CUDA-enabled computation of biomolecular configurational entropy from molecular dynamics trajectories
<br />  
<br />  
<br />  
# 0) QUICK AND DIRTY <br />  
Install the requirements. The program needs at least CUDA 9.0. Additionally, libgromacs-dev 
needs to be installes to read the trajectories. On Debian/Ubuntu/Linux Mint issue

    sudo apt install libgromacs-dev

inside a shell. Then unzip and navigate to the top folder. In the file run.sh change the parameter
IN_NAME to the path and name of your .top and .xtc files without those file extensions
(if those two files have different base names, consider changing them or using symbolic links).
Change OUTNAME to "output/{name_of_your_protein}". Review the commented configuration in the 
first ten lines of run.sh. Then run

    bash run.sh
    
inside your shell.
 


