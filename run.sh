#!/bin/bash

# Example script to operate the PAREN_GPU suite end-to-end
# Copyright (C) 2020  Markus Fleck
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as 
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

CPU_RAM="8.0" # The CPU RAM you would like to provide in GiB. Be somewhat conservative here. E. g., if your machine has 32 GB RAM, set 25  
GPU_RAM="2.0" # The GPU RAM you would like to provide in GiB. Be a litte conservative here. E. g., if your GPU has 6 GB RAM, set 5.25
CUDA_ARCH=61 # The CUDA capability of your GPU, normally the higher the better. In doubt, leave it as it is: 61 is the minimum required
export CUDA_VISIBLE_DEVICES=0 # If you have more than one GPUs, you can choose which one to use here
export OMP_NUM_THREADS=4 # Set the number of threads of your CPU.

IN_NAME="test_system/UBQ_UBM2" # The name and location of your trajectory (.xtc) and topology (.top) files
OUT_NAME="output/UBQ_UBM" # The name and location of the output
NBINS=50 # The number of bins used for the discretization of the trajectory. The 2D entropy values use the square of this nuber as the number of bins
BACKBONE_ATOMS="CA C N H1 O1" # The names of the backbone atoms. Adjusting them to the nomenclature of your .top file improves numerical accuracy by using phaseangles 




# make clean # remember to recompile if you change machines (e. g. when using cloud computing services)
make CUDA_ARCH=${CUDA_ARCH}
make checks; echo -e "\n\n\n" # run checks to make sure your system produces correct results. Consider commenting this line if your system configuration has not changed 

# rm -r output # remember that files will be overwritten without a warning
mkdir output

bin/BAT_builder -t ${IN_NAME}.top -x ${IN_NAME}.xtc -o ${OUT_NAME}.bat -bb "${BACKBONE_ATOMS}" # convert the trajectory from GROMACS .xtc to PARENT/PARENT_GPU .bat format. Remember that .bat files are generally large, so make sure you provide enough harddisk space  
bin/convert_BAT_to_GBAT -f ${OUT_NAME}.bat -o ${OUT_NAME}.gbat --ram $CPU_RAM # optionally, convert .bat to .gbat, tremendously enhancing harddisk reading times. Useful for large molecules/trajectories. 


./bin/PARENT_GPU -f ${OUT_NAME}.gbat -o ${OUT_NAME}.par -b $NBINS --cpu_ram $CPU_RAM --gpu_ram $GPU_RAM && bin/get_values_from_PAR -p ${OUT_NAME}.par --short | tee ${OUT_NAME}_MIE.txt # run the MIE calculation, the heart of PARENT_GPU 
echo -e "\n\n\n"

bin/MIST_GPU -f ${OUT_NAME}.par -o ${OUT_NAME}_MIST_GPU.par && bin/get_values_from_PAR -p ${OUT_NAME}_MIST_GPU.par --short | tee ${OUT_NAME}_MIST.txt # calculate the MIST approximation (in principle optional, but numerically mandatory)
# echo -e "\n\n\n"; bin/MIST_openMP -f ${OUT_NAME}.par -o ${OUT_NAME}_MIST_openMP.par && bin/get_values_from_PAR -p ${OUT_NAME}_MIST_openMP.par --short # no need to run this line unless for some exotic reason you don't want to use your GPU to calculate the MIST approximation as done just above


