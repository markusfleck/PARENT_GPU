#!/bin/bash
rm -r bin
rm -r obj
#~ rm -r output
mkdir bin
mkdir obj
#~ mkdir output

export OMP_NUM_THREADS=24
CPU_RAM="58"
GPU_RAM="7.5"
#~ PROFILE="nsys profile -o profiles/tmp --stats=true --force-overwrite true"

#~ IN_NAME="complexes/1UGH/1UGH"
#~ OUT_NAME_BAT="complexes/output/1UGH"
#~ OUT_NAME="complexes/output/double_prec/1UGH"

IN_NAME="complexes/2KTF/UBM2_1"
OUT_NAME_BAT="complexes/output/UBM2_1"
OUT_NAME="complexes/output/double_prec/UBM2_1"

#~ IN_NAME="complexes/1JIW/1JIW"
#~ OUT_NAME_BAT="complexes/output/1JIW"
#~ OUT_NAME="complexes/output/double_prec/1JIW"

g++ -std=c++11 -c src/util/classes/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o -Wall
g++ -std=c++11 -c src/util/classes/Bat.cpp -o obj/Bat_File.o -Wall
g++ -c src/util/io/io_binary.cpp -o obj/io_binary_tmp.o -Wall
g++ -c src/util/util.cpp -o obj/util.o -Wall
ld -r obj/io_binary_tmp.o obj/Entropy_Matrix.o obj/Bat_File.o -o obj/io_binary.o
g++ -c src/util/io/io_text.cpp -o obj/io_text.o -Wall
g++ -c src/util/classes/Arg_Parser.cpp -o obj/Arg_Parser.o -Wall
ld -r obj/io_binary.o obj/io_text.o obj/Arg_Parser.o -o obj/io.o


#~ g++ -c -std=c++11 src/BAT_builder/topology.cpp -o obj/topology.o -Wall
#~ g++ -c src/BAT_builder/BAT_topology.cpp -o obj/BAT_topology.o -Wall
#~ g++ -c src/BAT_builder/BAT_trajectory.cpp -o obj/BAT_trajectory.o -Wall
#~ g++ src/BAT_builder/bat.cpp obj/io.o obj/topology.o obj/BAT_topology.o obj/BAT_trajectory.o obj/util.o -lgromacs -o bin/BAT_builder -Wall


#~ bin/BAT_builder -t ${IN_NAME}.top -x ${IN_NAME}.xtc -o ${OUT_NAME_BAT}.bat -bb "CA C N H1 O1"


#~ g++ src/BAT_builder/convert_BAT_to_GBAT.cpp obj/util.o obj/Arg_Parser.o obj/Bat_File.o -o bin/convert_BAT_to_GBAT -Wall
#~ bin/convert_BAT_to_GBAT -f ${OUT_NAME_BAT}.bat -o ${OUT_NAME_BAT}.gbat --ram $CPU_RAM


g++ --std=c++11 -O3 src/process_output/get_values_from_PAR.cpp obj/io.o obj/util.o -o bin/get_values_from_PAR -Wall

nvcc --std=c++11 -O3 -Xptxas -O3 -Xcompiler -O3,-Wall,-fopenmp  -gencode=arch=compute_61,code=\"sm_61,compute_61\" src/PARENT_GPU/PARENT_GPU.cu obj/io.o obj/util.o -o bin/PARENT_GPU \
&& $PROFILE ./bin/PARENT_GPU -f ${OUT_NAME_BAT}.gbat -o ${OUT_NAME}.par -b 50 --cpu_ram $CPU_RAM --gpu_ram $GPU_RAM\
&& bin/get_values_from_PAR -p ${OUT_NAME}.par --short 


#~ nvcc --std=c++11 -O3 -Xptxas -O3 -Xcompiler -O3,-fopenmp -gencode=arch=compute_61,code=\"sm_61,compute_61\" src/MIST_GPU/MIST_GPU.cu obj/io.o obj/util.o -o bin/MIST_GPU \
#~ && bin/MIST_GPU -f ${OUT_NAME}.par -o ${OUT_NAME}_MIST_GPU.par \
#~ && bin/get_values_from_PAR -p ${OUT_NAME}_MIST_GPU.par --short


#~ g++ --std=c++11 -O3 -fopenmp src/MIST_GPU/MIST_openMP.cpp obj/io.o obj/util.o -o bin/MIST_openMP \
#~ && bin/MIST_openMP -f ${OUT_NAME}.par -o ${OUT_NAME}_MIST_openMP.par \
#~ && bin/get_values_from_PAR -p ${OUT_NAME}_MIST_openMP.par --short


