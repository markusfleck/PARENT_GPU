#!/bin/bash
rm -r bin
rm -r obj
rm -r output
mkdir bin
mkdir obj
mkdir output


g++ -std=c++11 -c src/util/classes/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o -Wall
g++ -std=c++11 -c src/util/classes/Arg_Parser.cpp -o obj/Arg_Parser.o -Wall
g++ -std=c++11 -c src/util/classes/Bat.cpp -o obj/Bat_File.o -Wall
g++ -c src/util/io/io_binary.cpp -o obj/io_binary_tmp.o -Wall
g++ -c src/util/util.cpp -o obj/util.o -Wall
ld -r obj/io_binary_tmp.o obj/Entropy_Matrix.o obj/Bat_File.o -o obj/io_binary.o
g++ -c src/util/io/io_text.cpp -o obj/io_text.o -Wall
ld -r obj/io_binary.o obj/io_text.o obj/Arg_Parser.o -o obj/io.o 
g++ --std=c++11 -O3 src/process_output/get_values_from_PAR.cpp obj/io.o obj/util.o -o bin/get_values_from_PAR

nvcc --std=c++11 -O3 -Xptxas -O3 -gencode=arch=compute_61,code=\"sm_61,compute_61\" src/PARENT_GPU/PARENT_GPU.cu obj/io.o obj/util.o -o bin/PARENT_GPU \
&& bin/PARENT_GPU -f devel/traj/UBM2_1.bat -o output/UBM2_1_MIE_GPU.par -b 50 \
&& bin/get_values_from_PAR -p output/UBM2_1_MIE_GPU.par --short 


