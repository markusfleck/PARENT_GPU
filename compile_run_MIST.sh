#!/bin/bash
rm -r bin
rm -r obj
#~ rm -r output
mkdir bin
mkdir obj
#~ mkdir output

g++ -std=c++11 -c src/util/classes/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o -Wall
g++ -std=c++11 -c src/util/classes/Bat.cpp -o obj/Bat_File.o -Wall
g++ -c src/util/io/io_binary.cpp -o obj/io_binary_tmp.o -Wall
g++ -c src/util/util.cpp -o obj/util.o -Wall
ld -r obj/io_binary_tmp.o obj/Entropy_Matrix.o obj/Bat_File.o -o obj/io_binary.o
g++ -c src/util/io/io_text.cpp -o obj/io_text.o -Wall
g++ -c src/util/classes/Arg_Parser.cpp -o obj/Arg_Parser.o -Wall
ld -r obj/io_binary.o obj/io_text.o obj/Arg_Parser.o -o obj/io.o

#~ gcc -c src/util/io/xdrfile/xdrfile.c -o obj/xdrfile.o -Wall -Wno-comment -Wno-unused-but-set-variable
#~ gcc -c src/util/io/xdrfile/xdrfile_xtc.c -o obj/xdrfile_xtc.o -Wall -Wno-unused-variable
#~ g++ -c -std=c++11 src/BAT_builder/topology.cpp -o obj/topology.o -Wall
#~ g++ -c src/BAT_builder/BAT_topology.cpp -o obj/BAT_topology.o -Wall
#~ g++ -c src/BAT_builder/BAT_trajectory.cpp -o obj/BAT_trajectory.o -Wall

#~ g++ src/BAT_builder/bat.cpp obj/xdrfile.o  obj/xdrfile_xtc.o obj/io.o obj/topology.o obj/BAT_topology.o obj/BAT_trajectory.o  obj/util.o -o bin/BAT_builder -Wall
g++ --std=c++11 -O3 src/process_output/get_values_from_PAR.cpp obj/io.o obj/util.o -o bin/get_values_from_PAR -Wall

#bin/BAT_builder -t complexes/1UGH/1UGH.top -x complexes/1UGH/1UGH.xtc -o complexes/output/1UGH.bat -bb "CA C N H1 O1"

#~ nvcc --std=c++11 -O3 -Xptxas -O3 -gencode=arch=compute_61,code=\"sm_61,compute_61\" src/PARENT_GPU/PARENT_GPU.cu obj/io.o obj/util.o -o bin/PARENT_GPU \
#~ && bin/PARENT_GPU -f devel/traj/UBM2_1.bat -o output/UBM2_1.par -b 50 \
#~ && bin/get_values_from_PAR -p output/UBM2_1.par --short 


nvcc --std=c++11 -O3 -Xptxas -O3 -gencode=arch=compute_61,code=\"sm_61,compute_61\" src/MIST_GPU/MIST_GPU.cu obj/io.o obj/util.o -o bin/MIST_GPU \
&& bin/MIST_GPU -f output/UBM2_1.par -o output/UBM2_1_MIST.par

