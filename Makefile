# Script for GNU Make to build the PARENT_GPU program suite
# Copyright (C) 2020  Markus Fleck

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

ifndef CUDA_ARCH
    CUDA_ARCH = 61
endif

CXX = g++
CXXFLAGS = -O3 -Wall

all : bin/BAT_builder bin/convert_BAT_to_GBAT bin/PARENT_GPU bin/MIST_openMP bin/MIST_GPU bin/get_values_from_PAR

clean :
	- rm -r bin
	- rm -r obj

obj :
	mkdir obj

bin :
	mkdir bin



bin/BAT_builder : src/BAT_builder/bat.cpp obj/io.o obj/topology.o obj/BAT_topology.o obj/BAT_trajectory.o obj/util.o | bin
	$(CXX) src/BAT_builder/bat.cpp obj/io.o obj/topology.o obj/BAT_topology.o obj/BAT_trajectory.o obj/util.o -lgromacs -o bin/BAT_builder $(CXXFLAGS)
    
obj/topology.o: src/BAT_builder/topology.cpp src/BAT_builder/topology.h| obj
	$(CXX) -c -std=c++11 src/BAT_builder/topology.cpp -o obj/topology.o $(CXXFLAGS)
    
obj/BAT_topology.o: src/BAT_builder/BAT_topology.cpp src/BAT_builder/BAT_topology.h | obj
	$(CXX) -c -std=c++11 src/BAT_builder/BAT_topology.cpp -o obj/BAT_topology.o $(CXXFLAGS)
    
obj/BAT_trajectory.o: src/BAT_builder/BAT_trajectory.cpp src/BAT_builder/BAT_trajectory.h | obj
	- $(CXX) -c -std=c++11 src/BAT_builder/BAT_trajectory.cpp -o obj/BAT_trajectory.o $(CXXFLAGS)
	if [ ! -f obj/BAT_trajectory.o ]; then echo "\nTRYING COMPATIBILITY COMPILATION.\n";$(CXX) -c -std=c++11 src/BAT_builder/BAT_trajectory.cpp -o obj/BAT_trajectory.o $(CXXFLAGS) -DCOMPAT; fi

obj/io.o: obj/io_binary.o obj/io_text.o obj/Arg_Parser.o | obj
	ld -r obj/io_binary.o obj/io_text.o obj/Arg_Parser.o -o obj/io.o
    
obj/io_binary.o: obj/io_binary_tmp.o obj/Entropy_Matrix.o obj/Bat_File.o | obj
	ld -r obj/io_binary_tmp.o obj/Entropy_Matrix.o obj/Bat_File.o -o obj/io_binary.o
    
obj/io_binary_tmp.o: src/util/io/io_binary.cpp src/util/io/io_binary.h | obj
	$(CXX) -c src/util/io/io_binary.cpp -o obj/io_binary_tmp.o $(CXXFLAGS)
    
obj/Entropy_Matrix.o: src/util/classes/Entropy_Matrix.cpp src/util/classes/Entropy_Matrix.h | obj    
	$(CXX) -std=c++11 -c src/util/classes/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o $(CXXFLAGS)

obj/Bat_File.o: src/util/classes/Bat.cpp src/util/classes/Bat.h | obj
	$(CXX) -std=c++11 -c src/util/classes/Bat.cpp -o obj/Bat_File.o $(CXXFLAGS)

obj/io_text.o: src/util/io/io_text.cpp src/util/io/io_text.h | obj
	$(CXX) -c src/util/io/io_text.cpp -o obj/io_text.o $(CXXFLAGS)

obj/Arg_Parser.o: src/util/classes/Arg_Parser.cpp src/util/classes/Arg_Parser.h | obj
	$(CXX) -c src/util/classes/Arg_Parser.cpp -o obj/Arg_Parser.o $(CXXFLAGS)

obj/util.o: src/util/util.cpp src/util/util.h | obj
	$(CXX) -c src/util/util.cpp -o obj/util.o $(CXXFLAGS)
    
    
    

bin/convert_BAT_to_GBAT: src/BAT_builder/convert_BAT_to_GBAT.cpp obj/util.o obj/Arg_Parser.o obj/Bat_File.o | bin
	$(CXX) src/BAT_builder/convert_BAT_to_GBAT.cpp obj/util.o obj/Arg_Parser.o obj/Bat_File.o -o bin/convert_BAT_to_GBAT $(CXXFLAGS)
    

bin/get_values_from_PAR: src/process_output/get_values_from_PAR.cpp obj/io.o obj/util.o | bin
	$(CXX) --std=c++11 -O3 src/process_output/get_values_from_PAR.cpp obj/io.o obj/util.o -o bin/get_values_from_PAR $(CXXFLAGS)



bin/PARENT_GPU: src/PARENT_GPU/PARENT_GPU.cu src/PARENT_GPU/PARENT_GPU_kernels.cu obj/io.o src/util/io/io.h obj/util.o src/util/types.h | bin
	nvcc --std=c++11 -O3 -Xptxas -O3 -Xcompiler -O3,-Wall,-fopenmp  -gencode=arch=compute_$(CUDA_ARCH),code=\"sm_$(CUDA_ARCH),compute_$(CUDA_ARCH)\" src/PARENT_GPU/PARENT_GPU.cu obj/io.o obj/util.o -o bin/PARENT_GPU
    


bin/MIST_openMP: src/MIST_GPU/MIST_openMP.cpp obj/io.o obj/util.o | bin
	$(CXX) --std=c++11 -O3 -fopenmp src/MIST_GPU/MIST_openMP.cpp obj/io.o obj/util.o -o bin/MIST_openMP



bin/MIST_GPU: src/MIST_GPU/MIST_GPU.cu obj/io.o obj/util.o | bin
	nvcc --std=c++11 -O3 -Xptxas -O3 -Xcompiler -O3,-Wall,-fopenmp -gencode=arch=compute_$(CUDA_ARCH),code=\"sm_$(CUDA_ARCH),compute_$(CUDA_ARCH)\" src/MIST_GPU/MIST_GPU.cu obj/io.o obj/util.o -o bin/MIST_GPU



RAND := $(shell openssl rand -hex 12)
IN_NAME = test_system/UBQ_UBM2
OUT_NAME = ./$(RAND)/UBQ_UBM2
CPU_RAM = 0.1
GPU_RAM = 0.1
checks: all
	mkdir ./$(RAND)
	bin/BAT_builder -t $(IN_NAME).top -x $(IN_NAME).xtc -o $(OUT_NAME).bat -bb "CA C N H1 O1"
	bin/convert_BAT_to_GBAT -f $(OUT_NAME).bat -o $(OUT_NAME).gbat --ram $(CPU_RAM)
	bin/PARENT_GPU -f $(OUT_NAME).gbat -o $(OUT_NAME).par -b 50 --cpu_ram $(CPU_RAM) --gpu_ram $(GPU_RAM)
	bin/MIST_GPU -f $(OUT_NAME).par -o $(OUT_NAME)_MIST_GPU.par
	bin/MIST_openMP -f $(OUT_NAME).par -o $(OUT_NAME)_MIST_openMP.par
	bin/get_values_from_PAR -p ${OUT_NAME}.par --short 2>&1 > $(OUT_NAME)_MIE.txt
	bin/get_values_from_PAR -p ${OUT_NAME}_MIST_GPU.par --short 2>&1 > $(OUT_NAME)_MIST_GPU.txt
	bin/get_values_from_PAR -p ${OUT_NAME}_MIST_openMP.par --short 2>&1 > $(OUT_NAME)_MIST_openMP.txt
	echo; echo; echo; \
    CHECK_MIE=$$(diff $(OUT_NAME)_MIE.txt test_system/sample_output/sample_output_MIE.txt); \
    CHECK_MIST_GPU=$$(diff $(OUT_NAME)_MIST_GPU.txt test_system/sample_output/sample_output_MIST.txt); \
    CHECK_MIST_OPENMP=$$(diff $(OUT_NAME)_MIST_openMP.txt test_system/sample_output/sample_output_MIST.txt); \
    if [ "$$CHECK_MIE" = "" ]; then\
        echo "PARENT_GPU: pass"; else\
        echo "PARENT_GPU: FAIL"; fi;\
    if [ "$$CHECK_MIST_GPU" = "" ]; then\
        echo "MIST_GPU: pass"; else\
        echo "MIST_GPU: FAIL"; fi;\
    if [ "$$CHECK_MIST_OPENMP" = "" ]; then\
        echo "MIST_openMP: pass"; else\
        echo "MIST_openMP: FAIL"; fi;\
    rm -r ./$(RAND)


