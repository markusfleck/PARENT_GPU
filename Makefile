COPS = -O3 -Wall

all : bin/BAT_builder bin/convert_BAT_to_GBAT bin/PARENT_GPU bin/MIST_openMP bin/MIST_GPU bin/get_values_from_PAR

clean :
	rm -r bin obj

obj :
	mkdir obj

bin :
	mkdir bin



bin/BAT_builder : src/BAT_builder/bat.cpp obj/io.o obj/topology.o obj/BAT_topology.o obj/BAT_trajectory.o obj/util.o | bin
	g++ src/BAT_builder/bat.cpp obj/io.o obj/topology.o obj/BAT_topology.o obj/BAT_trajectory.o obj/util.o -lgromacs -o bin/BAT_builder $(COPS)
    
obj/topology.o: src/BAT_builder/topology.cpp src/BAT_builder/topology.h| obj
	g++ -c -std=c++11 src/BAT_builder/topology.cpp -o obj/topology.o $(COPS)
    
obj/BAT_topology.o: src/BAT_builder/BAT_topology.cpp src/BAT_builder/BAT_topology.h | obj
	g++ -c -std=c++11 src/BAT_builder/BAT_topology.cpp -o obj/BAT_topology.o $(COPS)
    
obj/BAT_trajectory.o: src/BAT_builder/BAT_trajectory.cpp src/BAT_builder/BAT_trajectory.h | obj
	g++ -c -std=c++11 src/BAT_builder/BAT_trajectory.cpp -o obj/BAT_trajectory.o $(COPS)

obj/io.o: obj/io_binary.o obj/io_text.o obj/Arg_Parser.o | obj
	ld -r obj/io_binary.o obj/io_text.o obj/Arg_Parser.o -o obj/io.o
    
obj/io_binary.o: obj/io_binary_tmp.o obj/Entropy_Matrix.o obj/Bat_File.o | obj
	ld -r obj/io_binary_tmp.o obj/Entropy_Matrix.o obj/Bat_File.o -o obj/io_binary.o
    
obj/io_binary_tmp.o: src/util/io/io_binary.cpp src/util/io/io_binary.h | obj
	g++ -c src/util/io/io_binary.cpp -o obj/io_binary_tmp.o -Wall
    
obj/Entropy_Matrix.o: src/util/classes/Entropy_Matrix.cpp src/util/classes/Entropy_Matrix.h | obj    
	g++ -std=c++11 -c src/util/classes/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o -Wall

obj/Bat_File.o: src/util/classes/Bat.cpp src/util/classes/Bat.h | obj
	g++ -std=c++11 -c src/util/classes/Bat.cpp -o obj/Bat_File.o -Wall

obj/io_text.o: src/util/io/io_text.cpp src/util/io/io_text.h | obj
	g++ -c src/util/io/io_text.cpp -o obj/io_text.o -Wall

obj/Arg_Parser.o: src/util/classes/Arg_Parser.cpp src/util/classes/Arg_Parser.h | obj
	g++ -c src/util/classes/Arg_Parser.cpp -o obj/Arg_Parser.o -Wall

obj/util.o: src/util/util.cpp src/util/util.h | obj
	g++ -c src/util/util.cpp -o obj/util.o -Wall
    
    
    

bin/convert_BAT_to_GBAT: src/BAT_builder/convert_BAT_to_GBAT.cpp obj/util.o obj/Arg_Parser.o obj/Bat_File.o | bin
	g++ src/BAT_builder/convert_BAT_to_GBAT.cpp obj/util.o obj/Arg_Parser.o obj/Bat_File.o -o bin/convert_BAT_to_GBAT -Wall
    

bin/get_values_from_PAR: src/process_output/get_values_from_PAR.cpp obj/io.o obj/util.o | bin
	g++ --std=c++11 -O3 src/process_output/get_values_from_PAR.cpp obj/io.o obj/util.o -o bin/get_values_from_PAR -Wall



bin/PARENT_GPU: src/PARENT_GPU/PARENT_GPU.cu src/PARENT_GPU/PARENT_GPU_kernels.cu obj/io.o src/util/io/io.h obj/util.o src/util/types.h | bin
	nvcc --std=c++11 -O3 -Xptxas -O3 -Xcompiler -O3,-Wall,-fopenmp  -gencode=arch=compute_75,code=\"sm_75,compute_75\" src/PARENT_GPU/PARENT_GPU.cu obj/io.o obj/util.o -o bin/PARENT_GPU
    


bin/MIST_openMP: src/MIST_GPU/MIST_openMP.cpp obj/io.o obj/util.o | bin
	g++ --std=c++11 -O3 -fopenmp src/MIST_GPU/MIST_openMP.cpp obj/io.o obj/util.o -o bin/MIST_openMP



bin/MIST_GPU: src/MIST_GPU/MIST_GPU.cu obj/io.o obj/util.o | bin
	nvcc --std=c++11 -O3 -Xptxas -O3 -Xcompiler -O3,-fopenmp -gencode=arch=compute_75,code=\"sm_75,compute_75\" src/MIST_GPU/MIST_GPU.cu obj/io.o obj/util.o -o bin/MIST_GPU



