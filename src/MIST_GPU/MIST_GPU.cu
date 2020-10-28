#define PRECISION double
#define THREADS_PREPARATION_GPU 128

#include <iostream>
#include <sys/time.h>
#include <vector>

#include "../util/util.h"
#include "../util/io/io.h"
#include "../util/types.h"
#include "../util/classes/Entropy_Matrix.h"
#include "../util/classes/My_Error.cpp"

using namespace std;

struct Dof_Pair {
    unsigned int id1;
    unsigned int id2;
    PRECISION value;
};

struct Mutual {
    PRECISION value;
    unsigned int id;
};


// this function fetches the respective mutual information values to perform the reduction on
__global__ void prepare_data_gpu (PRECISION* mut_mat, unsigned int n_dofs, unsigned int* processed, unsigned int* unprocessed, 
                                    unsigned int n_unprocessed, unsigned int n_elements, Mutual* prepared_data_out){
                                    
    unsigned int idx_g = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx_g < n_elements){
        
        //Every possible pair of processed and unprocessed dofs needs to be treated. The mapping of thread/processed/unprocessed indices is 
        //idx_g = idx_processed * n_unprocessed + idx_unprocessed; Therefore:
        unsigned int idx_processed = idx_g / n_unprocessed;
        unsigned int idx_unprocessed = idx_g - idx_processed * n_unprocessed;
    
        idx_processed = processed[idx_processed]; // reuse the register variables to store the actual indices tied to the dofs
        idx_unprocessed = unprocessed[idx_unprocessed];
    
        // now we need to calculate the index of the corresponding value in the mutual information matrix, which is a symmetric matrix.
        // It is stored in a compact fashion. If n_dofs would be 4, then the matrix stores the values referring to the dof index pairs like 
        // (0,1), (0,2), (0,3), (1,2), (1,3), (2,3). The Formula we will use (after minor calculus) is 
        // idx_mut_mat = binomial(n_dofs, 2) - binomial(n_dofs - idx_smaller, 2) + idx_larger - idx_smaller - 1  
        unsigned int idx_larger = (idx_processed > idx_unprocessed) ? idx_processed : idx_unprocessed;
        unsigned int idx_smaller = (idx_processed > idx_unprocessed) ? idx_unprocessed : idx_processed;
    
        unsigned int idx_mut_mat = idx_smaller * (2 * n_dofs - 3 - idx_smaller) / 2 + idx_larger - 1;
        
        // Finally, we store the mutual information value as well the index in the (virtual) all-to-all processed-unprocessed matrix
        prepared_data_out[idx_g] = {mut_mat[idx_mut_mat], idx_g};
                                    
    }                                    
}

// function to reduce the maximum of the prepared mutual information matrix. blockDim.x always needs to be equal warpSize.
__global__ void reduce_gpu(Mutual* input, Mutual* output, unsigned int n_elements){
    
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x; // get the id in the whole grid
    
    if(idx < n_elements){ //  make sure not to process data which isn't there 
        Mutual max_mut = input[idx]; // set the maximum to the frist value read in
        unsigned int mask = __ballot_sync(0xffffffff, idx < n_elements); // sync the threads and get a bitmask which makes sure that only the correct threads participate in the following warp primitives 

        #pragma unroll
        for(unsigned int delta = warpSize / 2; delta > 0; delta >>= 1){ // reduce the warp
            
            // get a comparison value via warp primitives
            Mutual compare_mut;
            compare_mut.value = __shfl_down_sync(mask, max_mut.value, delta);
            compare_mut.id = __shfl_down_sync(mask, max_mut.id, delta);
            
            if (compare_mut.value > max_mut.value){ // if the mutual information is larger, set a new maximum
                max_mut.value = compare_mut.value;
                max_mut.id = compare_mut.id;
            }
        }
        if(!threadIdx.x){ // only the first thread in the warp (threadIdx.x == 0) has the correct value
            output[blockIdx.x] = max_mut; // and writes it to the output array. The max_mut.id value will be used for identifying the position in the cpu to shift the according dof from the unprocessed to the processed vector.
        }
    }
}


// this program uses Prim's algorithm to find the maximum information spanning tree
int main(int argc, char* argv[]){
    try{
        // start the stopwatch for the execution time
        timeval tv_start, tv_start_calc, tv_end, tv_end_calc;
        gettimeofday(&tv_start, NULL);
        
        // set the GPU up
        int device_count;
        gpuErrchk ( cudaGetDeviceCount(&device_count) );
        unsigned int device = 0; //TODO: implement choices for the GPU. For now, call the program via e. g. export CUDA_VISIBLE_DEVICES=1;
        cout << "Found " << device_count << " CUDA device(s). Chose CUDA device number "<< device << "." << endl;
        struct cudaDeviceProp prop;
        gpuErrchk (cudaGetDeviceProperties(&prop, device) );
        cout << "Device name: " << prop.name << endl;
        cout << "CUDA capability: " << prop.major << "." << prop.minor << endl;
        cout << "Global memory: " << prop.totalGlobalMem / 1024 / 1024 << " MiB" << endl;
        cout << "Shared memory per block: " << prop.sharedMemPerBlock / 1024 << " kiB" << endl;
        cout << "Maximum threads per block dimension: " << prop.maxThreadsDim[0] << " " << prop.maxThreadsDim[1] << " " << prop.maxThreadsDim[2] << endl;
        cout << "Maximum blocks per grid dimension: " << prop.maxGridSize[0] << " " << prop.maxGridSize[1] << " " << prop.maxGridSize[2] << endl;
        cout << "Warp size: "<< prop.warpSize << endl << endl;

        
        // parse the command line arguments
        if (argc != 5){
            cerr << "USAGE: " << argv[0] << " -f input.par -o output.par" << endl;
            exit(EXIT_FAILURE);
        }
        Arg_Parser arg_parser(argc, argv);
        if (!arg_parser.exists("-f") || !arg_parser.exists("-o")){
            cerr << "USAGE: " << argv[0] << " -f input.par -o output.par" << endl;
            exit(EXIT_FAILURE);
        }
        if (   strcmp( arg_parser.get_ext( arg_parser.get("-f") ), "par" ) 
            || strcmp( arg_parser.get_ext( arg_parser.get("-o") ), "par" ) ){
            cerr << "USAGE: " << argv[0] << " -f input.par -o output.par" << endl;
            exit(EXIT_FAILURE);
        }
        
        
        //initialize the variables
        Entropy_Matrix ent_mat( arg_parser.get("-f") );
        unsigned int n_bonds = ent_mat.getNBonds();
        unsigned int n_angles = ent_mat.getNAngles();
        unsigned int n_dihedrals = ent_mat.getNDihedrals();
        unsigned int n_type[3] = {n_bonds, n_angles, n_dihedrals};
        unsigned int n_atoms = n_dihedrals + 3;
        unsigned int n_dofs = 3 * n_atoms - 6; // the number of internal degrees of freedom in 3D space
        unsigned int n_mut = n_dofs * (n_dofs - 1) / 2; // all possible dof pairs
        
        vector<unsigned int> h_processed; // Using Prim's algorithm, we need to keep track of which dofs have already been processd 
        vector<unsigned int> h_unprocessed;
        unsigned int *d_processed, *d_unprocessed;
        gpuErrchk( cudaMalloc( (void**) &d_processed, (n_dofs - 1) * sizeof(unsigned int) ) ); // we start with one dof in processed, so n_dofs - 1 is the maximum array size required
        gpuErrchk( cudaMalloc( (void**) &d_unprocessed, (n_dofs - 1) * sizeof(unsigned int) ) );
        for(unsigned int id = 1; id < n_dofs; id++){
            h_unprocessed.push_back(id);
        }
        h_processed.push_back(0);
        
        // The maximum number of mutual information terms to be reduced is reached when there are equally many dofs
        // in h_processed and h_unprocessed. The following formula also works for odd n_dofs.
        unsigned int n_elements_max = (n_dofs * n_dofs) / 4; 
        
        Mutual *d_in_mat_original, *d_out_mat_original;
        gpuErrchk( cudaMalloc( (void**) &d_in_mat_original, n_elements_max * sizeof(Mutual) ) ); // allocate space for the maximum spaced input matrix
        gpuErrchk( cudaMalloc( (void**) &d_out_mat_original, ( (n_elements_max + (prop.warpSize - 1) ) / prop.warpSize ) * sizeof(Mutual) ) ); // as well as the maximum spaced output matrix
        Mutual *d_in_mat, *d_out_mat; // the reduction will be done in more than one step. These pointers will be set to the current input and output array d_in_mat_original d_out_mat_original in "leap-frog style", see below
        
        Dof_Pair* mi_mist = new Dof_Pair[n_dofs - 1]; // this array will be filled with the results of the reductions
        
        // put all mutual information terms in one array representing a symmetric matrix in the usual compact style
        unsigned int counter = 0;
        PRECISION* h_mut_mat = new PRECISION[n_mut];
        for (unsigned char type1 = 0; type1 < 3; type1++){ // for bonds, angles and dihedrals of the first member of the dof pair
            for(unsigned int idx1 = 0; idx1 < n_type[type1]; idx1++){ // and all dofs of the current type of the first member of the dof pair
                for (unsigned char type2 = type1; type2 < 3; type2++){ // and all "later" types for the second member of the dof pair
                    unsigned int start_idx2 = (type1 == type2) ? idx1 + 1 : 0; // for different types, all combinations are non-redundant. For equal types, only process "later" dofs for the second member of the dof pair
                    for(unsigned int idx2 = start_idx2; idx2 < n_type[type2]; idx2++ ){ // and all "later" dofs for the second member of the dof pair
                        h_mut_mat[counter] = ent_mat.getMutual(type1, type2, idx1 + 1, idx2 + 1); //fill the matrix //TODO: change indexing in Entropy_Matrix to start from 0
                        counter++;
                    }
                }
            }
        }
        PRECISION* d_mut_mat; // allocate the same array on the GPU and transfer the data
        gpuErrchk( cudaMalloc( (void**) &d_mut_mat, counter * sizeof(PRECISION) ) );
        gpuErrchk( cudaMemcpy(d_mut_mat, h_mut_mat, counter * sizeof(PRECISION), cudaMemcpyHostToDevice ) );
        
        
        unsigned int n_elements = n_dofs - 1; // during the first reduction, ther will be one element in h_processed and (n_dofs - 1) elements in h_unprocessed. Thus the first matrix to be reduced to be reduced contains 1 * (n_dofs - 1) == n_dofs - 1 elements
        gettimeofday(&tv_start_calc, NULL); // start the calculation timer
        while( h_unprocessed.size() > 0){ //until all dofs have been processed
        
            Mutual max; // start with a zero maximum (mutual information is non-negative definite)
            max.value = 0.0;
            max.id = 0;
        
            gpuErrchk( cudaMemcpy(d_processed, h_processed.data(), h_processed.size() * sizeof(unsigned int), cudaMemcpyHostToDevice ) ); //copy the current processed/unprocessed arrays to the GPU
            gpuErrchk( cudaMemcpy(d_unprocessed, h_unprocessed.data(), h_unprocessed.size() * sizeof(unsigned int), cudaMemcpyHostToDevice ) );
            
            // the preparation needs n_elements / THREADS_PREPARATION_GPU blocks plus one more block if this division has a remainder. THREADS_PREPARATION_GPU might be tweaked for performance, should be a power of 2, larger than the warpSize and smaller than the maximum number of threads of the GPU  
            unsigned int n_blocks_preparation = (n_elements + (THREADS_PREPARATION_GPU - 1) ) / THREADS_PREPARATION_GPU; 
            // prepare the current data in prcessed/unprocessed. The output is an array of Mutual structs, containing the mutual information values to be reduced. The returned ids are just the positions in the matrix to be reduced. 
            prepare_data_gpu<<<n_blocks_preparation, THREADS_PREPARATION_GPU>>>(d_mut_mat, n_dofs, d_processed, d_unprocessed, h_unprocessed.size(), n_elements, d_in_mat_original); 
        
            d_in_mat = d_in_mat_original; // set the input/output arrays for the first pass of the while loop
            d_out_mat = d_out_mat_original;

        
            while(n_elements > 1){ // reduce until there is only one element left
                
                unsigned int n_blocks = (n_elements + (prop.warpSize - 1) ) / prop.warpSize; // the reduction is done in a per-warp style, so there should be n_elements / warpSize blocks plus one additional block if this division has a remainder
                reduce_gpu<<< n_blocks, prop.warpSize >>> (d_in_mat, d_out_mat, n_elements); // perform the reduction
                gpuErrchk( cudaPeekAtLastError() ); // check for errors and synchronize
                gpuErrchk( cudaDeviceSynchronize() );
            
                n_elements = n_blocks; // in the next pass, there are as many elements left to be reduced as there have been blocks during the last reduction
                if(n_elements > 1){ // if there is is another reduction to be performed, swap the input and output arrays. Otherwise, leave the output as the ouput for the upcoming cudaMemcpy
                    Mutual* d_tmp_mat = d_in_mat;
                    d_in_mat = d_out_mat;
                    d_out_mat = d_tmp_mat;
                }
            }
            gpuErrchk( cudaMemcpy(&max, d_out_mat, sizeof(Mutual), cudaMemcpyDeviceToHost ) ); // copy the result when the reduction has finished
            
            
            unsigned int idx_processed = max.id / h_unprocessed.size(); // back-calculate the indices in the processed and unprocessed arrays from the id of the result
            unsigned int idx_unprocessed = max.id - h_unprocessed.size() * idx_processed;

            Dof_Pair tmp_dof_pair = {h_processed[idx_processed], h_unprocessed[idx_unprocessed], max.value}; // create a dof pair from the result
            mi_mist[h_processed.size() - 1] = tmp_dof_pair; // and store it in the final result array
            
            h_processed.push_back(h_unprocessed[idx_unprocessed]); // move the unprocessed dof of the result pair to the processed array 
            h_unprocessed.erase(h_unprocessed.begin() + idx_unprocessed);
            n_elements = h_processed.size() * h_unprocessed.size();
            
            if(h_unprocessed.size() % 100 == 0) cout<<"Degress of freedom processed: "<< (h_processed.size() - 1) / double(n_dofs - 1) * 100<<"%"<<endl; // write a progress message
        }
        gettimeofday(&tv_end_calc, NULL); // stop the calculation timer
        
        
        // zero out the mutual information terms of the input matrix (analogously to the filling of the h_mut_mat array above)
        for (unsigned char type1 = 0; type1 < 3; type1++){
            for(unsigned int idx1 = 0; idx1 < n_type[type1]; idx1++){
                for (unsigned char type2 = type1; type2 < 3; type2++){
                    unsigned int start_idx2 = (type1 == type2) ? idx1 + 1 : 0;
                    for(unsigned int idx2 = start_idx2; idx2 < n_type[type2]; idx2++ ){
                        ent_mat.setMutual(type1, type2, idx1 + 1, idx2 + 1, 0.0); //TODO: change indexing in Entropy_Matrix to start from 0
                    }
                }
            }
        }
        
        // then set only the mutual information values of the result in the matrix  
        for (unsigned int i = 0; i < n_dofs - 1; i++){ 
            Dof dof1 = get_dof_from_global_id(mi_mist[i].id1, n_dihedrals); // this helper function and struct is declared/defined in util.cpp and util.h 
            Dof dof2 = get_dof_from_global_id(mi_mist[i].id2, n_dihedrals); // the id of a Dof is local to its type, meaning e. g. angle indices start from 0 different to their global ids, which start at n_bonds
            ent_mat.setMutual(dof1.type, dof2.type, dof1.id + 1, dof2.id + 1, mi_mist[i].value);  //use the obtained Dof sructs to set the result in ent_mat//TODO: change indexing in Entropy_Matrix to start from 0
        }
        
        ent_mat.write( arg_parser.get("-o") ); // write the result to disk
        
        gettimeofday(&tv_end, NULL); // report calculation times
        cout << endl << endl;
        cout << "Calculation time: " << tv_end_calc.tv_sec + 1e-6 * tv_end_calc.tv_usec - tv_start_calc.tv_sec - 1e-6 * tv_start_calc.tv_usec << endl; // of the reduction part
        cout << "Total execution time: " << tv_end.tv_sec + 1e-6 * tv_end.tv_usec - tv_start.tv_sec - 1e-6 * tv_start.tv_usec << endl; // and the whole program
        cout << "PROGRAM FINISHED SUCCESSFULLY." << endl;
        
        cudaDeviceReset(); // clean up
        return 0;
    } 
    catch (My_Error my_error) { // catch errors using the convenient C++ style
        cerr << my_error.what() << endl;
        exit(EXIT_FAILURE);
    }
    catch (...) {
        cerr << "AN UNIDENTIFIED ERROR HAS OCCURRED! ABORTING." << endl;
        exit(EXIT_FAILURE);
    }
    
}


