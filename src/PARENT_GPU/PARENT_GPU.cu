#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>

#include "../util/io/io.h"
#include "../util/types.h"

#define PRECISION double
#define MODFITNBINS 100
#define WARPMULTIPLES 1
#define MEMORY_USAGE 0.8 //GPU
#define RAM_USAGE 0.8 //CPU
#define DEBUG false

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t return_code, const char *file, int line)
{
   if (return_code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(return_code), file, line);
      exit(return_code);
   }
}

using namespace std;



struct Bins {
	unsigned int bonds1D;
	unsigned int angles1D;
	unsigned int dihedrals1D;
	unsigned int bonds2D;
	unsigned int angles2D;
	unsigned int dihedrals2D;
};


unsigned char get_dof_type_from_id(unsigned int dof_id, unsigned int n_dihedrals)
{
	if(dof_id < n_dihedrals + 2) return TYPE_B;
	if(dof_id < 2 * n_dihedrals + 3) return TYPE_A;
	return TYPE_D;
}



unsigned int get_min_id_for_type(unsigned char type, unsigned int n_dihedrals)
{
	switch(type) {
  		case TYPE_B:
    			return 0;
  		case TYPE_A:
    			return n_dihedrals + 2;
		case TYPE_D:
    			return 2 * n_dihedrals + 3;
	}
	return 42;
}



unsigned int get_max_id_for_type(unsigned char type, unsigned int n_dihedrals)
{
	switch(type) {
  		case TYPE_B:
    			return n_dihedrals + 1;
  		case TYPE_A:
    			return 2 * n_dihedrals + 2;
		case TYPE_D:
    			return 3 * n_dihedrals + 2;
	}
	return 42;
}




class GPU_RAM_Block
{
	public:
		unsigned char type;
		unsigned int first_dof;
		unsigned int last_dof;
		unsigned int n_dofs;
		unsigned long long int n_bytes;
		PRECISION* cpu_ram_start;

		GPU_RAM_Block(PRECISION* cpu_ram_start, unsigned int first_dof, unsigned int last_dof, unsigned int n_frames, unsigned int n_dihedrals)
		{
			this->cpu_ram_start = cpu_ram_start;
			this->type = get_dof_type_from_id(first_dof, n_dihedrals);
			this->first_dof = first_dof;
			this->last_dof = last_dof;
			n_dofs = last_dof - first_dof + 1;
			n_bytes = n_dofs * n_frames * sizeof(PRECISION);
		}


		void deploy(PRECISION* gpu_ram_start)
		{
			gpuErrchk(cudaMemcpy(cpu_ram_start, gpu_ram_start, n_bytes, cudaMemcpyHostToDevice));
		}	
};



class CPU_RAM_Block
{
	public:
		unsigned int dof_id_start;
		unsigned int dof_id_end; //inclusive
		unsigned int n_dofs;
		unsigned int type_id_start[3];
		unsigned int type_id_end[3]; //inclusive
		unsigned int type_n_dofs[3];
		unsigned int gpu_ram_blocks_per_type[3];
		unsigned int n_dihedrals;
		unsigned int gpu_dofs_per_block;
		unsigned int n_frames;
		ifstream* bat_file;
		streamoff file_dofs_begin; 
		unsigned char precision_traj;
		vector<GPU_RAM_Block> blocks;
		PRECISION* block_start;
	
		CPU_RAM_Block(unsigned int dof_id_start, unsigned int dof_id_end, unsigned int gpu_dofs_per_block, unsigned int n_dihedrals, 
				unsigned int n_frames, ifstream* bat_file, streamoff file_dofs_begin, unsigned char precision_traj)
		{
			this->dof_id_start = dof_id_start;
			this->dof_id_end = dof_id_end;
			this->n_dofs = dof_id_end - dof_id_start + 1;
			this->n_dihedrals = n_dihedrals;
			this->gpu_dofs_per_block = gpu_dofs_per_block;
			this->n_frames = n_frames;
			this->bat_file = bat_file; 
			this->file_dofs_begin = file_dofs_begin;
			this->precision_traj = precision_traj;		
		
			for(unsigned char type = get_dof_type_from_id(dof_id_start, n_dihedrals); type <= get_dof_type_from_id(dof_id_end, n_dihedrals); type++){
				if(dof_id_start < get_min_id_for_type(type, n_dihedrals)) 
				{
					type_id_start[type] = get_min_id_for_type(type, n_dihedrals); 
				}
				else
				{
					type_id_start[type] = dof_id_start;
				}
				
				if(dof_id_end > get_max_id_for_type(type, n_dihedrals)) 
				{
					type_id_end[type] = get_max_id_for_type(type, n_dihedrals);
				}
				else
				{
					type_id_end[type] = dof_id_end;
				} 
				type_n_dofs[type] = type_id_end[type] - type_id_start[type] + 1;
				
				gpu_ram_blocks_per_type[type] = type_n_dofs[type] / gpu_dofs_per_block;
				if(type_n_dofs[type] % gpu_dofs_per_block > 0) gpu_ram_blocks_per_type[type] += 1;
 
			}
		}

		unsigned char load_dofs(PRECISION* dof_bank)
		{
    			unsigned char fail=0;
			bat_file->seekg(file_dofs_begin);
			for(unsigned int frame=0; frame<n_frames; frame++)
			{
			        //to read a frame from the .bat trajectory
			        double ddummy[6];
				float fdummy[11];
				unsigned int dof_counter = 0;		       	
			    
			        if(frame % 100000 == 0) 
				{
			            cout<<"Reading frame "<<frame<<" and the following.\n";   //every 10000 frames issue an information to stdout
			            cout.flush();
			        }
			        
			
			        bat_file->read((char*)fdummy, 11*sizeof(float));//read time, precision and box vectors to dummies
			        fail=fail | (bat_file->rdstate() & std::ifstream::failbit);	
				unsigned char inc;		
			        if(precision_traj==1) 
				{
					inc = sizeof(double);//if trajectory is in double precision
				}
				else
				{
					inc = sizeof(float);
				}	
			        bat_file->read((char*)ddummy, 6 * inc);//external coordinates to dummies
			        fail=fail | (bat_file->rdstate() & std::ifstream::failbit);
			            
			        bat_file->read((char*)ddummy, inc);//read the lengths of the two bonds connecting the root atoms (internal coordinates)
			        fail=fail | (bat_file->rdstate() & std::ifstream::failbit);
			        if( (dof_counter >= dof_id_start) && (dof_counter <= dof_id_end) ) dof_bank[(dof_counter - dof_id_start) * n_frames + frame] = ddummy[0];
				dof_counter++;

				bat_file->read((char*)ddummy, inc);//read the lengths of the two bonds connecting the root atoms (internal coordinates)
			        fail=fail | (bat_file->rdstate() & std::ifstream::failbit);
			        if( (dof_counter >= dof_id_start) && (dof_counter <= dof_id_end) ) dof_bank[(dof_counter - dof_id_start) * n_frames + frame] = ddummy[0];
				dof_counter++;
			            
			        bat_file->read((char*)ddummy, inc);//and the angle between the two rootbonds (internal coordinates)
			        fail=fail | (bat_file->rdstate() & std::ifstream::failbit);
			        if( (dof_counter >= dof_id_start) && (dof_counter <= dof_id_end) ) dof_bank[(dof_counter - dof_id_start) * n_frames + frame] = ddummy[0];
				dof_counter++;  
			        for(int i = 0; i < n_dihedrals; i++) 
			        { //then for all dihedrals in the system
					bat_file->read((char*)ddummy, inc);//read the bondlength between the last two atoms in the dihedral
			                fail=fail | (bat_file->rdstate() & std::ifstream::failbit);		              
			               	if( (dof_counter >= dof_id_start) && (dof_counter <= dof_id_end) ) dof_bank[(dof_counter - dof_id_start) * n_frames + frame] = ddummy[0];
					dof_counter++;
 
			                bat_file->read((char*)ddummy, inc);//read the angle between the last threee atoms of the dihedral#
			                fail=fail | (bat_file->rdstate() & std::ifstream::failbit);
			                if( (dof_counter >= dof_id_start) && (dof_counter <= dof_id_end) ) dof_bank[(dof_counter - dof_id_start) * n_frames + frame] = ddummy[0];
					dof_counter++;

			                bat_file->read((char*)ddummy, inc);//and the value of the dihedral itself
			                fail=fail | (bat_file->rdstate() & std::ifstream::failbit);
			                if( (dof_counter >= dof_id_start) && (dof_counter <= dof_id_end) ) dof_bank[(dof_counter - dof_id_start) * n_frames + frame] = ddummy[0];
					dof_counter++;
				}
			}
		   
		    return fail; //if anything failed return a 1, otherwise a 0
			

		}

		void deploy(PRECISION* block_start){
			this->block_start = block_start;
			load_dofs(block_start);
			
			blocks.clear();
			for(unsigned char type = get_dof_type_from_id(dof_id_start, n_dihedrals); type <= get_dof_type_from_id(dof_id_end, n_dihedrals); type++)
			{
				for(unsigned int i = 0; i < gpu_ram_blocks_per_type[type]; i++)
				{
					unsigned int block_id_start = type_id_start[type] + i *  gpu_dofs_per_block;
					unsigned int block_id_end = type_id_start[type] + (i + 1) *  gpu_dofs_per_block - 1; 
					if (block_id_end > type_id_end[type]) block_id_end = type_id_end[type];
					PRECISION* cpu_ram_start = block_start + (block_id_start - dof_id_start) * n_frames;	
					blocks.push_back( *new GPU_RAM_Block(cpu_ram_start, block_id_start, block_id_end, n_frames, n_dihedrals) ); 
				} 
			}		
		}

		void calculate_extrema(PRECISION* extrema)
		{
    			#pragma omp parallel
    			{
        			PRECISION tmpMin,tmpMax;

        			#pragma omp for
        			for(int j = 0; j < n_dofs; j++) 
				{ //for all dofs 
            				tmpMax = block_start[j * n_frames];
            				tmpMin = block_start[j * n_frames];
            				for(int i = 1; i < n_frames; i++) 
					{ //and all frames
                				if(block_start[j * n_frames + i] > tmpMax) 
						{
                    					tmpMax = block_start[j * n_frames + i];
                				}
                				if(block_start[j * n_frames + i] < tmpMin) 
						{
                    					tmpMin = block_start[j * n_frames + i];   //find the maximum and minmum values
                				}
            				}
           				if( (tmpMin < 0.0) || (tmpMax < 0.0) ) 
					{
                				cerr<<"ERROR: Degree of freedom "<< dof_id_start + j <<" is smaller than 0.0"<<endl;
                				exit(EXIT_FAILURE);
            				}
            				tmpMin -= 5e-9 * ( sizeof(PRECISION) == sizeof(float) ? 1e5 :1 );//and increase the boundaries a tiny bit
            				tmpMax += 5e-9 * ( sizeof(PRECISION) == sizeof(float) ? 1e5 :1 );
            				if(tmpMin < 0.0) 
					{
                				tmpMin = 0.0;
            				}
            				if ( (tmpMax-tmpMin) < 1.0e-4) 
					{
                				tmpMax+=0.05;
            				}
            				extrema[dof_id_start + j] = tmpMin;//and store the calculated values
            				extrema[dof_id_start + j + n_dofs] = tmpMax;
        			}
    			}
		}
};


class GPU_RAM_Layout
{
	public:
		unsigned int dofs_per_block;
		PRECISION* dof_block_1;
		PRECISION* dof_block_2;
		PRECISION* result_fwd_1;
		PRECISION* result_fwd_2;
		PRECISION* result_all2all;
		unsigned int* occupied_bins;
		unsigned int* histograms;

		GPU_RAM_Layout(unsigned int n_frames, unsigned int n_bins, unsigned long long int gpu_n_bytes, char* gpu_ram_start)
		{
			// calculate the maximum number of dofs (for one of the two dof_blocks) so that everything still fits into GPU RAM 
			double a = (2 * n_frames - 1) * sizeof(PRECISION) - sizeof(unsigned int);
			double b = (n_bins * n_bins + 2) * sizeof(unsigned int)	+ 2 * sizeof(PRECISION);
			this->dofs_per_block = (unsigned int)( (-a/2 + sqrt(a * a / 4 + gpu_n_bytes * b) ) /b);		

			// set the pointers for partitioning the GPU RAM according to the calculated dofs_per_block
			dof_block_1 = (PRECISION*) gpu_ram_start;
			dof_block_2 = dof_block_1 + dofs_per_block * n_frames;
			result_fwd_1 = dof_block_2 + dofs_per_block * n_frames;
			result_fwd_2 = result_fwd_1 + dofs_per_block * (dofs_per_block - 1);
			result_all2all =  result_fwd_2 + dofs_per_block * (dofs_per_block - 1);
			occupied_bins = (unsigned int*) (result_all2all + dofs_per_block * dofs_per_block);
			histograms = occupied_bins + (2 * dofs_per_block - 1) * dofs_per_block;
		}
};

class CPU_RAM_Layout
{
	public:
		unsigned int dofs_per_block;
		PRECISION* dof_block_1;
		PRECISION* dof_block_2;
		PRECISION* result_entropy1D;
		PRECISION* result_entropy2D;
		PRECISION* extrema;
		PRECISION* minima;
		PRECISION* maxima;
		PRECISION* tmp_result_entropy;
		unsigned int* tmp_result_occupied_bins;
		double* tmp_read;

		CPU_RAM_Layout(unsigned int n_frames, unsigned long long int cpu_n_bytes, char* cpu_ram_start, unsigned int gpu_dofs_per_block, unsigned int n_dofs_total)
		{
			// calculate the maximum number of dofs (for one of the two dof_blocks) so that everything still fits into CPU RAM
			dofs_per_block = (unsigned int)( ( cpu_n_bytes - n_dofs_total * ( ( n_dofs_total + 3 ) * sizeof(PRECISION) + sizeof(double) ) 
							+ ( 2 * gpu_dofs_per_block - 1 ) * gpu_dofs_per_block * ( sizeof(PRECISION) + sizeof(unsigned int) ) )
							/ ( 2 * n_frames * sizeof(PRECISION) ) );
			if(dofs_per_block < gpu_dofs_per_block)
			{
				cerr<<"WARNING: You probably have a GPU with a lot of RAM but your CPU RAM is rather small. ";
				cerr<<"I recommend to get more CPU RAM, as this should significantly enhance performance."<<endl;	
			}
			// if all dofs fit into RAM, still set up two blocks to be consistent with the algorithm
			if(2 * dofs_per_block >= n_dofs_total){
				dofs_per_block = n_dofs_total / 2;
				dofs_per_block += n_dofs_total % 2;
			}

			// set the pointers for partitioning the CPU RAM according to the calculated dofs_per_block
			dof_block_1 = (PRECISION*) cpu_ram_start;
			dof_block_2 = dof_block_1 + dofs_per_block * n_frames;
			result_entropy1D = dof_block_2 + dofs_per_block * n_frames;
			result_entropy2D = result_entropy1D + n_dofs_total;
			extrema = result_entropy2D + n_dofs_total * (n_dofs_total - 1) / 2;
			minima = extrema;
			maxima = extrema + n_dofs_total;
			tmp_result_entropy = maxima + n_dofs_total;
			tmp_result_occupied_bins = (unsigned int*)(tmp_result_entropy + (2 * gpu_dofs_per_block - 1) * gpu_dofs_per_block);
			tmp_read = (double*)(tmp_result_occupied_bins + (2 * gpu_dofs_per_block - 1) * gpu_dofs_per_block);
		}

};



class RAM{
	public:
		char *cpu_ram_start;
                char *cpu_ram_end;
                unsigned long long int cpu_n_bytes;
		char *gpu_ram_start;
                char *gpu_ram_end;
                unsigned long long int gpu_n_bytes;
		GPU_RAM_Layout* gpu_ram_layout;
		CPU_RAM_Layout* cpu_ram_layout;
		unsigned int n_dihedrals;
		unsigned int n_dofs_total;
		vector<CPU_RAM_Block> blocks;

		RAM(unsigned long long int cpu_n_bytes, unsigned long long int gpu_n_bytes, unsigned int n_dihedrals, unsigned int n_frames, 
			unsigned int n_bins, ifstream* bat_file, streamoff file_dofs_begin, unsigned char precision_traj)
		{
			cpu_ram_start = new char [cpu_n_bytes];
			cpu_ram_end = cpu_ram_start + cpu_n_bytes - 1;
                        this->cpu_n_bytes = cpu_n_bytes; 
			gpuErrchk( cudaMalloc((void**) &gpu_ram_start, gpu_n_bytes) );
                        gpu_ram_end = gpu_ram_start + gpu_n_bytes - 1;
                        this->gpu_n_bytes = gpu_n_bytes;
			gpu_ram_layout = new GPU_RAM_Layout(n_frames, n_bins, gpu_n_bytes, gpu_ram_start);
			this->n_dihedrals = n_dihedrals;
			n_dofs_total = 3 * (n_dihedrals + 1);
			cpu_ram_layout = new CPU_RAM_Layout(n_frames, cpu_n_bytes, cpu_ram_start, gpu_ram_layout->dofs_per_block, n_dofs_total);
			for(unsigned int i = 0; i < n_dofs_total; i+=cpu_ram_layout->dofs_per_block)
			{
				unsigned int end = i + cpu_ram_layout->dofs_per_block - 1;
				if(end > n_dofs_total - 1) end = n_dofs_total - 1;
				blocks.push_back(*new CPU_RAM_Block(i, end, gpu_ram_layout->dofs_per_block, n_dihedrals, n_frames, bat_file, file_dofs_begin, precision_traj));
			} 

			
		}


};



#include <algorithm>
char* getCmdOption(char ** begin, char ** end, const string & option)
{
    char ** itr = find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const string& option)
{
    return find(begin, end, option) != end;
}


int main(int argc, char *argv[]){

	//start the stopwatch for the execution time
	timeval tv_start,tv_end;
	gettimeofday (&tv_start, NULL);
	
	int deviceCount;    
	gpuErrchk(cudaGetDeviceCount(&deviceCount));
	
	
	unsigned int device=0;//TODO: implement choices for graphics card
	cout<<"Found "<<deviceCount<<" CUDA device(s). Chose CUDA device number "<<device<<"."<<endl;
	struct cudaDeviceProp prop;
	gpuErrchk(cudaGetDeviceProperties(&prop,device)); 	
	cout<<"Device name: "<<prop.name<<endl;
	cout<<"CUDA capability: "<<prop.major<<"."<<prop.minor<<endl;
	cout<<"Global memory: "<<prop.totalGlobalMem/1024/1024<<" MiB"<<endl;
	cout<<"Shared memory per block: "<<prop.sharedMemPerBlock/1024<<" kiB"<<endl;
	cout<<"Maximum threads per block dimension: "<<prop.maxThreadsDim[0]<<" "<<prop.maxThreadsDim[1]<<" "<<prop.maxThreadsDim[2]<<endl;
	cout<<"Maximum blocks per grid dimension: "<<prop.maxGridSize[0]<<" "<<prop.maxGridSize[1]<<" "<<prop.maxGridSize[2]<<endl;
	cout<<"Warp size: "<<prop.warpSize<<endl;
	cout<<endl<<endl;    
	
	//int threads_per_block = prop.warpSize*WARPMULTIPLES; 
	
	int precision_traj, n_frames;
	unsigned int n_bins;
	vector< vector <int> > dihedrals_top;
	vector <float> masses;
	vector <string> residues;
	vector <int> residueNumbers;
	vector <string> atomNames;
	vector <string> belongsToMolecule;
	
	if(argc!=7) 
	{
		cerr<<"USAGE: "<<argv[0]<<" -f input.bat -o entropy.par -b #bins\n";
	    	exit(EXIT_FAILURE);	
	}
	bool fail = false;
	fail |=  !cmdOptionExists(argv, argv+argc, "-f");
	fail |=  !cmdOptionExists(argv, argv+argc, "-o");
	fail |=  !cmdOptionExists(argv, argv+argc, "-b");
		          
	if(fail)
	{
		//check for correct command line options
		cerr<<"USAGE: "<<argv[0]<<" -f input.bat -o entropy.par -b #bins\n";
		exit(EXIT_FAILURE);
	}
	                
	string tmp1(getCmdOption(argv, argv+argc, "-f")); //first argument is the .bat trajectory file
	string tmp2(getCmdOption(argv, argv+argc, "-o")); //second argument is the .par output file
	char *ptr,*type1,*type2;
	char delimiter[] = ".";
	
	ptr = strtok((char*)tmp1.c_str(), delimiter);
	while(ptr != NULL) 
	{
		type1 = ptr;
		ptr = strtok(NULL, delimiter);
	}
	
	ptr = strtok((char*)tmp2.c_str(), delimiter);
	while(ptr != NULL) 
	{
		type2=ptr;
		ptr = strtok(NULL, delimiter);
	}
	if((strcmp(type1,"bat"))||(strcmp(type2,"par"))) 
	{
		//check for the extensions of the input and output file
		cerr<<"USAGE: "<<argv[0]<<" -f input.bat -o entropy.par -b #bins\n";
		exit(EXIT_FAILURE);
	}
	if(sscanf(getCmdOption(argv, argv+argc, "-b"),"%ud",&n_bins)!=1) 
	{
		 //read the number of bins and check for correctness
		cerr<<"ERROR: Could not read number of bins from command line! Aborting"<<endl;   		
		exit(EXIT_FAILURE);
	}
	
		
	ifstream bat_file;
	ofstream par_file;
	cout<<"Reading file "<<getCmdOption(argv, argv+argc, "-f")<<" .\n";
	//open the input/output files
	bat_file.open(getCmdOption(argv, argv+argc, "-f"), ios::binary | ios::in );
	par_file.open(getCmdOption(argv, argv+argc, "-o"),ios::binary | ios::out);
	if(!(bat_file.is_open())) 
	{
	    cerr<<"ERROR: Could not open file "<<getCmdOption(argv, argv+argc, "-f")<<" ! Aborting."<<endl;
	    exit(EXIT_FAILURE);
	}
	if(!(par_file.is_open())) 
	{
	    cerr<<"ERROR: Could not open file "<<getCmdOption(argv, argv+argc, "-f")<<" ! Aborting."<<endl;
	    exit(EXIT_FAILURE);
	}
	
	//and read the header of the trajectory
	if(read_BAT_header(&bat_file, &precision_traj, &n_frames, &dihedrals_top, &masses, &residues, &residueNumbers, &atomNames, &belongsToMolecule)!=0) 
	{
	    cerr<<"AN ERROR HAS OCCURED WHILE READING THE HEADER OF THE FILE " <<getCmdOption(argv, argv+argc, "-b")<<" . QUITTING PROGRAM.\n";
	    exit(EXIT_FAILURE);
	}
	unsigned int n_dihedrals=dihedrals_top.size();
	cout<<getCmdOption(argv, argv+argc, "-b")<<" specs:"<<endl;
	cout<<"Precision: "<<(precision_traj == 1 ? "double" : "single")<<" #Atoms: "<<n_dihedrals + 3<<" #Frames: "<<n_frames<<endl;  // ---------------------------------------------------
	streamoff file_dofs_begin = bat_file.tellg();
	
	//and write the header of the output (.par) file
	if(write_PAR_header(&par_file,n_dihedrals,precision_traj,n_frames,&dihedrals_top, &masses, 
				n_bins, n_bins, n_bins, n_bins, n_bins, n_bins, &residues, &residueNumbers, &atomNames, &belongsToMolecule) != 0) 
	{
	    cerr<<"AN ERROR HAS OCCURED WHILE WRITING THE HEADER OF THE FILE " <<getCmdOption(argv, argv+argc, "-o")<<" . QUITTING PROGRAM.\n";
	    exit(EXIT_FAILURE);
	}










	unsigned long long int cpu_ram_available = static_cast<unsigned long long int>(1024)*1024*1024*4;
	unsigned long long int gpu_ram_available = static_cast<unsigned long long int>(1024)*1024*1024*1;


	RAM ram(cpu_ram_available, gpu_ram_available, n_dihedrals, n_frames, n_bins, &bat_file, file_dofs_begin, precision_traj);
	//ram.cpu_ram_layout->result_entropy1D[0] = 42;
	//ram.cpu_ram_layout->result_entropy1D[-1] = 42;

	//cout<<ram.gpu_ram_layout->dofs_per_block<<" "<<ram.cpu_ram_layout->dofs_per_block<<endl;

	//ram.cpu_ram_layout->tmp_result_occupied_bins[0] = 1;


	for (unsigned int i = 0; i < ram.blocks.size() - 1; i++)
	{	cout<<"Deploying Block "<<i+1<<" to RAM bank 1."<<endl;
		ram.blocks[i].deploy(ram.cpu_ram_layout->dof_block_1);
		//cout<<ram.cpu_ram_layout->result_entropy1D[-1]<<" "<<ram.cpu_ram_layout->result_entropy1D[0]<<endl<<endl;
		for (unsigned int j = i + 1; j < ram.blocks.size(); j++)
		{
			//cout<<ram.blocks[i].dof_id_start<<" "<<ram.blocks[i].dof_id_end<<" "<<ram.blocks[j].dof_id_start<<" "<<ram.blocks[j].dof_id_end<<endl;
			cout<<"Deploy Block "<<j+1<<" to RAM bank 2."<<endl;
			ram.blocks[j].deploy(ram.cpu_ram_layout->dof_block_2);
			//cout<<ram.cpu_ram_layout->result_entropy1D[-1]<<" "<<ram.cpu_ram_layout->result_entropy1D[0]<<endl<<endl;

		}
	}


	//timings are written to stdout
	par_file.close();
	gettimeofday (&tv_end, NULL);
	cout<<endl<<endl;
	cout<<"Total execution time: "<<tv_end.tv_sec+1e-6 * tv_end.tv_usec-tv_start.tv_sec-1e-6 * tv_start.tv_usec<<endl;
 	cout<<"PROGRAM FINISHED SUCCESSFULLY."<<endl<<endl<<endl;
	
	cudaDeviceReset();
	return 0;
}



