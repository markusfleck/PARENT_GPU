#include <iostream>
#include <cmath>
#include <vector>
#include "../util/types.h"
#define PRECISION double


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
		unsigned int type_id_start[3];
		unsigned int type_id_end[3]; //inclusive
		unsigned int type_n_dofs[3];
		unsigned int gpu_ram_blocks_per_type[3];
		unsigned int n_dihedrals;
		unsigned int gpu_dofs_per_block;

		vector<GPU_RAM_Block> blocks;
	
		CPU_RAM_Block(unsigned int dof_id_start, unsigned int dof_id_end, unsigned int gpu_dofs_per_block, unsigned int n_dihedrals)
		{
			this->dof_id_start = dof_id_start;
			this->dof_id_end = dof_id_end;
			this->n_dihedrals = n_dihedrals;
			this->gpu_dofs_per_block = gpu_dofs_per_block;		
		
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

		void deploy(PRECISION* block_start, unsigned int n_frames){
			//TODO: load dofs from trajectory to block_start
			
			blocks.clear();
			for(unsigned char type = get_dof_type_from_id(dof_id_start, n_dihedrals); type <= get_dof_type_from_id(dof_id_end, n_dihedrals); type++){
				for(unsigned int i = 0; i < gpu_ram_blocks_per_type[type]; i++)
				{
					unsigned int block_id_start = type_id_start[type] + i *  gpu_dofs_per_block;
					unsigned int block_id_end = type_id_start[type] + (i+1) *  gpu_dofs_per_block - 1; 
					if (block_id_end > type_id_end[type]) block_id_end = type_id_end[type];
					PRECISION* cpu_ram_start = block_start + (block_id_start - dof_id_start) * n_frames;	
					blocks.push_back( *new GPU_RAM_Block(cpu_ram_start, block_id_start, block_id_end, n_frames, n_dihedrals) ); 
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
			tmp_result_entropy = extrema + 2 * n_dofs_total;
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

		RAM(unsigned long long int cpu_n_bytes, unsigned long long int gpu_n_bytes, unsigned int n_dihedrals, unsigned int n_frames, unsigned int n_bins)
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
		}


};




int main(){

unsigned long long int cpu_ram_available = static_cast<unsigned long long int>(1024)*1024*1024*4;
unsigned long long int gpu_ram_available = static_cast<unsigned long long int>(1024)*1024*1024*1;


RAM ram(cpu_ram_available, gpu_ram_available, 329, 8.0e5, 50);


cout<<ram.gpu_ram_layout->dofs_per_block<<" "<<ram.cpu_ram_layout->dofs_per_block<<endl;

ram.cpu_ram_layout->tmp_result_occupied_bins[0] = 1;


return 0;
}



