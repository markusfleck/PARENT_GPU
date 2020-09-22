#include <iostream>

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

class GPU_RAM_Block{
	public:
		unsigned char type;
		unsigned int first_dof;
		unsigned int last_dof;
		unsigned int n_dofs;
		PRECISION *cpu_ram_start;
		PRECISION *cpu_ram_end; 
		PRECISION *gpu_ram_start;
	       	PRECISION *gpu_ram_end;
		unsigned long long int n_bytes;

		GPU_RAM_Block(unsigned char type, unsigned int first_dof, unsigned int n_dofs, 
				PRECISION *cpu_ram_start, PRECISION *gpu_ram_start, unsigned long long int n_bytes)
		{
			this->type = type;
			this->first_dof = first_dof;
			this->n_dofs = n_dofs;
			this->last_dof = first_dof + n_dofs - 1;
			this->cpu_ram_start = cpu_ram_start;
			this->cpu_ram_end = cpu_ram_start + (n_bytes - 1) / sizeof(PRECISION);
			this->gpu_ram_start = gpu_ram_start;
			this->gpu_ram_end = gpu_ram_start + (n_bytes - 1) / sizeof(PRECISION);
			this->n_bytes = n_bytes;
		}


		void load_GPU()
		{
			gpuErrchk(cudaMemcpy(cpu_ram_start, gpu_ram_start, n_bytes, cudaMemcpyHostToDevice));
		}	



};



class CPU_RAM_Block{
	public:
		PRECISION *ram_start;
                PRECISION *ram_end;
                unsigned long long int n_bytes;
		GPU_RAM_Block *blocks;
	
		CPU_RAM_Block(PRECISION *ram_start, unsigned long long int n_bytes)
		{
			this->ram_start = ram_start;
			this->n_bytes = n_bytes;
			this->ram_end = ram_start + (n_bytes - 1) / sizeof(PRECISION);
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

		RAM(unsigned long long int cpu_n_bytes, unsigned long long int gpu_n_bytes, unsigned int n_dihedrals, unsigned int n_frames)
		{
			this->cpu_ram_start = new char [cpu_n_bytes];
			this->cpu_ram_end = cpu_ram_start + cpu_n_bytes - 1;
                        this->cpu_n_bytes = cpu_n_bytes; 
			gpuErrchk( cudaMalloc((void**) &gpu_ram_start, gpu_n_bytes) );
                        this->gpu_ram_end = gpu_ram_start + gpu_n_bytes - 1;
                        this->gpu_n_bytes = gpu_n_bytes;
		}


};

int main(){

long long unsigned int cpu_ram_available = (long long unsigned int)(1024)*1024*1024*4;
long long unsigned int gpu_ram_available = (long long unsigned int)(1024)*1024*1024*1;

RAM ram(cpu_ram_available, gpu_ram_available, 10, 10);

return 0;
}



