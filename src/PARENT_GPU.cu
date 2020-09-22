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
		PRECISION *ram_start;
		PRECISION *ram_end; 
		PRECISION *gpu_ram_start;
	       	PRECISION *gpu_ram_end;
		unsigned long long int n_bytes;

	GPU_RAM_Block(unsigned char type, unsigned int first_dof, unsigned int n_dofs, 
			PRECISION *ram_start, PRECISION *gpu_ram_start, unsigned long long int n_bytes){
		
		this->type = type;
		this->first_dof = first_dof;
		this->n_dofs = n_dofs;
		this->last_dof = first_dof + n_dofs - 1;
		this->ram_start = ram_start;
		this->ram_end = ram_start + (n_bytes - 1) / sizeof(PRECISION);
		this->gpu_ram_start = gpu_ram_start;
		this->gpu_ram_end = gpu_ram_start + (n_bytes - 1) / sizeof(PRECISION);
		this->n_bytes = n_bytes;

	};


	void load_GPU(){
		gpuErrchk(cudaMemcpy(ram_start, gpu_ram_start, n_bytes, cudaMemcpyHostToDevice));
	};	


};








int main(){


return 0;
}



