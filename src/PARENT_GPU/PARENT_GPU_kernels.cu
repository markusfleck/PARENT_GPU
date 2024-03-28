// The GPU kernels of PARENT_GPU
// Copyright (C) 2020  Markus Fleck
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3 as 
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.



// #define FAST_MATH
#ifdef FAST_MATH
    #define COS __cosf
    #define SIN __sinf
    #define EXP __expf
    #define LOG __logf
#else 
    #define COS cos
    #define SIN sin
    #define EXP exp
    #define LOG log
#endif


__global__ void histo2D(PRECISION const *const __restrict__ traj1, PRECISION const *const __restrict__ traj2, const size_t numFrames,
                        unsigned int *const __restrict__ histo, const size_t n_bins,
                        PRECISION const inv_binSize1, PRECISION const inv_binSize2, PRECISION const min1,
                        PRECISION const min2) {

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t const offset = blockDim.x * gridDim.x;
    
  while (idx < numFrames) {
    atomicAdd(&histo[size_t((traj1[idx] - min1) * inv_binSize1) * n_bins +
                     size_t((traj2[idx] - min2) * inv_binSize2)],
              1);
    idx += offset;
  }
}


__global__ void histo2D_shared_block(PRECISION const *const __restrict__ traj1, PRECISION const *const __restrict__ traj2, const size_t numFrames,
                        unsigned int *const __restrict__ histo, const size_t n_bins,
                        PRECISION const inv_binSize1, PRECISION const inv_binSize2, PRECISION const min1,
                        PRECISION const min2) {
                        
  extern __shared__ unsigned int histo_block[];

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t const offset = blockDim.x * gridDim.x;
  size_t tid = threadIdx.x;
  size_t const n_bins_total = n_bins * n_bins;
                      
  while(tid < n_bins_total){
      histo_block[tid] = 0;
      tid += blockDim.x;
  }
    
  while (idx * 4 < numFrames - 3) {
    PRECISION4 traj_tmp1 = reinterpret_cast<PRECISION4 const*>(traj1)[idx];
    PRECISION4 traj_tmp2 = reinterpret_cast<PRECISION4 const*>(traj2)[idx];
    
    atomicAdd(&histo_block[size_t((traj_tmp1.x - min1) * inv_binSize1) * n_bins + size_t((traj_tmp2.x - min2) * inv_binSize2)], 1);
    atomicAdd(&histo_block[size_t((traj_tmp1.y - min1) * inv_binSize1) * n_bins + size_t((traj_tmp2.y - min2) * inv_binSize2)], 1);
    atomicAdd(&histo_block[size_t((traj_tmp1.z - min1) * inv_binSize1) * n_bins + size_t((traj_tmp2.z - min2) * inv_binSize2)], 1);
    atomicAdd(&histo_block[size_t((traj_tmp1.w - min1) * inv_binSize1) * n_bins + size_t((traj_tmp2.w - min2) * inv_binSize2)], 1);
    
    idx += offset;
  }

  
  long long int leftover = numFrames - 4 *idx;
  while(leftover > 0){
    atomicAdd(&histo_block[size_t((traj1[numFrames - leftover] - min1) * inv_binSize1) * n_bins + size_t((traj2[numFrames - leftover] - min2) * inv_binSize2)], 1);
    leftover--;
  }
  
  
  __syncthreads();
  
  tid = threadIdx.x;
  while(tid < (n_bins_total >> 1)){
    atomicAdd(&(reinterpret_cast<unsigned long long int*>(histo)[tid]), reinterpret_cast<unsigned long long int*>(histo_block)[tid]);
    tid += blockDim.x;
  }
  
}



__global__ void cu_baEnt(unsigned int *histo, const size_t numFrames,
                         const size_t bins1, const size_t bins2, PRECISION binSize1,
                         PRECISION binSize2, PRECISION min1, PRECISION min2,
                         PRECISION *plnpsum, unsigned int *occupbins) {

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t offset = blockDim.x * gridDim.x;
  while (idx < bins1 * bins2) {
    PRECISION blen = min1 + binSize1 / 2.0 + binSize1 * (idx / bins2);
    PRECISION theta = min2 + binSize2 / 2.0 + binSize2 * (idx % bins2);
    PRECISION probDens = histo[idx] / (numFrames * binSize1 * binSize2 * blen *
                                       blen * SIN(theta));
    if (probDens > 0) {
      atomicAdd(plnpsum, blen * blen * SIN(theta) * probDens * LOG(probDens) *
                             binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}

__global__ void cu_bdEnt(unsigned int *histo, const size_t numFrames,
                         const size_t bins1, const size_t bins2, PRECISION binSize1,
                         PRECISION binSize2, PRECISION min1, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t offset = blockDim.x * gridDim.x;
  while (idx < bins1 * bins2) {
    PRECISION blen = min1 + binSize1 / 2.0 + binSize1 * (idx / bins2);
    PRECISION probDens =
        histo[idx] / (numFrames * binSize1 * binSize2 * blen * blen);
    if (probDens > 0) {
      atomicAdd(plnpsum,
                blen * blen * probDens * LOG(probDens) * binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}

__global__ void cu_adEnt(unsigned int *histo, const size_t numFrames,
                         const size_t bins1, const size_t bins2, PRECISION binSize1,
                         PRECISION binSize2, PRECISION min1, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t offset = blockDim.x * gridDim.x;
  while (idx < bins1 * bins2) {
    PRECISION theta = min1 + binSize1 / 2.0 + binSize1 * (idx / bins2);
    PRECISION probDens =
        histo[idx] / (numFrames * binSize1 * binSize2 * SIN(theta));
    if (probDens > 0) {
      atomicAdd(plnpsum,
                SIN(theta) * probDens * LOG(probDens) * binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}

__global__ void cu_bbEnt(unsigned int *histo, const size_t numFrames,
                         const size_t bins, PRECISION binSize1, PRECISION binSize2,
                         PRECISION min1, PRECISION min2, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t offset = blockDim.x * gridDim.x;
  while (idx < bins * bins) {
    PRECISION blen1 = min1 + binSize1 / 2.0 + binSize1 * (idx / bins);
    PRECISION blen2 = min2 + binSize2 / 2.0 + binSize2 * (idx % bins);
    PRECISION probDens = histo[idx] / (numFrames * binSize1 * binSize2 * blen1 *
                                       blen1 * blen2 * blen2);
    if (probDens > 0) {
      atomicAdd(plnpsum, blen1 * blen1 * blen2 * blen2 * probDens *
                             LOG(probDens) * binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}

__global__ void cu_aaEnt(unsigned int *histo, const size_t numFrames,
                         const size_t bins, PRECISION binSize1, PRECISION binSize2,
                         PRECISION min1, PRECISION min2, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t offset = blockDim.x * gridDim.x;
  while (idx < bins * bins) {
    PRECISION theta1 = min1 + binSize1 / 2.0 + binSize1 * (idx / bins);
    PRECISION theta2 = min2 + binSize2 / 2.0 + binSize2 * (idx % bins);
    PRECISION probDens = histo[idx] / (numFrames * binSize1 * binSize2 *
                                       SIN(theta1) * SIN(theta2));
    if (probDens > 0) {
      atomicAdd(plnpsum, SIN(theta1) * SIN(theta2) * probDens * LOG(probDens) *
                             binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}

__global__ void cu_ddEnt(unsigned int *histo, const size_t numFrames,
                         const size_t bins, PRECISION binSize1, PRECISION binSize2,
                         PRECISION *plnpsum, unsigned int *occupbins) {

  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  size_t offset = blockDim.x * gridDim.x;
  while (idx < bins * bins) {
    PRECISION probDens = histo[idx] / (numFrames * binSize1 * binSize2);
    if (probDens > 0) {
      atomicAdd(plnpsum, probDens * LOG(probDens) * binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}


