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


__global__ void histo2D(PRECISION* __restrict__ traj1, PRECISION* __restrict__ traj2, const int numFrames,
                        unsigned int* __restrict__ histo, const int n_bins,
                        PRECISION binSize1, PRECISION binSize2, PRECISION min1,
                        PRECISION min2) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
    
  while (idx < numFrames) {
    atomicAdd(&histo[int((traj1[idx] - min1) / binSize1) * n_bins +
                     int((traj2[idx] - min2) / binSize2)],
              1);
    idx += offset;
  }
}

__global__ void histo2D_shared_block(PRECISION* __restrict__ traj1, PRECISION* __restrict__ traj2, const int numFrames,
                        unsigned int* __restrict__ histo, const int n_bins,
                        PRECISION binSize1, PRECISION binSize2, PRECISION min1,
                        PRECISION min2) {
                        
    extern __shared__ unsigned int histo_block[];

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
    int tid = threadIdx.x;
    int n_bins_total = n_bins * n_bins;
                        
    while(tid < n_bins_total){
        histo_block[tid] = 0;
        tid += blockDim.x;
    }
    
  while (idx < numFrames) {
    atomicAdd(&histo_block[int((traj1[idx] - min1) / binSize1) * n_bins +
                     int((traj2[idx] - min2) / binSize2)],
              1);
    idx += offset;
  }
  
  
  __syncthreads();
  
  tid = threadIdx.x;
  while(tid < n_bins_total){
    atomicAdd(&(histo[tid]), histo_block[tid]);
    tid += blockDim.x;
  }
  
}

__global__ void cu_baEnt(unsigned int *histo, const int numFrames,
                         const int bins1, const int bins2, PRECISION binSize1,
                         PRECISION binSize2, PRECISION min1, PRECISION min2,
                         PRECISION *plnpsum, unsigned int *occupbins) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
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

__global__ void cu_bdEnt(unsigned int *histo, const int numFrames,
                         const int bins1, const int bins2, PRECISION binSize1,
                         PRECISION binSize2, PRECISION min1, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
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

__global__ void cu_adEnt(unsigned int *histo, const int numFrames,
                         const int bins1, const int bins2, PRECISION binSize1,
                         PRECISION binSize2, PRECISION min1, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
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

__global__ void cu_bbEnt(unsigned int *histo, const int numFrames,
                         const int bins, PRECISION binSize1, PRECISION binSize2,
                         PRECISION min1, PRECISION min2, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
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

__global__ void cu_aaEnt(unsigned int *histo, const int numFrames,
                         const int bins, PRECISION binSize1, PRECISION binSize2,
                         PRECISION min1, PRECISION min2, PRECISION *plnpsum,
                         unsigned int *occupbins) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
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

__global__ void cu_ddEnt(unsigned int *histo, const int numFrames,
                         const int bins, PRECISION binSize1, PRECISION binSize2,
                         PRECISION *plnpsum, unsigned int *occupbins) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
  while (idx < bins * bins) {
    PRECISION probDens = histo[idx] / (numFrames * binSize1 * binSize2);
    if (probDens > 0) {
      atomicAdd(plnpsum, probDens * LOG(probDens) * binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}


