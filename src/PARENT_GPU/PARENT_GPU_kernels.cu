

__global__ void cu_bondEnt(PRECISION *d_bondsEntropy1DMaster,
                           PRECISION *d_minBonds, PRECISION *d_maxBonds,
                           PRECISION *d_bondsChunk, const int nBonds,
                           const int bDens1D, const int numFrames) {

  int idx =
      blockIdx.x * blockDim.x + threadIdx.x; // points to the number of the bond
  if (idx > nBonds - 1)
    return;
  d_bondsEntropy1DMaster[idx] = 0;

  PRECISION binsize, probDens, blen, plnpsum;
  int occupbins;

  extern __shared__ int histo[];
  int offset_histo = threadIdx.x * bDens1D;

  for (int k = 0; k < bDens1D; k++) {
    histo[offset_histo + k] = 0; // initialize a histogram with zeros
  }

  binsize = (d_maxBonds[idx] - d_minBonds[idx]) /
            bDens1D; // and calculate the size of the bins
  for (int i = 0; i < numFrames;
       i++) { // and fill the histogram using all frames of the trajectory
    histo[offset_histo +
          int((d_bondsChunk[idx * numFrames + i] - d_minBonds[idx]) /
              binsize)] += 1;
  }

  occupbins = 0; // then use the histogram to calculate the (discretized)
                 // entropy, taking care of the Jacobian
  plnpsum = 0;
  blen = d_minBonds[idx] + (binsize / 2.0);
  for (int k = 0; k < bDens1D; k++) {
    probDens = histo[offset_histo + k] / (numFrames * binsize * blen * blen);
    if (probDens > 0) {
      plnpsum = plnpsum + blen * blen * probDens * log(probDens);
      occupbins = occupbins + 1;
    }
    blen += binsize;
  }
  plnpsum = -plnpsum * binsize;
  d_bondsEntropy1DMaster[idx] =
      plnpsum + (occupbins - 1.0) /
                    (2.0 * numFrames); // and apply Herzel entropy unbiasing
}

__global__ void cu_angleEnt(PRECISION *d_anglesEntropy1DMaster,
                            PRECISION *d_minAngles, PRECISION *d_maxAngles,
                            PRECISION *d_anglesChunk, const int nAngles,
                            const int aDens1D, const int numFrames) {

  int idx = blockIdx.x * blockDim.x +
            threadIdx.x; // points to the number of the angle
  if (idx > nAngles - 1)
    return;
  d_anglesEntropy1DMaster[idx] = 0;

  PRECISION binsize, probDens, theta, plnpsum;
  int occupbins;

  extern __shared__ int histo[];
  int offset_histo = threadIdx.x * aDens1D;

  for (int k = 0; k < aDens1D; k++) {
    histo[offset_histo + k] = 0; // initialize a histogram with zeros
  }

  binsize = (d_maxAngles[idx] - d_minAngles[idx]) /
            aDens1D; // and calculate the size of the bins
  for (int i = 0; i < numFrames;
       i++) { // and fill the histogram using all frames of the trajectory
    histo[offset_histo +
          int((d_anglesChunk[idx * numFrames + i] - d_minAngles[idx]) /
              binsize)] += 1;
  }

  occupbins = 0; // then use the histogram to calculate the (discretized)
                 // entropy, taking care of the Jacobian
  plnpsum = 0;
  theta = d_minAngles[idx] + (binsize / 2.0);
  for (int k = 0; k < aDens1D; k++) {
    probDens = histo[offset_histo + k] / (numFrames * binsize * sin(theta));
    if (probDens > 0) {
      plnpsum = plnpsum + sin(theta) * probDens * log(probDens);
      occupbins = occupbins + 1;
    }
    theta += binsize;
  }
  plnpsum = -plnpsum * binsize;
  d_anglesEntropy1DMaster[idx] =
      plnpsum + (occupbins - 1.0) /
                    (2.0 * numFrames); // and apply Herzel entropy unbiasing
}

__global__ void
cu_dihedralEnt(PRECISION *d_dihedralsEntropy1DMaster, PRECISION *d_minDihedrals,
               PRECISION *d_maxDihedrals, PRECISION *d_dihedralsChunk,
               const int nDihedrals, const int dDens1D, const int numFrames) {

  int idx = blockIdx.x * blockDim.x +
            threadIdx.x; // points to the number of the dihedral
  if (idx > nDihedrals - 1)
    return;
  d_dihedralsEntropy1DMaster[idx] = 0;

  PRECISION binsize, probDens, plnpsum;
  int occupbins;

  extern __shared__ int histo[];
  int offset_histo = threadIdx.x * dDens1D;

  for (int k = 0; k < dDens1D; k++) {
    histo[offset_histo + k] = 0; // initialize a histogram with zeros
  }

  binsize = (d_maxDihedrals[idx] - d_minDihedrals[idx]) /
            dDens1D; // and calculate the size of the bins
  for (int i = 0; i < numFrames;
       i++) { // and fill the histogram using all frames of the trajectory
    histo[offset_histo +
          int((d_dihedralsChunk[idx * numFrames + i] - d_minDihedrals[idx]) /
              binsize)] += 1;
  }

  occupbins = 0; // then use the histogram to calculate the (discretized)
                 // entropy, taking care of the Jacobian
  plnpsum = 0;
  for (int k = 0; k < dDens1D; k++) {
    probDens = histo[offset_histo + k] / (numFrames * binsize);
    if (probDens > 0) {
      plnpsum = plnpsum + probDens * log(probDens);
      occupbins = occupbins + 1;
    }
  }
  plnpsum = -plnpsum * binsize;
  d_dihedralsEntropy1DMaster[idx] =
      plnpsum + (occupbins - 1.0) /
                    (2.0 * numFrames); // and apply Herzel entropy unbiasing
}

__global__ void histo2D(PRECISION *traj1, PRECISION *traj2, const int numFrames,
                        unsigned int *histo, const int bins2,
                        PRECISION binSize1, PRECISION binSize2, PRECISION min1,
                        PRECISION min2) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int offset = blockDim.x * gridDim.x;
  //~ if(idx==0)printf("%d %p\n",idx,histo);
  while (idx < numFrames) {
    if (int((traj1[idx] - min1) / binSize1) < 0)
      printf("H1 %d %d %f %f %f %d\n", idx, int((traj1[idx] - min1) / binSize1),
             traj1[idx], min1, binSize1, numFrames);
    if (int((traj2[idx] - min2) / binSize2) < 0)
      printf("H2 %d %d %f %f %f %d\n", idx, int((traj2[idx] - min2) / binSize2),
             traj2[idx], min2, binSize2, numFrames);
    if (int((traj1[idx] - min1) / binSize1) > 49)
      printf("H3 %d %d %f %f %f %d\n", idx, int((traj1[idx] - min1) / binSize1),
             traj1[idx], min1, binSize1, numFrames);
    if (int((traj2[idx] - min2) / binSize2) > 49)
      printf("H4 %d %d %f %f %f %d\n", idx, int((traj2[idx] - min2) / binSize2),
             traj2[idx], min2, binSize2, numFrames);
    //~ atomicAdd(&histo[int((traj1[idx]-min1)/binSize1) * bins2 +
    //int((traj2[idx]-min2)/binSize2)], 1); ~
    //printf("%d\n",int((traj1[idx]-min1)/binSize1) * bins2 +
    //int((traj2[idx]-min2)/binSize2),int((traj1[idx]-min1)/binSize1) *
    //bins2,int((traj2[idx]-min2)/binSize2));
    atomicAdd(&histo[int((traj1[idx] - min1) / binSize1) * bins2 +
                     int((traj2[idx] - min2) / binSize2)],
              1);
    idx += offset;
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
                                       blen * sin(theta));
    if (probDens > 0) {
      atomicAdd(plnpsum, blen * blen * sin(theta) * probDens * log(probDens) *
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
                blen * blen * probDens * log(probDens) * binSize1 * binSize2);
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
        histo[idx] / (numFrames * binSize1 * binSize2 * sin(theta));
    if (probDens > 0) {
      atomicAdd(plnpsum,
                sin(theta) * probDens * log(probDens) * binSize1 * binSize2);
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
                             log(probDens) * binSize1 * binSize2);
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
                                       sin(theta1) * sin(theta2));
    if (probDens > 0) {
      atomicAdd(plnpsum, sin(theta1) * sin(theta2) * probDens * log(probDens) *
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
      atomicAdd(plnpsum, probDens * log(probDens) * binSize1 * binSize2);
      atomicAdd(occupbins, 1);
    }
    idx += offset;
  }
}
