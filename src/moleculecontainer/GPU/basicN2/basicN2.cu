/*
 * Copyright (c) 2012-2014 Mohamad Amer Wafai <amerwafai@gmail.com>
 * Amer Wafai
 * 07.08.2012
 */
#include "moleculecontainer/GPU/basicN2/basicN2.cuh"
__global__ void basicN2_reset_forces_and_momenta_kernel(real *fX,real *fY,real *fZ, long DN);
__global__ void basicN2_calc_forces_kernel(real *rX,real *rY,real *rZ,real *fX,real *fY, real *fZ, real *hrX,real *hrY, real *hrZ,real *U_pot, long DN, long DHN);
inline __device__ double atomicAdd(double *address, double inc );

__constant__ real Dsigma2;
__constant__ real Depsilon24;
__constant__ real Dcutoff_radius_sq;

#ifdef USE_STREAMS
cudaStream_t *streams;
#endif

#define BSX 128

dim3 BlockSz(BSX, 1, 1);
dim3 GridSz(1, 1, 1);



__global__ void basicN2_calc_forces_kernel(real *rX,real *rY,real *rZ,real *fX,real *fY, real *fZ, real *hrX,real *hrY, real *hrZ,real *U_pot, long DN, long DHN){
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    extern __shared__ real sh_buf[];
    for (; i < DN; i += blockDim.x * gridDim.x) {
        real U_pot_tmp = 0.;
        real dr[3];
        real dr2;
        real lj6;
        real lj12;
        real lj12m6;
        real u_pot;
        real factor;
        real invdr2;
        sh_buf[threadIdx.x] = rX[i];
        sh_buf[threadIdx.x + blockDim.x] = rY[i];
        sh_buf[threadIdx.x + 2 * blockDim.x] = rZ[i];
        __syncthreads();
        /* Calculate all interactions between real molecules. */
        for (int j = i+1; j < DN; j++) {
            dr[0] = rX[j] - sh_buf[threadIdx.x ];
            dr[1] = rY[j] - sh_buf[threadIdx.x + blockDim.x];
            dr[2] = rZ[j] - sh_buf[threadIdx.x + 2 * blockDim.x];

            dr2 = 0.;
            dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

            if (dr2 > Dcutoff_radius_sq) {
                continue;
            }

            invdr2 = 1. / dr2;

            /* Lennard Jones interaction forces */
            lj6 = Dsigma2 * invdr2;
            lj6 = lj6 * lj6 * lj6;
            lj12 = lj6 * lj6;
            lj12m6 = lj12 - lj6;
            u_pot  = Depsilon24 * lj12m6;
            factor = Depsilon24 * ( lj12 + lj12m6 ) * invdr2;

            atomicAdd(&fX[i],-1* factor * dr[0]);
            atomicAdd(&fY[i],-1* factor * dr[1]);
            atomicAdd(&fZ[i],-1* factor * dr[2]);
            atomicAdd(&fX[j], factor * dr[0]);
            atomicAdd(&fY[j], factor * dr[1]);
            atomicAdd(&fZ[j], factor * dr[2]);
            U_pot_tmp += u_pot;
        }
    
        for (int j = 0; j < DHN; j++) {
            dr[0] = hrX[j] - sh_buf[threadIdx.x];
            dr[1] = hrY[j] - sh_buf[threadIdx.x + blockDim.x];
            dr[2] = hrZ[j] - sh_buf[threadIdx.x + 2 * blockDim.x];

            dr2 = 0.;
            dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if (dr2 > Dcutoff_radius_sq) {
                continue;
            }

            invdr2 = 1. / dr2;

            /* Lennard Jones interaction forces */

            lj6 = Dsigma2 * invdr2;
            lj6 = lj6 * lj6 * lj6;
            lj12 = lj6 * lj6;
            lj12m6 = lj12 - lj6;
            u_pot  = Depsilon24 * lj12m6;
            factor = Depsilon24 * ( lj12 + lj12m6 ) * invdr2;

            atomicAdd(&fX[i],-1* factor * dr[0]);
            atomicAdd(&fY[i],-1* factor * dr[1]);
            atomicAdd(&fZ[i],-1* factor * dr[2]);
            U_pot_tmp += u_pot * .5;
        }
        U_pot_tmp/=6.0;
        atomicAdd(U_pot, U_pot_tmp);
    }
}

__global__ void basicN2_reset_forces_and_momenta_kernel(real *fX,real *fY,real *fZ, long DN){
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < DN; i += blockDim.x * gridDim.x) {
        fX[i] = 0.;
        fY[i] = 0.;
        fZ[i] = 0.;
    }
}
inline __device__ double atomicAdd(double *address, double value)
{
    unsigned long long oldval, newval, assumed;

    oldval = __double_as_longlong(*address);
    do {
        assumed = oldval;
        newval = __double_as_longlong(__longlong_as_double(oldval) + value);
        oldval = atomicCAS((unsigned long long *)address, assumed, newval);
    } while (assumed != oldval);
    return __longlong_as_double(oldval);
}


void cuBasicN2_set_constants(long N){
    BlockSz.x = BSX;
    GridSz.x = ( N + BlockSz.x - 1 ) / BlockSz.x;
    cudaVerify(cudaMemcpyToSymbolAsync(Dsigma2, &sigma2, sizeof( real ),0,cudaMemcpyHostToDevice));
    cudaVerify(cudaMemcpyToSymbolAsync(Depsilon24, &epsilon24, sizeof( real ),0, cudaMemcpyHostToDevice));
    cudaVerify(cudaMemcpyToSymbolAsync(Dcutoff_radius_sq, &config.cutoff_radius_sq, sizeof( real ),0, cudaMemcpyHostToDevice));
}

void cuBasicN2_copyin(real *rX,real *rY,real *rZ,real *hrX,real *hrY,real *hrZ,real *drX,real *drY,real *drZ,real *dhrX,real *dhrY,real *dhrZ,long DN,long DHN, int stream_id){
#ifdef USE_STREAMS
	cudaVerify(cudaMemcpyAsync(drX,rX,DN*sizeof(real),cudaMemcpyHostToDevice, streams[stream_id]));
	cudaVerify(cudaMemcpyAsync(drY,rY,DN*sizeof(real),cudaMemcpyHostToDevice, streams[stream_id]));
	cudaVerify(cudaMemcpyAsync(drZ,rZ,DN*sizeof(real),cudaMemcpyHostToDevice, streams[stream_id]));
	cudaVerify(cudaMemcpyAsync(dhrX,hrX,DHN*sizeof(real),cudaMemcpyHostToDevice, streams[stream_id]));
	cudaVerify(cudaMemcpyAsync(dhrY,hrY,DHN*sizeof(real),cudaMemcpyHostToDevice, streams[stream_id]));
	cudaVerify(cudaMemcpyAsync(dhrZ,hrZ,DHN*sizeof(real),cudaMemcpyHostToDevice, streams[stream_id]));
#else
	cudaVerify(cudaMemcpyAsync(drX,rX,DN*sizeof(real),cudaMemcpyHostToDevice));
	cudaVerify(cudaMemcpyAsync(drY,rY,DN*sizeof(real),cudaMemcpyHostToDevice));
	cudaVerify(cudaMemcpyAsync(drZ,rZ,DN*sizeof(real),cudaMemcpyHostToDevice));
	cudaVerify(cudaMemcpyAsync(dhrX,hrX,DHN*sizeof(real),cudaMemcpyHostToDevice));
	cudaVerify(cudaMemcpyAsync(dhrY,hrY,DHN*sizeof(real),cudaMemcpyHostToDevice));
	cudaVerify(cudaMemcpyAsync(dhrZ,hrZ,DHN*sizeof(real),cudaMemcpyHostToDevice));
#endif
}
void cuBasicN2_reset_forces(real *f,long *count,long N, int stream_id){
#ifdef USE_STREAMS
    cudaVerifyKernel((basicN2_reset_forces_and_momenta_kernel<<<GridSz, BlockSz,0,streams[stream_id]>>>( &f[stream_id*3*N],&f[stream_id*3*N+N],&f[stream_id*3*N+N*2], count[stream_id])));
#else
    cudaVerifyKernel((basicN2_reset_forces_and_momenta_kernel<<<GridSz, BlockSz>>>( &f[stream_id*3*N], &f[stream_id*3*N+N], &f[stream_id*3*N+N*2], count[stream_id])));
#endif
}
void cuBasicN2_calc_forces(real *r,real *f,real *hr, real *U_pot, long *count,long *hcount,long N, long hN, int stream_id){
#ifdef USE_STREAMS
//    cudaVerifyKernel((basicN2_reset_forces_and_momenta_kernel<<<GridSz, BlockSz,0,streams[stream_id]>>>( &fX[stream_id*N],&fY[stream_id*N],&fZ[stream_id*N], count[stream_id])));
    cudaVerifyKernel(( basicN2_calc_forces_kernel << < GridSz, BlockSz,3*BlockSz.x*sizeof(real), streams[stream_id] >> >( &r [stream_id*3* N], &r [stream_id*3*N + N], &r [stream_id*3*N + N*2], 
															  &f [stream_id*3* N], &f [stream_id*3*N + N], &f [stream_id*3*N + N*2], 
                                                                                                                          &hr[stream_id*3*hN], &hr[stream_id*3*hN+hN], &hr[stream_id*3*hN+hN*2], 
															  &U_pot[stream_id], count[stream_id], hcount[stream_id] )));
#else
//    cudaVerifyKernel((basicN2_reset_forces_and_momenta_kernel<<<GridSz, BlockSz>>>( &fX[stream_id*N], &fY[stream_id*N], &fZ[stream_id*N], count[stream_id])));
    cudaVerifyKernel(( basicN2_calc_forces_kernel << < GridSz, BlockSz,3*BlockSz.x*sizeof(real) >> >( &r [stream_id*3*N ], &r [stream_id*3*N +N], &r [stream_id*3*N +N*2], 
												      &f [stream_id*3*N ], &f [stream_id*3*N +N], &f [stream_id*3*N +N*2], 
                                                                                                      &hr[stream_id*3*hN], &hr[stream_id*3*hN+N], &hr[stream_id*3*hN+N*2], 
												      &U_pot[stream_id], count[stream_id], hcount[stream_id] )));
#endif
}
void cuBasicN2_Sync(){
    cudaDeviceSynchronize();
}
#ifdef USE_STREAMS
void cuBasicN2_create_streams(int N){
    if (( streams = (cudaStream_t *)malloc(N * sizeof( cudaStream_t ))) == NULL) {
        printf("Sorry there is no enough memory for stream's array");
    }
    for (int i = 0; i < N; i++) {
        cudaStreamCreate(&streams[i]);
    }
}

void cuBasicN2_destroy_streams(int N){
    for (int i = 0; i < N; i++) {
        cudaStreamDestroy(streams[i]);
    }
}

void HDcopyAsyncStreams(void * device, void * host, size_t size,int stream_id){
    cudaVerify(cudaMemcpyAsync(device, host, size, cudaMemcpyHostToDevice,streams[stream_id]));
}

void DHcopyAsyncStreams(void * host, void * device, size_t size, int stream_id){
    cudaVerify(cudaMemcpyAsync(host, device, size, cudaMemcpyDeviceToHost,streams[stream_id]));
}

#endif

