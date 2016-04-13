/*
 * Copyright (c) 2012-2014 Mohamad Amer Wafai <amerwafai@gmail.com>
 * Mhd. Amer Wafai
 * 30.05.2012
 */

#ifndef CUDA_UTILS_CUH
#define CUDA_UTILS_CUH
#ifndef NDEBUG
#define cudaVerify(x) do {                                               \
        cudaError_t __cu_result = x;                                         \
        if (__cu_result != cudaSuccess) {                                      \
            fprintf(stderr, "%s:%i: error: cuda function call failed:\n"        \
                            "  %s;\nmessage: %s\n",                                    \
                    __FILE__, __LINE__, # x, cudaGetErrorString(__cu_result));     \
            exit(1);                                                           \
        }                                                                    \
} \
    while (0)
#define cudaVerifyKernel(x) do {                                         \
        x;                                                                   \
        cudaError_t __cu_result = cudaGetLastError();                        \
        if (__cu_result != cudaSuccess) {                                      \
            fprintf(stderr, "%s:%i: error: cuda function call failed:\n"        \
                            "  %s;\nmessage: %s\n",                                    \
                    __FILE__, __LINE__, # x, cudaGetErrorString(__cu_result));     \
            exit(1);                                                           \
        }                                                                    \
} \
    while (0)
#else
#define cudaVerify(x) do {                                               \
        x;                                                                   \
} \
    while (0)
#define cudaVerifyKernel(x) do {                                         \
        x;                                                                   \
} \
    while (0)
#endif

#ifdef __cplusplus
extern "C" {
#endif

void Dmalloc(void ** device, size_t size);
void HmallocP(void ** host, size_t size);
void HDcopy(void * device, void * host, size_t size);
void DHcopy(void * host, void * device, size_t size);
void HDcopyAsync(void * device, void * host, size_t size);
void DHcopyAsync(void * host, void * device, size_t size);
void DDcopy(void * dev2, void * dev1, size_t size);
void HHcopy(void * host2, void * host1, size_t size);
void Dfree(void * device);
void Hfree(void *);
void HDSync();

#ifdef __cplusplus
}
#endif


#endif /*CUDA_UTILS_CUH*/
