/*
 * Copyright (c) 2012-2014 Mohamad Amer Wafai <amerwafai@gmail.com>
 * Amer Wafai
 * 04.06.2012
 */
#include <stdio.h>
#include <stdlib.h>
#include "GPU/cuda-utils.cuh"

void Dmalloc(void ** device, size_t size){
    cudaVerify(cudaMalloc(device, size));
}

void HmallocP(void ** host, size_t size){
    cudaVerify(cudaMallocHost(host, size));
}

void HDcopy(void * device, void * host, size_t size){
    cudaVerify(cudaMemcpy(device, host, size, cudaMemcpyHostToDevice));
}

void DHcopy(void * host, void * device, size_t size){
    cudaVerify(cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost));
}

void HDcopyAsync(void * device, void * host, size_t size){
    cudaVerify(cudaMemcpyAsync(device, host, size, cudaMemcpyHostToDevice));
}

void DHcopyAsync(void * host, void * device, size_t size){
    cudaVerify(cudaMemcpyAsync(host, device, size, cudaMemcpyDeviceToHost));
}

void DDcopy(void * dev2, void * dev1, size_t size){
    cudaVerify(cudaMemcpy(dev2, dev1, size, cudaMemcpyDeviceToDevice));
}

void HHcopy(void * host2, void * host1, size_t size){
    cudaVerify(cudaMemcpy(host2, host1, size, cudaMemcpyHostToHost));
}

void Dfree(void * device){
    cudaVerify(cudaFree(device));
    device = NULL;
}

void Hfree(void * host){
    cudaVerify(cudaFreeHost(host));
    host = NULL;
}

void HDSync(){
    cudaDeviceSynchronize();
}
