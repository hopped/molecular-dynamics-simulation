/*
 * Copyright (c) 2012-2014 Mohamad Amer Wafai <amerwafai@gmail.com>
 * Amer Wafai
 * 07.08.2012
 */

#ifndef BASICN2_CUH
#define BASICN2_CUH

#include "GPU/cuda-utils.cuh"
#include "config.h"
#include "lj_parameter.h"

#define XP(mc, member) ( mc->member )
#define YP(mc, member, i) ( mc->member + i )
#define ZP(mc, member, i) ( mc->member + 2 * i )

#ifdef __cplusplus
extern "C"
{
#endif

void cuBasicN2_create_streams(int N);
void HDcopyAsyncStreams(void * device, void * host, size_t size,int stream_id);
void DHcopyAsyncStreams(void * host, void * device, size_t size,int stream_id);
void cuBasicN2_destroy_streams(int N);
void cuBasicN2_set_DHN(long N);
void cuBasicN2_set_constants(long N);
void cuBasicN2_calc_forces(real *r,real *f,real *hr, real *U_pot, long *count,long *hcount, long N, long hN, int stream_id);
void cuBasicN2_reset_forces(real *f,long *count,long DN, int stream_id);
void cuBasicN2_Sync();
#ifdef __cplusplus
}
#endif

#endif /*BASICN2_CUH*/
