/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 * Copyright (c) 2014      Mohamad Amer Wafai <amerwafai@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#define BASICN2_C
#include "basicN2.h"
#include "lj_parameter.h"
#include "utils/logger.h"
#include "constants.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>


#ifdef MPI
#include <string.h>
#include <mpi.h>
MPI_Datatype  MPI_MOL_T;
#endif



#ifdef GPU
#include "moleculecontainer/GPU/basicN2/basicN2.cuh"
static real * dr;
static real * dhr;
static real * df;
static real * dU_pot;
static long * dcount;
static long * dhcount;
#else
#ifdef _OPENACC
#include <openacc.h>
long prev_halo_count = 0;
int  not_first       = 0; 
long Num             = 0;
#endif
#endif

static real * sr;
static real * hsr;
static real * sf;
static long * indx;
static long * count;
static long * hcount;
static real * Upot_buffer;


int nfactor = 9000; 
int hfactor = 9000; 
int halo_count = 0; 
int SSx = 1; 
int SSy = 1; 
int SSz = 1; 
int NSSD = 1; 
int SizeOfSubDom = 200; 
int CPU_PART= 0;

static void set_subdomains(long N){
    long temp = ( N + SizeOfSubDom - 1 ) / SizeOfSubDom;
    SSx = ceil(pow(temp, 1. / 3)); //adding small number to prevent rounding errors
    SSz = sqrt((temp+SSx-1)/SSx);
    long temp1 = SSz * SSx;
    SSy = ceil(temp / temp1);
    NSSD = SSx * SSy * SSz;
    //printf("SSx=%d, SSy=%d, SSz=%d\n",SSx,SSy,SSz);
}


void basicN2_init(void *container, ensemble_t *ensemble , domain_t * domain) {
    basicN2 *mc = (basicN2 *) container;
    long N = ensemble->N;
    mc->capacity = N;
    mc->N  = 0;
    if(ensemble->NPS > 0) 
      SizeOfSubDom = ensemble->NPS; 
    else 
      ensemble->NPS = SizeOfSubDom;
    set_subdomains(N);
    hsr = NULL;
    sr = NULL;
    sf = NULL;
    indx = NULL;
    count = NULL;
    hcount = NULL;
    while((domain->L[0]/(real)SSx)<mc->cutoff_radius||(domain->L[1]/(real)SSy)<mc->cutoff_radius||(domain->L[2]/(real)SSz)<mc->cutoff_radius){
	  SizeOfSubDom++;
          set_subdomains(N);
          LOG(Debug,"Changing NPS to %lu\n",SizeOfSubDom);
    }
    ensemble->NPS = SizeOfSubDom;
#ifdef MPI
    const int tcount = 4;
    int i, array_of_lens[4] = {1, 3, 3, 4};
    mpi_buffer_t *temp = (mpi_buffer_t*)malloc(sizeof(mpi_buffer_t));
     /* compute displacements of structure components */
    MPI_Aint disp[4] = {0, 0, 0, 0};
    MPI_Address(&temp[0].id, disp  );
    MPI_Address( temp[0].r , disp+1);
    MPI_Address( temp[0].v , disp+2);
    MPI_Address( temp[0].t , disp+3);
    int base = disp[0];
    for (i = 0; i < 4; i++) disp[i] -= base;
 
    MPI_Datatype array_of_types[4] = {MPI_LONG, MPIREAL, MPIREAL, MPI_INT};
    MPI_Type_create_struct(tcount, array_of_lens, disp, array_of_types, &MPI_MOL_T);
    MPI_Type_commit(&MPI_MOL_T);
    free(temp);
    
    mc->send_buffer = NULL;
    mc->recv_buffer = NULL;
    mc->Max_transfer_sz = 0;
    mc->Max_receive_sz = 0;
    mc->t = (int *) realloc(mc->t, 3 * N * sizeof( int ));
    mc->Max_receive_sz = 0;
    mc->Max_receive_sz = 0;
#endif

#ifdef GPU
    if(NSSD > 1){
      if(ensemble->CPU_PART <=  NSSD)
         CPU_PART = ensemble->CPU_PART;
      else{
         CPU_PART = NSSD;
         ensemble->CPU_PART = NSSD;
      }
    }else{
      if(ensemble->CPU_PART >  NSSD)
         ensemble->CPU_PART = NSSD;
         CPU_PART = ensemble->CPU_PART;
    }
#ifdef USE_STREAMS
    cuBasicN2_create_streams(NSSD);
#endif
    Hfree(mc->r);
    HmallocP((void **)&(mc->r), 3 * N * sizeof(real));
    Hfree(mc->F);
    HmallocP((void **)&(mc->F), 3 * N * sizeof(real));
    Dmalloc((void **)&dU_pot, (NSSD - CPU_PART) * sizeof( real ));
    HmallocP((void **)&Upot_buffer, NSSD * sizeof(real));
    cuBasicN2_set_constants(ensemble->NPS);
#else
    ensemble->CPU_PART = NSSD;
    mc->r  = (real *) realloc(mc->r, 3 * N * sizeof(real));
    mc->F  = (real *) realloc(mc->F, 3 * N * sizeof(real));
    Upot_buffer  = (real *) realloc(Upot_buffer, NSSD * sizeof(real));
#endif

    mc->v  = (real *) realloc(mc->v, 3*N * sizeof(real));
    mc->id = (long *) realloc(mc->id, N * sizeof(long));
    mc->halo = (basicN2 *) basicN2_alloc();
    mc->cutoff_radius = config.cutoff_radius;
    count = (long *) realloc(count, NSSD * sizeof( long ));
    hcount = (long *) realloc(hcount, NSSD * sizeof( long ));
    if (NSSD > 1) {
        indx  = (long *) realloc(indx, NSSD * (( N / NSSD ) + nfactor) * sizeof( long ));
#ifdef GPU 
        HmallocP((void **)&sr, 3 * NSSD * (( N / NSSD ) + nfactor) * sizeof(real)); 
        HmallocP((void **)&sf, 3 * NSSD * (( N / NSSD ) + nfactor) * sizeof(real)); 
#else 
        sr  = (real *) realloc(sr, 3 * NSSD * (( N / NSSD ) + nfactor) * sizeof( real )); 
        sf  = (real *) realloc(sf, 3 * NSSD * (( N / NSSD ) + nfactor) * sizeof( real )); 
#endif 
    } else { 
        sr = &X(mc, r, 0);
        sf = &X(mc, F, 0);
        count[0] = N;
        nfactor = 0;
        hfactor = 0;
    }
#ifdef GPU
    Dmalloc((void **)&dr, 3 * (NSSD - CPU_PART) * (( N / NSSD ) + nfactor) * sizeof( real ));
    Dmalloc((void **)&df, 3 * (NSSD - CPU_PART) * (( N / NSSD ) + nfactor) * sizeof( real ));
#else
#ifdef _OPENACC
    acc_create(sr, 3 * NSSD * (( N / NSSD ) + nfactor) * sizeof( real ));
    acc_create(sf, 3 * NSSD * (( N / NSSD ) + nfactor) * sizeof( real ));
    acc_create(Upot_buffer,  NSSD * sizeof( real ));
#endif
#endif
}

static void* basicN2_alloc_N(long N) {
    basicN2 *mc;
    mc = (basicN2 *) malloc(sizeof(basicN2));
    mc->capacity = N;
    mc->N  = 0;
    if(N > 0) {
#ifdef GPU 
        HmallocP((void **)&(mc->r), 3 * N * sizeof(real)); 
        HmallocP((void **)&(mc->F), 3 * N * sizeof(real)); 
#else 
        mc->r  = (real *) malloc(3 * N * sizeof(real)); 
        mc->F  = (real *) malloc(3 * N * sizeof(real)); 
#endif 
        mc->v  = (real *) malloc(3*N * sizeof(real));
        mc->id = (long *) malloc(N * sizeof(long));
#ifdef MPI
        mc->t = (int *) malloc(3 * N * sizeof( int ));
        mc->send_buffer = NULL;
        mc->recv_buffer = NULL;
        mc->Max_transfer_sz = 0;
        mc->Max_receive_sz = 0;
#endif
    } else {
#ifdef MPI
        mc->Max_transfer_sz = 0;
        mc->Max_receive_sz = 0;
        mc->t = NULL;
        mc->send_buffer = NULL;
        mc->recv_buffer = NULL;
#endif
        mc->r  = NULL;
        mc->v  = NULL;
        mc->F  = NULL;
        mc->id = NULL;
    }
    mc->cutoff_radius = config.cutoff_radius;
    mc->halo = NULL;
    return mc;
}

void* basicN2_alloc() {
    basicN2 *mc;
    mc = (basicN2 *) basicN2_alloc_N(0);
    return mc;
}

void basicN2_free(void *container) {
    basicN2 *mc = (basicN2 *) container;
    assert(mc != NULL);
    if(mc->halo) {
        basicN2_free(mc->halo);
        mc->halo = NULL;
    }
    mc->capacity = 0;
    mc->N = 0;
#ifdef GPU
    Hfree(mc->r); mc->r  = NULL;
    Hfree(mc->F); mc->F  = NULL;
#else
    free(mc->r);  mc->r  = NULL;
    free(mc->F);  mc->F  = NULL;
#endif
    free(mc->v);  mc->v  = NULL;
    free(mc->id); mc->id = NULL;
#ifdef MPI
    mc->Max_receive_sz = 0;
    mc->Max_receive_sz = 0;
    free (mc->send_buffer);
    mc->send_buffer = NULL;
    free (mc->recv_buffer);
    mc->recv_buffer = NULL;
    free (mc->t);
    mc->t = NULL;
#endif
    free(mc);
}

inline static long fetch_and_add_ld(long *var, long value) {
/* @todo: Implement as atomic operation. */
    long old_val = *var;
    *var += value;
    return old_val;
}

void basicN2_add_molecule(void *container, molecule_t *m) {
    basicN2 *mc = (basicN2 *) container;
    /* get insert position and increment insert position to next position. */
    long i = fetch_and_add_ld(&mc->N, 1); 
    X(mc, r, i) = m->r[0];
    Y(mc, r, i) = m->r[1];
    Z(mc, r, i) = m->r[2];
    X(mc, v, i) = m->v[0];
    Y(mc, v, i) = m->v[1];
    Z(mc, v, i) = m->v[2];
    LOG(Debug, "%lu\n", m->id);
    ID(mc, i)   = m->id;
}

inline static void copy_molecule(basicN2 *fromMc, long fromId, basicN2 *toMc, long toId) {
    assert( fromId < fromMc->capacity);
    assert( toId < toMc->capacity);
    X(toMc, r, toId) = X(fromMc, r, fromId);
    Y(toMc, r, toId) = Y(fromMc, r, fromId);
    Z(toMc, r, toId) = Z(fromMc, r, fromId);
    X(toMc, v, toId) = X(fromMc, v, fromId);
    Y(toMc, v, toId) = Y(fromMc, v, fromId);
    Z(toMc, v, toId) = Z(fromMc, v, fromId);
    X(toMc, F, toId) = X(fromMc, F, fromId);
    Y(toMc, F, toId) = Y(fromMc, F, fromId);
    Z(toMc, F, toId) = Z(fromMc, F, fromId);
    ID(toMc, toId)   = ID(fromMc, fromId);
#ifdef MPI
    X(toMc, t, toId) = X(fromMc, t, fromId);
    Y(toMc, t, toId) = Y(fromMc, t, fromId);
    Z(toMc, t, toId) = Z(fromMc, t, fromId);
#endif
}
static void basicN2_del_molecule(basicN2 *mc, long index) {
    long i = fetch_and_add_ld(&mc->N, -1) - 1; 
    copy_molecule(mc, i, mc, index);
    X (mc, r, i) = 0;
    Y (mc, r, i) = 0;
    Z (mc, r, i) = 0;
    X (mc, v, i) = 0;
    Y (mc, v, i) = 0;
    Z (mc, v, i) = 0;
    X (mc, F, i) = 0;
    Y (mc, F, i) = 0;
    Z (mc, F, i) = 0;
    ID(mc,    i) = 0;
#ifdef MPI
    X(mc, t, i) = 0;
    Y(mc, t, i) = 0;
    Z(mc, t, i) = 0;
#endif
}

void basicN2_del_molecule_by_id(void *container, long id) {
    basicN2 *mc = (basicN2 *) container;
    long i;
    for(i = 0; i < mc->N - 1; i++) {
        if( mc->id[i] == id)  break; 
    }
    basicN2_del_molecule(mc, i);
}
static void calc_forces(long N, long hN, int dom){
    long i, j;
    real U_pot_tmp = 0.;
    int size = count[dom];
    int hsize= hcount[dom];

    LOG(Info, "LJ parameters: epsilon = %" PRIreal ", sigma^2 = %" PRIreal "\n",
        epsilon24 / 24.,
        sigma2);
    /* Calculate all interactions between real molecules. */
#ifdef GPU
#pragma omp parallel default(shared) private(i,j) reduction(+:U_pot_tmp) if(NSSD==1)
{
#pragma omp for 
#else
#ifdef _OPENACC
#pragma acc parallel reduction(+:U_pot_tmp) present(sr[dom * 3 * N : 3 * N], sf[dom * 3 * N : 3 * N], hsr[dom * 3 * hN : 3 * hN], Upot_buffer[0:NSSD]) vector_length(256) async(dom) wait
{
#pragma acc loop reduction(+:U_pot_tmp) gang
#else
#pragma omp parallel default(shared) private(i,j) reduction(+:U_pot_tmp) if(NSSD==1)
{
#pragma omp for 
#endif
#endif
    for (i = 0; i < size; i++) {
        for (j = i + 1; j < size; j++) {
            real dr[3];
            dr[0] = sr[dom * 3 * N + j         ] - sr[dom * 3 * N + i        ];
            dr[1] = sr[dom * 3 * N + j + N     ] - sr[dom * 3 * N + i + N    ];
            dr[2] = sr[dom * 3 * N + j + N * 2 ] - sr[dom * 3 * N + i + N * 2];

            real dr2 = 0.;
            dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

            if (dr2 > config.cutoff_radius_sq) {
                continue;
            }
#ifndef GPU
#ifndef _OPENACC
#endif
            assert(dr2 != 0.0);
            assert(dr2 > 0.0);
#ifndef GPU
#endif
#endif

            real invdr2 = 1. / dr2;

            /* Lennard Jones interaction forces */
            real lj6 = sigma2 * invdr2;
            lj6 = lj6 * lj6 * lj6;
            real lj12 = lj6 * lj6;
            real lj12m6 = lj12 - lj6;
            real u_pot = epsilon24 * lj12m6;
            real factor = epsilon24 * ( lj12 + lj12m6 ) * invdr2;
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif
            sf[dom * 3 * N + j        ] += factor * dr[0];
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + j + N    ] += factor * dr[1];
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + j + N * 2] += factor * dr[2];
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + i        ] -= factor * dr[0];
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + i + N    ] -= factor * dr[1];
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + i + N * 2] -= factor * dr[2];
            U_pot_tmp += u_pot;
        }
    }

    /* Calculate all interactions of real with halo molecules. */
    /* Loops reordered as the i loop is longer. */
#ifdef GPU
#pragma omp for 
#else
#ifdef _OPENACC 
#pragma acc loop reduction(+:U_pot_tmp) gang
#else
#pragma omp for 
#endif
#endif 
    for (i = 0; i < size; i++) {
        for (j = 0; j < hsize; j++) {
            
            real dr[3];
            dr[0] = hsr[dom * 3 * hN + j         ] - sr[dom * 3 * N + i        ];
            dr[1] = hsr[dom * 3 * hN + j + hN    ] - sr[dom * 3 * N + i + N    ];
            dr[2] = hsr[dom * 3 * hN + j + hN * 2] - sr[dom * 3 * N + i + N * 2];
            real dr2 = 0.;
            dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

            if (dr2 > config.cutoff_radius_sq) {
                continue;
            }
#ifndef GPU
#ifndef _OPENACC
#endif
            assert(dr2 != 0.0);
            assert(dr2 > 0.0);
#ifndef GPU
#endif
#endif
            real invdr2 = 1. / dr2;

            /* Lennard Jones interaction forces */
            real lj6 = sigma2 * invdr2;
            lj6 = lj6 * lj6 * lj6;
            real lj12 = lj6 * lj6;
            real lj12m6 = lj12 - lj6;
            real u_pot = epsilon24 * lj12m6;
            real factor = epsilon24 * ( lj12 + lj12m6 ) * invdr2;

#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + i        ] -= factor * dr[0];
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + i + N    ] -= factor * dr[1];
#ifdef GPU
#pragma omp atomic update
#else
#ifdef _OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif 
#endif 
            sf[dom * 3 * N + i + N * 2] -= factor * dr[2];
            U_pot_tmp += u_pot * .5;
        }
    }
#ifdef GPU
}
#else
#ifndef _OPENACC
}
#endif
#endif


//#pragma acc wait(dom)
    U_pot_tmp /= 6.0;  /* newtown 3, using epsilon24 instead of epsilon4 in potential calculation */
    Upot_buffer[dom] = U_pot_tmp;
#ifndef GPU
#ifdef _OPENACC
}
#endif
#endif

}
void basicN2_copyBack(void * container){
    basicN2 * mc = (basicN2 *)container;
    long i,j;
    long N = ( mc->capacity / NSSD ) + nfactor;
#ifdef GPU
  cuBasicN2_Sync();
#else
#ifdef _OPENACC
#pragma acc wait
#endif
#endif
#pragma omp parallel default(shared) private(i,j) if(NSSD>1) 
#pragma omp for schedule(runtime)
    for(j=0; j<NSSD; j++){
       for (i = 0; i < count[j]; i++) {
           X(mc, F, indx[j * N + i]) = sf[j * 3 * N + i        ];
           Y(mc, F, indx[j * N + i]) = sf[j * 3 * N + i + N    ];
           Z(mc, F, indx[j * N + i]) = sf[j * 3 * N + i + N * 2];
       }
    }
}

void basicN2_calc_forces(void *container, real *U_pot) {
    basicN2 *mc = (basicN2 *) container;
    long N = ( mc->capacity / NSSD ) + nfactor;
    long hN = halo_count + hfactor;
    *U_pot = 0;
    int i;
#ifdef GPU
    cuBasicN2_Sync();
    for (i = 0; i < NSSD - CPU_PART; i++) {
        cuBasicN2_calc_forces(dr, df, dhr, dU_pot, count, hcount, N, hN, i);
        DHcopyAsyncStreams(&Upot_buffer[i],&dU_pot[i], sizeof(real), i);
    }
    //DHcopy(Upot_buffer,dU_pot, (NSSD - CPU_PART) * sizeof(real));
    DHcopyAsync(sf,df, (NSSD - CPU_PART) * 3 * N  * sizeof(real));
#pragma omp parallel if (CPU_PART>1)
#pragma omp for schedule(runtime)
    for (i = NSSD - CPU_PART; i < NSSD; i++) {
#pragma omp task 
{
	Upot_buffer[i]=0;
        calc_forces(N , hN, i);
}
    }
    cuBasicN2_Sync();
#else
#ifndef _OPENACC
#pragma omp parallel if (NSSD>1)
#pragma omp for schedule(runtime)
#endif
    for (i = 0; i < NSSD; i++) {
#ifndef _OPENACC
#pragma omp task 
{
#endif
	Upot_buffer[i]=0;
        calc_forces(N , hN, i);
#ifndef _OPENACC
}
#endif
   }
#ifdef _OPENACC
#pragma acc wait 
#pragma acc update host(Upot_buffer[0:NSSD])
#pragma acc update host(sf[0: 3 * N * NSSD]) async
#endif
#endif

    for (i = 0; i < NSSD; i++) {
	*U_pot+=Upot_buffer[i];
    }

}

void basicN2_clear(void *container) {
    basicN2 *mc = (basicN2 *) container;
    mc->N = 0;
    /* @todo: Do we need a memset 0 here? */
}

void basicN2_integrate_pref(void *container, simulation_t *simulation) {
    basicN2 *mc = (basicN2 *) container;
    long i;
    real dt = simulation->dt;
    assert(dt > 0.);
    real m  = config.m;
    assert(m > 0.);
    real a = 0.5 * dt / m;
#pragma omp parallel default(shared) private(i) if (NSSD>1)
{
    /* first sumand of v(t+dt) = v(t) + (F(t) + F(t+dt))/2m (old forces) */
#pragma omp for schedule(runtime)
    for(i = 0; i < mc->N; i++) {
        X(mc, v, i) += a * X(mc, F, i);
        Y(mc, v, i) += a * Y(mc, F, i);
        Z(mc, v, i) += a * Z(mc, F, i);
    }

#pragma omp for schedule(runtime)
    for(i = 0; i < mc->N; i++) {
        X(mc, r, i) += dt * X(mc, v, i);
        Y(mc, r, i) += dt * Y(mc, v, i);
        Z(mc, r, i) += dt * Z(mc, v, i);
    }
}

}
void basicN2_integrate_postf(void *container, simulation_t *simulation) {
    basicN2 *mc = (basicN2 *) container;
    long i;
    real dt = simulation->dt;
    real m  = config.m;
    real a = 0.5 * dt / m;
#ifdef GPU
  cuBasicN2_Sync();
#else
#ifdef _OPENACC
#pragma acc wait
#endif
#endif
    if (NSSD > 1) {
       basicN2_copyBack(container);
    }
    /* second sumand of v(t+dt) = v(t) + (F(t) + F(t+dt))/2m (new forces) */
#pragma omp parallel default(shared) private(i)
#pragma omp for schedule(runtime)
    for(i = 0; i < mc->N; i++) {
        X(mc, v, i) += a * X(mc, F, i);
        Y(mc, v, i) += a * Y(mc, F, i);
        Z(mc, v, i) += a * Z(mc, F, i);
    }
}

#ifdef MPI
static inline long find_max(long *X, long *Y, long *Z){
    long temp = X[p];
    int i = 0;
    if (X[n] > temp) {temp = X[n];}
    for (i = 0; i < 2; i++) {
        if (Y[i] > temp) {temp = Y[i];}
    }
    for (i = 0; i < 2; i++) {
        if (Z[i] > temp) {temp = Z[i];}
    }
    return temp;
}

static void basicN2_realloc(basicN2 *mc, long newCap){
    basicN2 *mc_new = basicN2_alloc_N(newCap);
    memcpy(mc_new->r, mc->r, mc->N * sizeof( real ));
    memcpy(mc_new->r + newCap, mc->r + mc->capacity, mc->N * sizeof( real ));
    memcpy(mc_new->r + 2 * newCap, mc->r + 2 * mc->capacity, mc->N * sizeof( real ));
    memcpy(mc_new->v, mc->v, mc->N * sizeof( real ));
    memcpy(mc_new->v + newCap, mc->v + mc->capacity, mc->N * sizeof( real ));
    memcpy(mc_new->v + 2 * newCap, mc->v + 2 * mc->capacity, mc->N * sizeof( real ));
    memcpy(mc_new->F, mc->F, mc->N * sizeof( real ));
    memcpy(mc_new->F + newCap, mc->F + mc->capacity, mc->N * sizeof( real ));
    memcpy(mc_new->F + 2 * newCap, mc->F + 2 * mc->capacity, mc->N * sizeof( real ));
    memcpy(mc_new->id, mc->id, mc->N * sizeof( long ));
    memcpy(mc_new->id + newCap, mc->id + mc->capacity, mc->N * sizeof( long ));
    memcpy(mc_new->id + 2 * newCap, mc->id + 2 * mc->capacity, mc->N * sizeof( long ));
    memcpy(mc_new->t, mc->t, mc->N * sizeof( int ));
    memcpy(mc_new->t + newCap, mc->t + mc->capacity, mc->N * sizeof( int ));
    memcpy(mc_new->t + 2 * newCap, mc->t + 2 * mc->capacity, mc->N * sizeof( int ));
    mc_new->capacity = newCap;
    free(mc);
    mc = mc_new;
    if (NSSD > 1) {
        free(indx);
        free(sr);
        free(sf);
        indx= (long *)malloc(    NSSD * (( newCap / NSSD ) + nfactor) * sizeof( long ));
        sr  = (real *)malloc(3 * NSSD * (( newCap / NSSD ) + nfactor) * sizeof( real ));
        sf  = (real *)malloc(3 * NSSD * (( newCap / NSSD ) + nfactor) * sizeof( real ));
    } else {
        sr = &X(mc, r, 0);
        sf = &X(mc, F, 0);
    }
}
static inline int find_hops(real dist, real dim){
    return ceil(( fabs(dist) / dim ));
}

static void basicN2_recv_buffer_update(basicN2 * mc, long recv_buffer_sz){
    if (mc->Max_receive_sz < recv_buffer_sz) {
        mc->Max_receive_sz = recv_buffer_sz;
        mc->recv_buffer = realloc(mc->recv_buffer, mc->Max_receive_sz * sizeof(mpi_buffer_t));
    }
    long newCap = mc->N + recv_buffer_sz;
    if (mc->capacity < newCap) {
        basicN2_realloc(mc, newCap);
    }
}

static void basicN2_extract_buffer(basicN2 * mc, long N){
    long i;
    for (i = 0; i < N; i++) {
        X(mc, t, mc->N ) = mc->recv_buffer[i].t[0];
        Y(mc, t, mc->N ) = mc->recv_buffer[i].t[1];
        Z(mc, t, mc->N ) = mc->recv_buffer[i].t[2];
        X(mc, r, mc->N ) = mc->recv_buffer[i].r[0];
        Y(mc, r, mc->N ) = mc->recv_buffer[i].r[1];
        Z(mc, r, mc->N ) = mc->recv_buffer[i].r[2];
        X(mc, v, mc->N ) = mc->recv_buffer[i].v[0];
        Y(mc, v, mc->N ) = mc->recv_buffer[i].v[1];
        Z(mc, v, mc->N ) = mc->recv_buffer[i].v[2];
        ID(mc, mc->N++ ) = mc->recv_buffer[i].id;
    }
}

static void basicN2_buffer_exchange(basicN2 *mc, int myrank, int direct, int sign){
   MPI_Status status;
   switch (direct){
   case X:
       if (myrank % Sx == 0) {
           MPI_Send(mc->send_buffer, mc->TXs[sign], MPI_MOL_T, mc->Nx[sign], 0, MPI_COMM_WORLD);
           MPI_Recv(mc->recv_buffer, mc->TXr[(sign+1)%2], MPI_MOL_T, mc->Nx[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
       } else {
           MPI_Recv(mc->recv_buffer, mc->TXr[(sign+1)%2], MPI_MOL_T, mc->Nx[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
           MPI_Send(mc->send_buffer, mc->TXs[sign], MPI_MOL_T, mc->Nx[sign], 0, MPI_COMM_WORLD);
       }
       basicN2_extract_buffer(mc, mc->TXr[(sign+1)%2]);
       break;
   case Y:
       if (myrank % ( Sx * Sy ) < Sx) {
           MPI_Send(mc->send_buffer, mc->TYs[sign], MPI_MOL_T, mc->Ny[sign], 0, MPI_COMM_WORLD);
           MPI_Recv(mc->recv_buffer, mc->TYr[(sign+1)%2], MPI_MOL_T, mc->Ny[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
       } else {
           MPI_Recv(mc->recv_buffer, mc->TYr[(sign+1)%2], MPI_MOL_T, mc->Ny[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
           MPI_Send(mc->send_buffer, mc->TYs[sign], MPI_MOL_T, mc->Ny[sign], 0, MPI_COMM_WORLD);
       }
       basicN2_extract_buffer(mc, mc->TYr[(sign+1)%2]);
       break;
   case Z:
       if (myrank / ( Sx * Sy ) == 0) {
           MPI_Send(mc->send_buffer, mc->TZs[sign], MPI_MOL_T, mc->Nz[sign], 0, MPI_COMM_WORLD);
           MPI_Recv(mc->recv_buffer, mc->TZr[(sign+1)%2], MPI_MOL_T, mc->Nz[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
       } else {
           MPI_Recv(mc->recv_buffer, mc->TZr[(sign+1)%2], MPI_MOL_T, mc->Nz[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
           MPI_Send(mc->send_buffer, mc->TZs[sign], MPI_MOL_T, mc->Nz[sign], 0, MPI_COMM_WORLD);
       }
       basicN2_extract_buffer(mc, mc->TZr[(sign+1)%2]);
       break;
   }
}

inline static void exchange_buffer_sizes(basicN2 *mc , int myrank, int direct, int sign){ 
    MPI_Status status; 
    switch (direct){ 
    case X: 
        if (myrank % Sx == 0) { 
            MPI_Send(&mc->TXs[sign], 1, MPI_LONG, mc->Nx[sign], 0, MPI_COMM_WORLD); 
            MPI_Recv(&mc->TXr[(sign+1)%2], 1, MPI_LONG, mc->Nx[(sign+1)%2], 0, MPI_COMM_WORLD, &status); 
        } else { 
            MPI_Recv(&mc->TXr[(sign+1)%2], 1, MPI_LONG, mc->Nx[(sign+1)%2], 0, MPI_COMM_WORLD, &status); 
            MPI_Send(&mc->TXs[sign], 1, MPI_LONG, mc->Nx[sign], 0, MPI_COMM_WORLD); 
        } 
        break; 
    case Y: 
        if (myrank % ( Sx * Sy ) < Sx) { 
            MPI_Send(&mc->TYs[sign], 1, MPI_LONG, mc->Ny[sign], 0, MPI_COMM_WORLD);
            MPI_Recv(&mc->TYr[(sign+1)%2], 1, MPI_LONG, mc->Ny[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&mc->TYr[(sign+1)%2], 1, MPI_LONG, mc->Ny[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
            MPI_Send(&mc->TYs[sign], 1, MPI_LONG, mc->Ny[sign], 0, MPI_COMM_WORLD);
        }
        break;
    case Z:
        if (myrank / ( Sx * Sy ) == 0) {
            MPI_Send(&mc->TZs[sign], 1, MPI_LONG, mc->Nz[sign], 0, MPI_COMM_WORLD);
            MPI_Recv(&mc->TZr[(sign+1)%2], 1, MPI_LONG, mc->Nz[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&mc->TZr[(sign+1)%2], 1, MPI_LONG, mc->Nz[(sign+1)%2], 0, MPI_COMM_WORLD, &status);
            MPI_Send(&mc->TZs[sign], 1, MPI_LONG, mc->Nz[sign], 0, MPI_COMM_WORLD);
        }
        break;
    }
}

static void basicN2_send_buffer_update(basicN2 *mc, int myrank, int direct, int sign){
    long i, count = 0;
    int struct_size = 6 * sizeof(real) + 3 * sizeof(int) + sizeof(long); 
    int cond;
    switch (direct){
    case X:
        for (i = 0; i < mc->N; i++) {
            if (sign == p) cond = (X(mc, t, i) >= 1);
            else cond = (X(mc, t, i) < 0);
            if (cond) {
                assert(count <= mc->TXs[sign]);
                if (sign == p)mc->send_buffer[count  ].t[0] = --X(mc,t,i);
                else          mc->send_buffer[count  ].t[0] = ++X(mc,t,i);
                              mc->send_buffer[count  ].t[1] =   Y(mc,t,i);
                              mc->send_buffer[count  ].t[2] =   Z(mc,t,i);
                              mc->send_buffer[count  ].r[0] =   X(mc,r,i);
                              mc->send_buffer[count  ].r[1] =   Y(mc,r,i);
                              mc->send_buffer[count  ].r[2] =   Z(mc,r,i);
                              mc->send_buffer[count  ].v[0] =   X(mc,v,i);
                              mc->send_buffer[count  ].v[1] =   Y(mc,v,i);
                              mc->send_buffer[count  ].v[2] =   Z(mc,v,i);
		              mc->send_buffer[count++].id   =   ID(mc, i); 
                basicN2_del_molecule(mc, i--);  //decremented to test the last element which is copied
            }
        }
        exchange_buffer_sizes(mc , myrank, X, sign); 
        break;
    case Y:
        for (i = 0; i < mc->N; i++) {
            if (sign == p) cond = (Y(mc, t, i) >= 1);
            else cond = (Y(mc, t, i) < 0);
            if (cond) {
                assert(count <= mc->TYs[sign]);
                              mc->send_buffer[count  ].t[0] =   X(mc,t,i);
                if(sign == p) mc->send_buffer[count  ].t[1] = --Y(mc,t,i);
                else          mc->send_buffer[count  ].t[1] = ++Y(mc,t,i);
                              mc->send_buffer[count  ].t[2] =   Z(mc,t,i);
                              mc->send_buffer[count  ].r[0] =   X(mc,r,i);
                              mc->send_buffer[count  ].r[1] =   Y(mc,r,i);
                              mc->send_buffer[count  ].r[2] =   Z(mc,r,i);
                              mc->send_buffer[count  ].v[0] =   X(mc,v,i);
                              mc->send_buffer[count  ].v[1] =   Y(mc,v,i);
                              mc->send_buffer[count  ].v[2] =   Z(mc,v,i);
		              mc->send_buffer[count++].id   =   ID(mc, i); 
                basicN2_del_molecule(mc, i--);  //decremented to test the last element which is copied
            }
        }
        exchange_buffer_sizes(mc , myrank, Y, sign); 
        break;
    case Z:
        for (i = 0; i < mc->N; i++) {
            if (sign == p) cond = (Z(mc, t, i) >= 1);
            else cond = (Z(mc, t, i) < 0);
            if (cond) {
                assert(count <= mc->TZs[sign]);
                              mc->send_buffer[count  ].t[0] =   X(mc,t,i);
                              mc->send_buffer[count  ].t[1] =   Y(mc,t,i);
                if(sign == p) mc->send_buffer[count  ].t[2] = --Z(mc,t,i);
                else          mc->send_buffer[count  ].t[2] = ++Z(mc,t,i);
                              mc->send_buffer[count  ].r[0] =   X(mc,r,i);
                              mc->send_buffer[count  ].r[1] =   Y(mc,r,i);
                              mc->send_buffer[count  ].r[2] =   Z(mc,r,i);
                              mc->send_buffer[count  ].v[0] =   X(mc,v,i);
                              mc->send_buffer[count  ].v[1] =   Y(mc,v,i);
                              mc->send_buffer[count  ].v[2] =   Z(mc,v,i);
		              mc->send_buffer[count++].id   =   ID(mc, i); 
                basicN2_del_molecule(mc, i--);  //decremented to test the last element which is copied
            }
        }
        exchange_buffer_sizes(mc , myrank, Z, sign); 
        break;
    }
}

static void apply_periodic_boundary_mpi(basicN2 *mc, domain_t *domain){
    long i;
    int hops;
    mc->TXs[p] = 0;
    mc->TYs[p] = 0;
    mc->TZs[p] = 0;
    mc->TXs[n] = 0;
    mc->TYs[n] = 0;
    mc->TZs[n] = 0;
    mc->TXr[p] = 0;
    mc->TYr[p] = 0;
    mc->TZr[p] = 0;
    mc->TXr[n] = 0;
    mc->TYr[n] = 0;
    mc->TZr[n] = 0;

    /* Apply periodic boundary conditions */
    for (i = 0; i < mc->N; i++) {
        X(mc, t, i) = 0;
        Y(mc, t, i) = 0;
        Z(mc, t, i) = 0;

        if (X(mc, r, i) < 0) {
            X(mc, r, i) += domain->L[0];
            hops = find_hops(X(mc, r, i), domain->L[0]);
            X(mc, t, i) -= hops;
            mc->TXs[n]++;
        } else if (X(mc, r, i) >= domain->L[0]) {
            X(mc, r, i) -= domain->L[0];
            hops = find_hops(X(mc, r, i), domain->L[0]);
            X(mc, t, i) += hops;
            mc->TXs[p]++;
        }


        if (Y(mc, r, i) < 0) {
            Y(mc, r, i) += domain->L[1];
            hops = find_hops(Y(mc, r, i), domain->L[1]);
            Y(mc, t, i) -= hops;
            mc->TYs[n]++;
        } else if (Y(mc, r, i) >= domain->L[1]) {
            Y(mc, r, i) -= domain->L[1];
            hops = find_hops(Y(mc, r, i), domain->L[1]);
            Y(mc, t, i) += hops;
            mc->TYs[p]++;
        }


        if (Z(mc, r, i) < 0) {
            Z(mc, r, i) += domain->L[2];
            hops = find_hops(Z(mc, r, i), domain->L[2]);
            Z(mc, t, i) -= hops;
            mc->TZs[n]++;
        } else if (Z(mc, r, i) >= domain->L[2]) {
            Z(mc, r, i) -= domain->L[2];
            hops = find_hops(Z(mc, r, i), domain->L[2]);
            Z(mc, t, i) += hops;
            mc->TZs[p]++;
        }
    }
    long max = find_max(mc->TXs, mc->TYs, mc->TZs);
    if (max > mc->Max_transfer_sz) {
        mc->Max_transfer_sz = max;
        free(mc->send_buffer);
        mc->send_buffer = (mpi_buffer_t *) malloc(mc->Max_transfer_sz * sizeof(mpi_buffer_t));
    }
}

#else
static void apply_periodic_boundary(basicN2 *mc, domain_t *domain) {
    long i;

    /* Apply periodic boundary conditions */
    for (i = 0; i < mc->N; i++) {
        if (X(mc, r, i) < 0) {
            X(mc, r, i) += domain->L[0];
        } else if (X(mc, r, i) >= domain->L[0]) {
            X(mc, r, i) -= domain->L[0];
        }

        if (Y(mc, r, i) < 0) {
            Y(mc, r, i) += domain->L[1];
        } else if (Y(mc, r, i) >= domain->L[1]) {
            Y(mc, r, i) -= domain->L[1];
        }

        if (Z(mc, r, i) < 0) {
            Z(mc, r, i) += domain->L[2];
        } else if (Z(mc, r, i) >= domain->L[2]) {
            Z(mc, r, i) -= domain->L[2];
        }
    }
}

#endif

real basicN2_calc_E_kin(void *container) {
    basicN2 *mc = (basicN2 *) container;
    long i;
    real m  = config.m;

    real v2 = 0.;
    /* second sumand of v(t+dt) = v(t) + (F(t) + F(t+dt))/2m (new forces) */
    for(i = 0; i < mc->N; i++) {
        v2 += X(mc, v, i)*X(mc, v, i) + Y(mc, v, i)*Y(mc, v, i) + Z(mc, v, i)*Z(mc, v, i);
    }
    real E_kin = 0.5*m * v2;
    return E_kin;
}

void basicN2_update_halo(void *container, domain_t *domain) {
    basicN2 *mc = (basicN2 *) container;
    /* Delete old halo molecules */
    basicN2_clear(mc->halo);

    /* Determine the number of halo molecules. */
    long i;
    long N_halo = 0;
    for(i = 0; i < mc->N; i++) {
        int ximg = ( X(mc, r, i) < mc->cutoff_radius  ||  X(mc, r, i) >= domain->L[0] - mc->cutoff_radius );
        int yimg = ( Y(mc, r, i) < mc->cutoff_radius  ||  Y(mc, r, i) >= domain->L[1] - mc->cutoff_radius );
        int zimg = ( Z(mc, r, i) < mc->cutoff_radius  ||  Z(mc, r, i) >= domain->L[2] - mc->cutoff_radius );
        if(ximg) N_halo++;
        if(yimg) N_halo++;
        if(zimg) N_halo++;
        if(ximg && yimg) N_halo++;
        if(ximg && zimg) N_halo++;
        if(yimg && zimg) N_halo++;
        if(ximg && yimg && zimg) N_halo++;

    }
    LOG(Debug, "Number of halo molecules: %ld\n", N_halo);
    if( mc->halo->capacity < N_halo ) {
        DEBUG("Reallocating halo container with capacity %lu (old capacity was %lu)\n", N_halo, mc->halo->capacity);
        if (NSSD == 1) {
#ifdef GPU
            Dfree(dhr);
#else
#ifdef _OPENACC
            if (not_first){
               acc_delete(hsr, 3 * prev_halo_count * sizeof( real ));
            }
            prev_halo_count = N_halo;
            not_first = 1;
#endif
#endif
        }
        basicN2_free(mc->halo);
        mc->halo = NULL;

        mc->halo = basicN2_alloc_N(N_halo);
#ifdef MPI
        mc->halo->Nx[p] = mc->Nx[p];
        mc->halo->Ny[p] = mc->Ny[p];
        mc->halo->Nz[p] = mc->Nz[p];
        mc->halo->Nx[n] = mc->Nx[n];
        mc->halo->Ny[n] = mc->Ny[n];
        mc->halo->Nz[n] = mc->Nz[n];
        mc->halo->myrank = mc->myrank;
#endif
        if (NSSD == 1) {
            hsr = &X(mc->halo, r, 0);
            hcount[0] = N_halo;
            halo_count = N_halo;
#ifdef GPU
            Dmalloc((void **)&dhr, 3 * halo_count * sizeof( real ));
#else
#ifdef _OPENACC
            acc_create(hsr,3 * halo_count * sizeof( real ));
#endif
#endif
        }
    }
#ifdef MPI
    mc->halo->TXs[p] = 0;
    mc->halo->TYs[p] = 0;
    mc->halo->TZs[p] = 0;
    mc->halo->TXs[n] = 0;
    mc->halo->TYs[n] = 0;
    mc->halo->TZs[n] = 0;
    mc->halo->TXr[p] = 0;
    mc->halo->TYr[p] = 0;
    mc->halo->TZr[p] = 0;
    mc->halo->TXr[n] = 0;
    mc->halo->TYr[n] = 0;
    mc->halo->TZr[n] = 0;
#endif
    /* Copy molecules to halo container. */
    LOG(Debug, "Domain: %"PRIreal" x %"PRIreal" x %"PRIreal"\n", domain->L[0], domain->L[1], domain->L[2]);
    LOG(Debug, "Cutoff %"PRIreal"\n", mc->cutoff_radius);
    for(i = 0; i < mc->N; i++) {
        assert(X(mc, r, i) >= 0);
        assert(Y(mc, r, i) >= 0);
        assert(Z(mc, r, i) >= 0);
        assert(X(mc, r, i) < domain->L[0]);
        assert(Y(mc, r, i) < domain->L[1]);
        assert(Z(mc, r, i) < domain->L[2]);
#ifdef MPI
        X(mc->halo, t, mc->halo->N) = 0;
        Y(mc->halo, t, mc->halo->N) = 0;
        Z(mc->halo, t, mc->halo->N) = 0;
#endif
        int ximg = ( X(mc, r, i) < mc->cutoff_radius  ||  X(mc, r, i) >= domain->L[0] - mc->cutoff_radius );
        int yimg = ( Y(mc, r, i) < mc->cutoff_radius  ||  Y(mc, r, i) >= domain->L[1] - mc->cutoff_radius );
        int zimg = ( Z(mc, r, i) < mc->cutoff_radius  ||  Z(mc, r, i) >= domain->L[2] - mc->cutoff_radius );

      if (ximg) {
            copy_molecule(mc, i, mc->halo, mc->halo->N);
            if (X(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) += domain->L[0];
#ifdef MPI
                mc->halo->TXs[n]++;
                X(mc->halo, t, mc->halo->N)--;
#endif
            } else if (X(mc->halo, r, mc->halo->N) >= domain->L[0] - mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) -= domain->L[0];
#ifdef MPI
                mc->halo->TXs[p]++;
                X(mc->halo, t, mc->halo->N)++;
#endif
            }
            mc->halo->N++;
        }
        if (yimg) {
            copy_molecule(mc, i, mc->halo, mc->halo->N);
            if (Y(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) += domain->L[1];
#ifdef MPI
                mc->halo->TYs[n]++;
                Y(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Y(mc->halo, r, mc->halo->N) >= domain->L[1] - mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) -= domain->L[1];
#ifdef MPI
                mc->halo->TYs[p]++;
                Y(mc->halo, t, mc->halo->N)++;
#endif
            }
            mc->halo->N++;
        }
        if (zimg) {
            copy_molecule(mc, i, mc->halo, mc->halo->N);
            if (Z(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) += domain->L[2];
#ifdef MPI
                mc->halo->TZs[n]++;
                Z(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Z(mc->halo, r, mc->halo->N) >= domain->L[2] - mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) -= domain->L[2];
#ifdef MPI
                mc->halo->TZs[p]++;
                Z(mc->halo, t, mc->halo->N)++;
#endif
            }
            mc->halo->N++;
        }
        if (ximg && yimg) {
            copy_molecule(mc, i, mc->halo, mc->halo->N);
            if (X(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) += domain->L[0];
#ifdef MPI
                mc->halo->TXs[n]++;
                X(mc->halo, t, mc->halo->N)--;
#endif
            } else if (X(mc->halo, r, mc->halo->N) >= domain->L[0] - mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) -= domain->L[0];
#ifdef MPI
                mc->halo->TXs[p]++;
                X(mc->halo, t, mc->halo->N)++;
#endif
            }
            if (Y(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) += domain->L[1];
#ifdef MPI
                mc->halo->TYs[n]++;
                Y(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Y(mc->halo, r, mc->halo->N) >= domain->L[1] - mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) -= domain->L[1];
#ifdef MPI
                mc->halo->TYs[p]++;
                Y(mc->halo, t, mc->halo->N)++;
#endif
            }
            mc->halo->N++;
        }
        if (ximg && zimg) {
            copy_molecule(mc, i, mc->halo, mc->halo->N);
            if (X(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) += domain->L[0];
#ifdef MPI
                mc->halo->TXs[n]++;
                X(mc->halo, t, mc->halo->N)--;
#endif
            } else if (X(mc->halo, r, mc->halo->N) >= domain->L[0] - mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) -= domain->L[0];
#ifdef MPI
                mc->halo->TXs[p]++;
                X(mc->halo, t, mc->halo->N)++;
#endif
            }
            if (Z(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) += domain->L[2];
#ifdef MPI
                mc->halo->TZs[n]++;
                Z(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Z(mc->halo, r, mc->halo->N) >= domain->L[2] - mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) -= domain->L[2];
#ifdef MPI
                mc->halo->TZs[p]++;
                Z(mc->halo, t, mc->halo->N)++;
#endif
            }
            mc->halo->N++;
        }
        if (yimg && zimg) {
            copy_molecule(mc, i, mc->halo, mc->halo->N);
            if (Y(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) += domain->L[1];
#ifdef MPI
                mc->halo->TYs[n]++;
                Y(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Y(mc->halo, r, mc->halo->N) >= domain->L[1] - mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) -= domain->L[1];
#ifdef MPI
                mc->halo->TYs[p]++;
                Y(mc->halo, t, mc->halo->N)++;
#endif
            }
            if (Z(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) += domain->L[2];
#ifdef MPI
                mc->halo->TZs[n]++;
                Z(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Z(mc->halo, r, mc->halo->N) >= domain->L[2] - mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) -= domain->L[2];
#ifdef MPI
                mc->halo->TZs[p]++;
                Z(mc->halo, t, mc->halo->N)++;
#endif
            }
            mc->halo->N++;
        }
        if (ximg && yimg && zimg) {
            copy_molecule(mc, i, mc->halo, mc->halo->N);
            if (X(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) += domain->L[0];
#ifdef MPI
                mc->halo->TXs[n]++;
                X(mc->halo, t, mc->halo->N)--;
#endif
            } else if (X(mc->halo, r, mc->halo->N) >= domain->L[0] - mc->cutoff_radius) {
                X(mc->halo, r, mc->halo->N) -= domain->L[0];
#ifdef MPI
                mc->halo->TXs[p]++;
                X(mc->halo, t, mc->halo->N)++;
#endif
            }
            if (Y(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) += domain->L[1];
#ifdef MPI
                mc->halo->TYs[n]++;
                Y(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Y(mc->halo, r, mc->halo->N) >= domain->L[1] - mc->cutoff_radius) {
                Y(mc->halo, r, mc->halo->N) -= domain->L[1];
#ifdef MPI
                mc->halo->TYs[p]++;
                Y(mc->halo, t, mc->halo->N)++;
#endif
            }
            if (Z(mc->halo, r, mc->halo->N) < mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) += domain->L[2];
#ifdef MPI
                mc->halo->TZs[n]++;
                Z(mc->halo, t, mc->halo->N)--;
#endif
            } else if (Z(mc->halo, r, mc->halo->N) >= domain->L[2] - mc->cutoff_radius) {
                Z(mc->halo, r, mc->halo->N) -= domain->L[2];
#ifdef MPI
                mc->halo->TZs[p]++;
                Z(mc->halo, t, mc->halo->N)++;
#endif
            }
            mc->halo->N++;
        }
    }
#ifdef MPI
    long max = find_max(mc->halo->TXs, mc->halo->TYs, mc->halo->TZs);

    if (max > mc->halo->Max_transfer_sz) {
        mc->halo->Max_transfer_sz = max;
        free(mc->halo->send_buffer);
        mc->halo->send_buffer = (mpi_buffer_t*)malloc(mc->halo->Max_transfer_sz * sizeof(mpi_buffer_t));
    }
    if (Sx > 1) {
        basicN2_send_buffer_update(mc->halo, mc->myrank, X, p);
        basicN2_recv_buffer_update(mc->halo, mc->halo->TXr[n]);
        basicN2_buffer_exchange(mc->halo, mc->myrank, X, p);

        basicN2_send_buffer_update(mc->halo, mc->myrank, X, n);
        basicN2_recv_buffer_update(mc->halo, mc->halo->TXr[p]);
        basicN2_buffer_exchange(mc->halo, mc->myrank, X, n);
    }
    if (Sy > 1) {
        basicN2_send_buffer_update(mc->halo, mc->myrank, Y, p);
        basicN2_recv_buffer_update(mc->halo, mc->halo->TYr[n]);
        basicN2_buffer_exchange(mc->halo, mc->myrank, Y, p);

        basicN2_send_buffer_update(mc->halo, mc->myrank, Y, n);
        basicN2_recv_buffer_update(mc->halo, mc->halo->TYr[p]);
        basicN2_buffer_exchange(mc->halo, mc->myrank, Y, n);
    }
    if (Sz > 1) {
        basicN2_send_buffer_update(mc->halo, mc->myrank, Z, p);
        basicN2_recv_buffer_update(mc->halo, mc->halo->TZr[n]);
        basicN2_buffer_exchange(mc->halo, mc->myrank, Z, p);

        basicN2_send_buffer_update(mc->halo, mc->myrank, Z, n);
        basicN2_recv_buffer_update(mc->halo, mc->halo->TZr[p]);
        basicN2_buffer_exchange(mc->halo, mc->myrank, Z, n);
    }
#endif
}

void basicN2_print_ascii(void *container, char *filename) {
    basicN2 *mc = (basicN2 *) container;
    FILE *fh;
    fh = fopen(filename, "w+");
    
    long i;
    for(i = 0; i < mc->N; i++) {
        fprintf(fh, "%lu\t%"PRIreal", %"PRIreal", %"PRIreal"\t%"PRIreal", %"PRIreal", %"PRIreal"\t%"PRIreal", %"PRIreal", %"PRIreal"\n",
                ID(mc,i), X(mc,r,i), Y(mc,r,i), Z(mc,r,i), X(mc,v,i), Y(mc,v,i), Z(mc,v,i), X(mc,F,i), Y(mc,F,i), Z(mc,F,i) );
    }
    fclose(fh);
}

void basicN2_apply_thermostat(void *container, real T_cur, real T_target) {
    assert(container != NULL);
    assert( T_cur > 0.0 );
    assert( T_target > 0.0 );
    
    long i;
    basicN2 *mc = (basicN2 *) container;    
    real f = sqrt(T_target / T_cur);
    LOG(Info, "Thermorstat scaling factor: f = %"PRIreal"\n", f);
    if(f > VELOCITY_SCALING_FACTOR_MAX || f < VELOCITY_SCALING_FACTOR_MIN) {
        LOG(Warning, "Thermostat scaling factor outside of range [%"PRIreal",%"PRIreal"]: %"PRIreal"\n",
            VELOCITY_SCALING_FACTOR_MIN, VELOCITY_SCALING_FACTOR_MAX, f);
        real E_kin_limit = VELOCITY_SCALING_E_KIN_LIMIT_FACTOR * T_target;
        for(i = 0; i < mc->N; i++) {
            real E_kin = 0.5*config.m* (X(mc, v, i)*X(mc, v, i) + Y(mc, v, i)*Y(mc, v, i) + Z(mc, v, i)*Z(mc, v, i));
            if(E_kin > E_kin_limit) {
                LOG(Warning, "E_kin: %"PRIreal"/%"PRIreal"\n", E_kin, E_kin_limit);
                real ff = sqrt(E_kin_limit/E_kin);
                X(mc, v, i) *= ff;
                Y(mc, v, i) *= ff;
                Z(mc, v, i) *= ff;
            }
        }
    }
    else {
        for(i = 0; i < mc->N; i++) {
            X(mc, v, i) *= f;
            Y(mc, v, i) *= f;
            Z(mc, v, i) *= f;
        }
    }
}

long basicN2_get_num_molecules(void *container) {
    basicN2 *mc = (basicN2 *) container;
    return mc->N;
}

void basicN2_update(void *container, domain_t *domain){
    basicN2 * mc = (basicN2 *) container;
#ifdef MPI
    apply_periodic_boundary_mpi(mc, domain);
    if (Sx > 1) {
        basicN2_send_buffer_update(mc, mc->myrank, X, p);
        basicN2_recv_buffer_update(mc, mc->TXr[n]);
        basicN2_buffer_exchange(mc, mc->myrank, X, p);

        basicN2_send_buffer_update(mc, mc->myrank, X, n);
        basicN2_recv_buffer_update(mc, mc->TXr[p]);
        basicN2_buffer_exchange(mc, mc->myrank, X, n);
    }
    if (Sy > 1) {
        basicN2_send_buffer_update(mc, mc->myrank, Y, p);
        basicN2_recv_buffer_update(mc, mc->TYr[n]);
        basicN2_buffer_exchange(mc, mc->myrank, Y, p);

        basicN2_send_buffer_update(mc, mc->myrank, Y, n);
        basicN2_recv_buffer_update(mc, mc->TYr[p]);
        basicN2_buffer_exchange(mc, mc->myrank, Y, n);
    }
    if (Sz > 1) {
        basicN2_send_buffer_update(mc, mc->myrank, Z, p);
        basicN2_recv_buffer_update(mc, mc->TZr[n]);
        basicN2_buffer_exchange(mc, mc->myrank, Z, p);

        basicN2_send_buffer_update(mc, mc->myrank, Z, n);
        basicN2_recv_buffer_update(mc, mc->TZr[p]);
        basicN2_buffer_exchange(mc, mc->myrank, Z, n);
    }
#else
    apply_periodic_boundary(mc, domain);
#endif
}

void basicN2_reset_forces_and_momenta(void *container) {
   basicN2 *mc = (basicN2 *) container;
    long i,j;
    long N = (mc->capacity / NSSD) + nfactor;

#ifdef GPU
    for (i = 0; i < NSSD - CPU_PART; i++) {
        Upot_buffer[i]=0;
        HDcopyAsyncStreams(&dU_pot[i], &Upot_buffer[i], sizeof(real),i);
        cuBasicN2_reset_forces(df, count, N, i);
    }
#pragma omp parallel if (CPU_PART>1)
#pragma omp for schedule(runtime)
    for(i = NSSD - CPU_PART; i < NSSD; i++) {
       long size = count[i];
#pragma omp task 
{
       for (j = 0; j < size; j++) {
           sf[i * 3 * N + j        ] = 0;
           sf[i * 3 * N + j + N    ] = 0;
           sf[i * 3 * N + j + N * 2] = 0;
       }
}
    }

#else
#ifndef _OPENACC
#pragma omp parallel if (NSSD>1)
#pragma omp for schedule(runtime)
#endif
    for(i = 0; i < NSSD; i++) {
       long size = count[i];
#ifdef _OPENACC
#pragma acc parallel loop present(sf[i * 3 * N:3 * N]) vector_length(256) async(i)
#else
#pragma omp task 
{
#endif
       for (j = 0; j < size; j++) {
           sf[i * 3 * N + j        ] = 0;
           sf[i * 3 * N + j + N    ] = 0;
           sf[i * 3 * N + j + N * 2] = 0;
       }
#ifndef _OPENACC
}
#endif
    }
#endif
}

void basicN2_shift_velocity(void *container, real v[3]) {
    basicN2 *mc = (basicN2 *) container;
    long i;
    for(i = 0; i < mc->N; i++) {
        X(mc, v, i) -= v[0];
        Y(mc, v, i) -= v[1];
        Z(mc, v, i) -= v[2];
    }
}

void basicN2_print_stats(void *container, ensemble_t *ensemble) { 
    basicN2 *mc = (basicN2 *) container;
    LOG(Info, "Capacity: %ld, Size: %ld\n", mc->capacity, mc->N);
}

void basicN2_print_pov(void *container, char *filename) {
    basicN2 *mc = (basicN2 *) container;
    FILE *fh;
    fh = fopen(filename, "w+");
    
    fprintf(fh, "#include \"psp-header.inc\"\n");
    long i;
    for(i = 0; i < mc->N; i++) {
        fprintf(fh, "sphere { <%"PRIreal", %"PRIreal", %"PRIreal">, 0.15   pigment { color rgb <1, 0, 0> } }\n", 
                X(mc,r,i), Y(mc,r,i), Z(mc,r,i) );
    }
    fclose(fh);
}

void basicN2_clear_halo(void *container) {
    basicN2 *mc = (basicN2 *) container;
    mc->halo->N = 0;
}

inline static void halo_counter(long x, long y, long z, int * hcounter){
#pragma omp atomic update
    hcounter[z * SSx * SSy + y * SSx + x ]++;	
}

inline static void halo_add(basicN2 * mc, long x, long y, long z, long i, int offset){
    long old =0;
#pragma omp atomic capture
    old = (hcount[z * SSx * SSy + y * SSx + x]++);
    assert(hcount[z * SSx * SSy + y * SSx + x]<halo_count+hfactor);
    hsr[old + ( z * SSx * SSy + y * SSx + x ) * 3 * offset             ] = X(mc, r, i);
    hsr[old + ( z * SSx * SSy + y * SSx + x ) * 3 * offset + offset    ] = Y(mc, r, i);
    hsr[old + ( z * SSx * SSy + y * SSx + x ) * 3 * offset + offset * 2] = Z(mc, r, i);
}

inline static void halo_check(basicN2 *mc, int * hcounter, long m, long k, long j, real DomLengthX, real DomLengthY, real DomLengthZ, real cutoff_radius, int hN, int i){
    long _m = m;
    long m_ = m;
    long _k = k;
    long k_ = k;
    long _j = j;
    long j_ = j;

    int ximg_ = 0;
    int yimg_ = 0;
    int zimg_ = 0;
    int _ximg = 0;
    int _yimg = 0;
    int _zimg = 0;
    if (m < SSx - 1){ ximg_ = ((X(mc,r,i) - (m + 1) * DomLengthX + cutoff_radius) >=  0 && (X(mc,r,i) - (m + 1) * DomLengthX ) < 0);}
    if (k < SSy - 1){ yimg_ = ((Y(mc,r,i) - (k + 1) * DomLengthY + cutoff_radius) >=  0 && (Y(mc,r,i) - (k + 1) * DomLengthY ) < 0);}
    if (j < SSz - 1){ zimg_ = ((Z(mc,r,i) - (j + 1) * DomLengthZ + cutoff_radius) >=  0 && (Z(mc,r,i) - (j + 1) * DomLengthZ ) < 0);}
    if (ximg_)                  { m_++; if (hcounter == NULL) halo_add(mc, m_, k , j , i, hN); else halo_counter(m_, k , j , hcounter);}
    if (yimg_)                  { k_++; if (hcounter == NULL) halo_add(mc, m , k_, j , i, hN); else halo_counter(m , k_, j , hcounter);}
    if (zimg_)                  { j_++; if (hcounter == NULL) halo_add(mc, m , k , j_, i, hN); else halo_counter(m , k , j_, hcounter);}
    if (ximg_ && yimg_)         {       if (hcounter == NULL) halo_add(mc, m_, k_, j , i, hN); else halo_counter(m_, k_, j , hcounter);}
    if (ximg_ && zimg_)         {       if (hcounter == NULL) halo_add(mc, m_, k , j_, i, hN); else halo_counter(m_, k , j_, hcounter);}
    if (yimg_ && zimg_)         {       if (hcounter == NULL) halo_add(mc, m , k_, j_, i, hN); else halo_counter(m , k_, j_, hcounter);}
    if (ximg_ && yimg_ && zimg_){       if (hcounter == NULL) halo_add(mc, m_, k_, j_, i, hN); else halo_counter(m_, k_, j_, hcounter);}
    if (m > 0){ _ximg = ((X(mc,r,i) - m * DomLengthX) >= 0 && (X(mc,r,i) - m * DomLengthX) < cutoff_radius);}
    if (k > 0){ _yimg = ((Y(mc,r,i) - k * DomLengthY) >= 0 && (Y(mc,r,i) - k * DomLengthY) < cutoff_radius);}
    if (j > 0){ _zimg = ((Z(mc,r,i) - j * DomLengthZ) >= 0 && (Z(mc,r,i) - j * DomLengthZ) < cutoff_radius);}
    if (_ximg)                  { _m--; if (hcounter == NULL) halo_add(mc, _m, k , j , i, hN); else halo_counter(_m, k , j , hcounter);}
    if (_yimg)                  { _k--; if (hcounter == NULL) halo_add(mc, m , _k, j , i, hN); else halo_counter(m , _k, j , hcounter);}
    if (_zimg)                  { _j--; if (hcounter == NULL) halo_add(mc, m , k , _j, i, hN); else halo_counter(m , k , _j, hcounter);}
    if (_ximg && _yimg)         {       if (hcounter == NULL) halo_add(mc, _m, _k, j , i, hN); else halo_counter(_m, _k, j , hcounter);}
    if (_ximg && _zimg)         {       if (hcounter == NULL) halo_add(mc, _m, k , _j, i, hN); else halo_counter(_m, k , _j, hcounter);}
    if (_yimg && _zimg)         {       if (hcounter == NULL) halo_add(mc, m , _k, _j, i, hN); else halo_counter(m , _k, _j, hcounter);}
    if (_ximg && _yimg && _zimg){       if (hcounter == NULL) halo_add(mc, _m, _k, _j, i, hN); else halo_counter(_m, _k, _j, hcounter);}

    if (ximg_ && _yimg)         {       if (hcounter == NULL) halo_add(mc, m_, _k, j , i, hN); else halo_counter(m_, _k, j , hcounter);}
    if (ximg_ && _zimg)         {       if (hcounter == NULL) halo_add(mc, m_, k , _j, i, hN); else halo_counter(m_, k , _j, hcounter);}
    if (_ximg && yimg_)         {       if (hcounter == NULL) halo_add(mc, _m, k_, j , i, hN); else halo_counter(_m, k_, j , hcounter);}
    if (_ximg && zimg_)         {       if (hcounter == NULL) halo_add(mc, _m, k , j_, i, hN); else halo_counter(_m, k , j_, hcounter);}
    if (yimg_ && _zimg)         {       if (hcounter == NULL) halo_add(mc, m , k_, _j, i, hN); else halo_counter(m , k_, _j, hcounter);}
    if (_yimg && zimg_)         {       if (hcounter == NULL) halo_add(mc, m , _k, j_, i, hN); else halo_counter(m , _k, j_, hcounter);}
    if (ximg_ && _yimg && _zimg){       if (hcounter == NULL) halo_add(mc, m_, _k, _j, i, hN); else halo_counter(m_, _k, _j, hcounter);}
    if (_ximg && yimg_ && zimg_){       if (hcounter == NULL) halo_add(mc, _m, k_, j_, i, hN); else halo_counter(_m, k_, j_, hcounter);}
    if (_ximg && yimg_ && _zimg){       if (hcounter == NULL) halo_add(mc, _m, k_, _j, i, hN); else halo_counter(_m, k_, _j, hcounter);}
    if (ximg_ && _yimg && zimg_){       if (hcounter == NULL) halo_add(mc, m_, _k, j_, i, hN); else halo_counter(m_, _k, j_, hcounter);}
    if (ximg_ && yimg_ && _zimg){       if (hcounter == NULL) halo_add(mc, m_, k_, _j, i, hN); else halo_counter(m_, k_, _j, hcounter);}
    if (_ximg && _yimg && zimg_){       if (hcounter == NULL) halo_add(mc, _m, _k, j_, i, hN); else halo_counter(_m, _k, j_, hcounter);}
}

inline static void domain_scan(basicN2 * mc, domain_t * domain, long i){
    long old, k, j, m;
    real DomLengthX = (domain->L[0]  / (real)SSx);
    real DomLengthY = (domain->L[1]  / (real)SSy);
    real DomLengthZ = (domain->L[2]  / (real)SSz);
    m = (long) (X(mc,r,i)/DomLengthX);
    k = (long) (Y(mc,r,i)/DomLengthY);
    j = (long) (Z(mc,r,i)/DomLengthZ);
    long N = ( mc->capacity / NSSD ) + nfactor;
#pragma omp atomic capture
    old = (count[j * SSx * SSy + k * SSx + m]++);
    assert(count[j * SSx * SSy + k * SSx + m]< (mc->capacity/NSSD)+nfactor);
    sf  [old + ( j * SSx * SSy + k * SSx + m ) * N * 3        ] = X(mc, F, i);
    sf  [old + ( j * SSx * SSy + k * SSx + m ) * N * 3 + N    ] = Y(mc, F, i);
    sf  [old + ( j * SSx * SSy + k * SSx + m ) * N * 3 + N * 2] = Z(mc, F, i);
    sr  [old + ( j * SSx * SSy + k * SSx + m ) * N * 3        ] = X(mc, r, i);
    sr  [old + ( j * SSx * SSy + k * SSx + m ) * N * 3 + N    ] = Y(mc, r, i);
    sr  [old + ( j * SSx * SSy + k * SSx + m ) * N * 3 + N * 2] = Z(mc, r, i);
    indx[old + ( j * SSx * SSy + k * SSx + m ) * N            ] = i;
}

inline static void internal_halo_scan(basicN2 * mc, domain_t * domain, real cutoff_radius, long i, int * hcounter){
    long k, j, m;
    real DomLengthX = (domain->L[0]  / (real)SSx);
    real DomLengthY = (domain->L[1]  / (real)SSy);
    real DomLengthZ = (domain->L[2]  / (real)SSz);
    m = (long) (X(mc,r,i)/DomLengthX);
    k = (long) (Y(mc,r,i)/DomLengthY);
    j = (long) (Z(mc,r,i)/DomLengthZ);
    long hN = halo_count + hfactor;
    halo_check(mc, hcounter, m, k, j, DomLengthX, DomLengthY, DomLengthZ, cutoff_radius, hN, i);
}

inline static void external_halo_scan(basicN2 * mc, domain_t * domain, real cutoff_radius, long i, int * hcounter){
    long k, j, m;
    real DomLengthX = (domain->L[0]  / (real)SSx);
    real DomLengthY = (domain->L[1]  / (real)SSy);
    real DomLengthZ = (domain->L[2]  / (real)SSz);
    m = (long) (X(mc,r,i)/DomLengthX);
    k = (long) (Y(mc,r,i)/DomLengthY);
    j = (long) (Z(mc,r,i)/DomLengthZ);
    long hN = halo_count + hfactor;
    if (m == SSx) m--;
    if (k == SSy) k--;
    if (j == SSz) j--;
    if (hcounter == NULL) halo_add(mc, m, k, j, i, hN); else  halo_counter(m, k, j, hcounter);
    halo_check(mc, hcounter, m, k, j, DomLengthX, DomLengthY, DomLengthZ, cutoff_radius, hN, i);
}

void basicN2_estimate_halos( void * container, domain_t * domain){
    basicN2 *mc = (basicN2 *) container;
    long  i;
    if(NSSD > 1){
      int *counter=(int*) malloc(NSSD*sizeof(int));
      for (i = 0; i < NSSD; i++){
  	counter[i]=0;
      }
#pragma omp parallel default(shared) private (i) 
{ 
#pragma omp for schedule(runtime) 
      for (i = 0; i < mc->N; i++) {
          internal_halo_scan( mc, domain, config.cutoff_radius, i, counter);
      }
#pragma omp for schedule(runtime)
      for (i = 0; i < mc->halo->N; i++) {
          external_halo_scan(mc->halo, domain, config.cutoff_radius, i, counter);
      }
}
      for( i = 0; i < NSSD; i++){
  	if (halo_count < counter[i] ){
  	   halo_count= counter[i];
  	}
      }
      free(counter);
#ifdef GPU
      Dmalloc ((void **)&dhr, 3 * (NSSD - CPU_PART) * (halo_count + hfactor) * sizeof( real ));
      HmallocP((void **)&hsr, 3 * NSSD * (halo_count + hfactor) * sizeof( real ));
#else
      hsr = (real *)   malloc(3 * NSSD * (halo_count + hfactor) * sizeof( real ));
#ifdef _OPENACC
      acc_create(hsr, 3 * NSSD * (halo_count + hfactor) * sizeof( real ));
#endif
#endif

   }
#ifdef GPU
   HDcopyAsync(dr, sr, 3 * (NSSD - CPU_PART) * ((mc->capacity / NSSD ) + nfactor) * sizeof( real ));
#else
#ifdef _OPENACC
   Num =  3 * NSSD * ((mc->capacity / NSSD ) + nfactor);
#pragma acc update device(sr[0:Num]) async
#endif
#endif
}

void basicN2_decompose(void * container, domain_t *domain){
    basicN2 * mc = (basicN2 *)container;
    if (NSSD > 1) {
        long i;
        for (i = 0; i < NSSD; i++) {
            count[i] = 0;
            hcount[i] = 0;
        }
#pragma omp parallel default(shared) private (i) 
{ 
#pragma omp for schedule(runtime) 
        for (i = 0; i < mc->N; i++) {
            domain_scan(mc, domain, i);
            internal_halo_scan(mc, domain, config.cutoff_radius, i, NULL);
        }
#pragma omp for schedule(runtime)
        for (i = 0; i < mc->halo->N; i++) {
            external_halo_scan(mc->halo, domain, config.cutoff_radius, i, NULL);
        }
}
    } else {
       count[0] = mc->N;
       hcount[0] = mc->halo->N;
    }
#ifdef GPU
    long hN = halo_count + hfactor;
    long  N = (mc->capacity / NSSD) + nfactor;
    HDcopyAsync(dhr, hsr, 3 * (NSSD - CPU_PART) * hN * sizeof( real )); 
    HDcopyAsync(dr,  sr, (NSSD - CPU_PART) * 3 * N  * sizeof( real )); 
#else
#ifdef _OPENACC
    long hNum = 3 * NSSD * (halo_count + hfactor);
    long  Num = 3 * NSSD * ((mc->capacity / NSSD) + nfactor);
#pragma acc update device(hsr[0:hNum]) async
#pragma acc update device(sr[0 : Num]) async
#endif
#endif
}

void basicN2_destroy(void *container){
#ifdef GPU
#ifdef USE_STREAMS
    cuBasicN2_destroy_streams(NSSD);
#endif
    Dfree(dr);
    Dfree(dhr);
    Dfree(df);
    Dfree(dcount);
    Dfree(dhcount);
    Dfree(dU_pot);
    Hfree(Upot_buffer);
    if(NSSD > 1){
      Hfree(sr);
      Hfree(sf);
      Hfree(hsr);
    }
#else
    if (NSSD > 1){
       free(indx);
       free(sr);
       free(hsr);
       free(sf);
    }
    free(Upot_buffer);
#ifdef _OPENACC
    basicN2 * mc = (basicN2 *)container;
    acc_delete(hsr,3 * NSSD * (halo_count + hfactor) * sizeof( real ));
    acc_delete(sr, 3 * NSSD * (( mc->N / NSSD ) + nfactor) * sizeof( real ));
    acc_delete(sf, 3 * NSSD * (( mc->N / NSSD ) + nfactor) * sizeof( real ));
    acc_delete(Upot_buffer,  NSSD * sizeof( real ));
#endif
#endif
    free(count);
    free(hcount);
#ifdef MPI
   MPI_Type_free(&MPI_MOL_T);
#endif
}


