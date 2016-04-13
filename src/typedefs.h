/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <stdio.h>
//#define USE_FLOAT
#ifdef USE_FLOAT
typedef float real;
#ifdef MPI
#define  MPIREAL MPI_FLOAT
#endif
#define PRIreal "f"
#else
typedef double real;
#ifdef MPI
#define  MPIREAL MPI_DOUBLE
#endif
#define PRIreal "lf"
#endif


#endif /* TYPEDEFS_H */
