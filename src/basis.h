/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef BASIS_H
#define BASIS_H

#include "config.h"
#include "typedefs.h"

/** Structure holding the basis used within a unit cell */
typedef struct {
    int N;      /**< number of atoms/molecules in the basis */
    int *cid;   /**< component ID of the atoms/molecules */
    real *r[3]; /**< relative coordinates of the atoms/molecules */
} basis_t;

void basis_init(int N, basis_t *basis);
void basis_set_site(int i, int compoentn_id, real r[3], basis_t *basis);
void basis_print_info(basis_t *basis);
void basis_free(basis_t *basis);

#endif /* BASIS_H */

