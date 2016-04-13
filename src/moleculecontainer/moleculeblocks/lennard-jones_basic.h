/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef LENNARD_JONES_BASIC_H
#define LENNARD_JONES_BASIC_H

#include "typedefs.h"
#include "molecule_block.h"
#include "moleculeblocks.h"
#include "ensemble.h"

void lj_basic_block_calc(molecule_block_t *block1, molecule_block_t *blokc2, real *U_pot);

void lj_basic_diagblock_calc(molecule_block_t *block, real *U_pot);

#endif /* LENNARD_JONES_BASIC_H */
