/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef GRID_GENERATOR_H
#define GRID_GENERATOR_H

#include "domain.h"
#include "ensemble.h"
#include "phasespace.h"
#include "utils/generator/Lattice.h"

void grid_generator(ensemble_t *ensemble, LatticeSystem system, LatticeCentering centering, domain_t *domain, phasespace_t *psp);

#endif /* GRID_GENERATOR_H */
