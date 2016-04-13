/*/*
 * Copyright (c) 2014      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef DROPOVERBASIN_GENERATOR_H
#define DROPOVERBASIN_GENERATOR_H

#include "domain.h"
#include "ensemble.h"
#include "molecule.h"
#include "phasespace.h"
#include "utils/generator/Lattice.h"

void dropOverBasin_generator(ensemble_t* ensemble, double hBasin, double dropPos[3], double rDrop, domain_t* domain, phasespace_t* psp);

#endif /* DROPOVERBASIN_GENERATOR_H */
