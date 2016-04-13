/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef MOLECULE_H
#define MOLECULE_H

#include "config.h"
#include "typedefs.h"

/** molecule type
 *  
 *  The molecule_t can store the data for a single molecule.
 */
typedef struct {
    long cid; /**< component id */
    real r[3]; /**< coordinate of center of mass (x,y,z) */
    real F[3]; /**< force actnig on the molecule */
    real v[3]; /**< velocity of center of mass (x,y,z) */
    long id; /**< molecule id */
} molecule_t;

#endif /* MOLECULE_H */
