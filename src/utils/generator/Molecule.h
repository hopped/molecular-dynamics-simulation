/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef GENERATOR_MOLECULE_H
#define GENERATOR_MOLECULE_H

/** molecule type
 *
 *  The molecule_t can store the data for a single molecule.
 */
typedef struct generator_molecule_t {
    long cid; /**< component id */
    double r[3]; /**< coordinate of center of mass (x,y,z) */
} generator_molecule_t;


#endif /* MOLECULE_H */
