/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef DOMAIN_H
#define DOMAIN_H

#include "config.h"
#include "typedefs.h"

/** domain types/shapes */
typedef enum {
	UNKNOWN,
	CUBE,       /**< cubic domain */
	CUBOID      /**< cuboidal domain */
} domain_types_t;

/** domain 
 *
 * stores information about the simulation area
 */
typedef struct {
	domain_types_t type; /**< domain type */
    real L[3];           /**< dimemsions of the domain */
} domain_t;

#ifndef DOMAIN_SRC

extern domain_t domain; /**< global accessible domain object */
#endif

/** set domain type 
 * param[in]   domain_type  domain type (cube,rectangle)
 * param[out]  domain       domain whose type shall be set
 */
void domain_set_type(domain_t *domain, char *domain_type);

/** set domain size
 * param[in]   L            size {x,y,z}
 * param[out]  domain       domain whose size shall be set
 */
void domain_set_size(domain_t *domain, real L[3]);

real domain_get_volume(domain_t *domain);

#endif /* DOMAIN_H */
