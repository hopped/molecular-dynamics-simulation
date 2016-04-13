/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#define DOMAIN_SRC
#include "domain.h"

#include "utils/logger.h"

domain_t domain;


#include <stdlib.h>
#include <string.h>

void domain_set_type(domain_t *domain, char *domain_type) {
	if( strcmp("cube", domain_type) == 0) {
		domain->type = CUBE;
	} 
	else if ( strcmp("cuboid", domain_type) == 0 ) {
		domain->type = CUBOID;
	}
	else {
		LOG(Error, "Unknown domain type: %s\n", domain_type);
		exit(1);
	}
}

void domain_set_size(domain_t *domain, real L[3]) {
    domain->L[0] = L[0];
    domain->L[1] = L[1];
    domain->L[2] = L[2];
}

real domain_get_volume(domain_t *domain) {
    real V = 0;

    real *L = domain->L;
    switch(domain->type) {
        case CUBE:
            V = L[0] * L[0] * L[0];
            break;
        case CUBOID:
            V = L[0] * L[1] * L[2];
            break;
        default:
            LOG(Error, "Unknown domain type: %u\n", domain->type);
            exit(1);
    }
    return V;
}
