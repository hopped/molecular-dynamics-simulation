/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef UNIT_CELL_H
#define UNIT_CELL_H

#include "unit_cell_site.h"

#include <stdio.h>

typedef struct unit_cell_t {
    char name[256];
    int num_sites;      /**< number of unit cell sites */
    unit_cell_site_t sites[4];     /* unit cell sites */
} unit_cell_t;

void unit_cell_site_alloc(unit_cell_t *unit_cell);
void unit_cell_site_free(unit_cell_t *unit_cell);
void unit_cell_site_init(unit_cell_t *unit_cell);
void unit_cell_site_print_ascii(unit_cell_t *unit_cell, FILE *fh);

void unit_cell_site_add_site(unit_cell_t *unit_cell, unit_cell_site_t *site);


#endif /* UNIT_CELL_H */
