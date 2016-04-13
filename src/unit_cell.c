/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "unit_cell.h"
#include "unit_cell_site.h"

#include <stdlib.h>

void unit_cell_site_alloc(unit_cell_t *unit_cell) { }

void unit_cell_site_free(unit_cell_t *unit_cell) {
	free(unit_cell->sites);
}

void unit_cell_site_init(unit_cell_t *unit_cell){}
void unit_cell_site_print_ascii(unit_cell_t *unit_cell, FILE *fh){}

void unit_cell_site_add_site(unit_cell_t *unit_cell, unit_cell_site_t *site){}
