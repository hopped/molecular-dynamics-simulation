/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#define SIMULATION_SRC

#include "simulation.h"
#include "utils/logger.h"

#include <stdlib.h>
#include <string.h>

simulation_t simulation;

void simulation_alloc(simulation_t *simulation) {
    simulation = (simulation_t *) malloc(sizeof(simulation_t));
}

void simulation_free(simulation_t *simulation) {
    free(simulation);
}

void simulation_init(simulation_t *simulation){}

void simulation_read_xml(simulation_t *simulation){}

void simulation_set_thermostat(simulation_t *simulation, char *thermostat) {
	if( strcmp( "velocity-scaling", thermostat ) == 0 ) {
		simulation->thermostat_type = VELOCITY_SCALING;
	} else {
		LOG(Error, "Unknown thermostat: %s\n", thermostat);
		exit(1);
	}
}
