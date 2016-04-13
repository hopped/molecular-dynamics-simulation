/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include "config.h"
#include "typedefs.h"
typedef enum {
	VELOCITY_SCALING = 0,
	THERMOSTAT_MAX
} thermostat_type_t;

typedef struct {
	real start_time;
	real end_time;
	real dt;
    real equilibration_time;
    thermostat_type_t thermostat_type;
} simulation_t;

#ifndef SIMULATION_SRC
extern simulation_t simulation;
#endif

void simulation_alloc(simulation_t *simulation);
void simulation_free(simulation_t *simulation);

void simulation_init(simulation_t *simulation);

void simulation_read_xml(simulation_t *simulation /* , xmlconfig */);

void simulation_set_thermostat(simulation_t *simulation, char *thermostat);

#endif
