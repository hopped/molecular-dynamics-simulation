/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

/**
 * @file
 *
 * Timer.
 */

#ifndef TIMER_H
#define TIMER_H

#include "sys/time.h"

static inline double timer(){
  struct timeval timer;
  gettimeofday(&timer, NULL);
  return (double) timer.tv_sec + (double) timer.tv_usec / 1.0e6;
}

#endif
