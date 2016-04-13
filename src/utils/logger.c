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
 * Logger.
 */

#define TIMER_SRC

#include "logger.h"
#include "timer.h"

#include <stdio.h>  
#include <stdarg.h>

static char *logLevelNames[] = {
    [None]    = "None",
    [Fatal]   = "Fatal",
    [Error]   = "Error",
    [Warning] = "Warning",
    [Info]    = "Info",
    [Debug]   = "Debug",
    [All]     = "All"
};

/* The global logger used in logger macros. */
Logger global_log;

void logger_init(Logger *logger, char *filename) {
    logger->start_time = timer();
    logger->fh = fopen(filename, "w+");
    logger->output_level = Warning;
}

void logger_finalize(Logger *logger) {
    fclose(logger->fh);
}

void logger_msg(Logger *logger, LogLevel level, char *filename, int line, char *format, ...) {
    if(logger->output_level < level) { /* skip if log level is lower then message level */
        return ;
    }
    double time = timer();
    time -= logger->start_time;

    va_list va_args;
    va_start(va_args, format);

    char msg[MAX_LOG_MSG_SIZE];
    vsnprintf(msg, MAX_LOG_MSG_SIZE, format, va_args);
#ifndef NDEBUG
    fprintf(logger->fh, "%lf [%s] %s[%d]:\t%s", time, logLevelNames[level], filename, line, msg);
#else
    fprintf(logger->fh, "%lf [%s]:\t%s", time, logLevelNames[level], msg);
#endif
    va_end(va_args);
}


void logger_set_output_level(Logger *logger, LogLevel level) {
    logger->output_level = level;
}

LogLevel logger_get_output_level(Logger *logger) {
    return logger->output_level;
}
