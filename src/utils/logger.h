/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>

#define GLOBAL_LOG (&global_log)

#define INIT_GLOBAL_LOG(FILENAME) do { \
    logger_init(GLOBAL_LOG, FILENAME); \
} while(0);
#define FINALIZE_GLOBAL_LOG() do { \
    logger_finalize(GLOBAL_LOG); \
} while(0);

#define LOG(level, ...) logger_msg(GLOBAL_LOG, level, __FILE__, __LINE__, __VA_ARGS__)

#ifndef NDEBUG
#define DEBUG(...)      logger_msg(GLOBAL_LOG, Debug, __FILE__, __LINE__, __VA_ARGS__)
#else
#define DEBUG(...) 
#endif

#define MAX_LOG_MSG_SIZE 4096

typedef enum {
    None    = 0,  /* supress output */
    Fatal   = 1,  /* program exit */
    Error   = 2,  /* program corrected */
    Warning = 3,  /* perhaps wrong */
    Info    = 4,  /* user info */
    Debug   = 5,  /* detailed info for debugging */
    All
} LogLevel; /**< Log message levels.*/

typedef struct Logger {
    FILE *fh;
    double start_time;
    LogLevel output_level;
} Logger; /**< Logger object */


extern Logger global_log; /**< global logger object used by logger macors. */

void logger_init(Logger *logger, char *filename); /**< Initialize a new logger. */
void logger_finalize(Logger *logger); /**< Finalize logger. */

/** Create New log message entry. 
 *
 * @param[in]  logger    Logger to use for output.
 * @param[in]  level     Log level at which this message should occure
 * @param[in]  filename  Filename corresponding to the message. (e.g. __FILE__)
 * @param[in]  line      Line number corresponding to the message. (e.g. __LINE__)
 * @param[in]  format    Format string as used in printf.
 * @param[in]  ...       Argumentlist as used in printf format.
 */
void logger_msg(Logger *logger, LogLevel level, char *filename, int line, char *format, ...);

void logger_set_output_level(Logger *logger, LogLevel level);
LogLevel logger_get_output_level(Logger *logger);

#endif /* LOGGER_H */
