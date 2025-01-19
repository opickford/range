#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <stdarg.h>

// TODO: Allow for formatting like with printf.
inline void log_info(const char* format, ...)
{
    va_list args;
    va_start(args, format);

    fprintf(stderr, "[INFO] : ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");

    va_end(args);
}


inline void log_warn(const char* format, ...)
{
    va_list args;
    va_start(args, format);

    fprintf(stderr, "[WARN] : ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");

    va_end(args);
}


inline void log_error(const char* format, ...)
{
    va_list args;
    va_start(args, format);

    fprintf(stderr, "[ERROR] : ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");

    va_end(args);
}

#endif