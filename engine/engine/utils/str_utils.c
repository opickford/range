#include "str_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

char* format_str(const char* format, ...)
{
    // TODO: I don't like this function at all....

    // TODO: Surelyyyy i dont need to use malloc here. Just define some buffer
    //       and return if requested string is too long.

    // Initialize the variable argument list.
    va_list args;
    va_start(args, format);

    // Calculate the length of the formatted string.
    int len = vsnprintf(NULL, 0, format, args) + 1;

    // Reset va_list to be used again
    va_end(args);
    va_start(args, format);

    // Allocate memory for the formatted string.
    char* buffer = malloc(len * sizeof(char));
    if (buffer == NULL)
    {
        printf("Memory allocation failed!\n");
        va_end(args);
        return NULL;
    }

    // Write the formatted string to the buffer.
    vsnprintf(buffer, len, format, args);

    // Clean up the variable argument list.
    va_end(args);

    return buffer;
}