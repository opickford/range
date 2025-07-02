#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include "common/status.h"

#include <stdint.h>

// Helpers for resizing arrays.
Status resize_int_array(int** out_buffer, unsigned int len);
Status resize_float_array(float** out_buffer, unsigned int len);
Status resize_uint8_array(uint8_t** out_buffer, unsigned int len);

// TODO: A vector implementation would probably be nicer.
// TODO: Not sure how to handle invalid args or realloc fail. - how to return status...
// TODO: Also confusing as we don't pass in **.
// TODO: I think resize array would be nicer if it called a function with a size of bytes.
#define resize_array(T, arr, len)                                              \
    do                                                                         \
    {                                                                          \
        if (len == 0) break;                                                   \
        T* temp = realloc((arr), (len) * sizeof(T));                           \
        if (0 == temp)                                                         \
        {                                                                      \
            log_error("Failed to resize_array for " #T);                       \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            (arr) = temp;                                                      \
        }                                                                      \
    } while (0)

#endif