#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include "common/status.h"

#include <stdint.h >

// Helpers for resizing arrays.
Status resize_int_array(int** out_buffer, unsigned int len);
Status resize_float_array(float** out_buffer, unsigned int len);
Status resize_uint8_array(uint8_t** out_buffer, unsigned int len);

#endif