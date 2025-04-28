#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include "common/status.h"

#include <stdint.h >

// Helpers for resizing the models buffers.
Status resize_int_buffer(int** out_buffer, unsigned int len);
Status resize_float_buffer(float** out_buffer, unsigned int len);
Status resize_uint8_buffer(uint8_t** out_buffer, unsigned int len);

#endif