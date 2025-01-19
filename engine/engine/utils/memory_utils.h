#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include "common/status.h"

// Helpers for resizing the models buffers.
Status resize_int_buffer(int** out_buffer, const unsigned int len);
Status resize_float_buffer(float** out_buffer, const unsigned int len);

#endif