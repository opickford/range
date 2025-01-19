#include "memory_utils.h"

#include "logger.h"

#include <stdlib.h>

Status resize_int_buffer(int** out_buffer, const unsigned int len)
{
	// Realloc with size = 0 fails.
	if (len == 0)
	{
		return STATUS_INVALID_ARGUMENT;
	}

	int* temp_ptr = realloc(*out_buffer, len * sizeof(int));
	if (0 == temp_ptr)
	{
		log_error("Failed to resize_int_buffer.");
		return STATUS_ALLOC_FAILURE;
	}

	*out_buffer = temp_ptr;
	return STATUS_OK;
}

Status resize_float_buffer(float** out_buffer, const unsigned int len)
{
	// Realloc with size = 0 fails.
	if (len == 0)
	{
		return STATUS_INVALID_ARGUMENT;
	}

	float* temp_ptr = realloc(*out_buffer, len * sizeof(float));
	if (0 == temp_ptr)
	{
		log_error("Failed to resize_float_buffer.");
		return STATUS_ALLOC_FAILURE;
	}

	*out_buffer = temp_ptr;
	return STATUS_OK;
}