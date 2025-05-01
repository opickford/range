#include "memory_utils.h"

#include "logger.h"

#include <stdlib.h>
#include <stdint.h>

Status resize_int_array(int** out_buffer, unsigned int len)
{
	// Realloc with size = 0 fails.
	if (len == 0)
	{
		return STATUS_INVALID_ARGUMENT;
	}

	int* temp_ptr = realloc(*out_buffer, len * sizeof(int));
	if (0 == temp_ptr)
	{
		log_error("Failed to resize_int_array.");
		return STATUS_ALLOC_FAILURE;
	}

	*out_buffer = temp_ptr;
	return STATUS_OK;
}

Status resize_float_array(float** out_buffer, unsigned int len)
{
	// Realloc with size = 0 fails.
	if (len == 0)
	{
		return STATUS_INVALID_ARGUMENT;
	}

	float* temp_ptr = realloc(*out_buffer, len * sizeof(float));
	if (0 == temp_ptr)
	{
		log_error("Failed to resize_float_array.");
		return STATUS_ALLOC_FAILURE;
	}

	*out_buffer = temp_ptr;
	return STATUS_OK;
}


Status resize_uint8_array(uint8_t** out_buffer, unsigned int len)
{
	// Realloc with size = 0 fails.
	if (len == 0)
	{
		return STATUS_INVALID_ARGUMENT;
	}

	uint8_t* temp_ptr = realloc(*out_buffer, len * sizeof(uint8_t));
	if (0 == temp_ptr)
	{
		log_error("Failed to resize_uint8_array.");
		return STATUS_ALLOC_FAILURE;
	}

	*out_buffer = temp_ptr;
	return STATUS_OK;
}

