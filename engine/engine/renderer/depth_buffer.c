#include "depth_buffer.h"

#include "common/status.h"
#include "common/colour.h"

#include "utils/logger.h"

#include <Windows.h>

#include <stdlib.h>
#include <string.h>

status_t depth_buffer_init(depth_buffer_t* depth_buffer, int width, int height)
{
	memset(depth_buffer, 0, sizeof(depth_buffer_t));

	depth_buffer->width = width;
	depth_buffer->height = height;
	depth_buffer->data = malloc((size_t)width * height * sizeof(float));

	if (!depth_buffer->data)
	{
		log_error("Failed to allocate memory for depth_buffer data.");
		return STATUS_ALLOC_FAILURE;
	}

	return STATUS_OK;
}

status_t depth_buffer_rev3_size(depth_buffer_t* depth_buffer, int width, int height)
{
	// Check the size has changed.
	if (depth_buffer->width == width && depth_buffer->height == height)
	{
		return STATUS_OK;
	}

	// Allocate memory for the new array.
	// TODO: Use my memory allocating helpers for this.
	float* new_data = realloc(depth_buffer->data, (size_t)width * height * sizeof(float));

	// Check the allocation worked.
	if (!new_data)
	{
		log_error("Failed to reallocate memory for depth_buffer data on resize.");
		return STATUS_ALLOC_FAILURE;
	}

	// Update the canvas.
	depth_buffer->data = new_data;
	depth_buffer->width = width;
	depth_buffer->height = height;

	return STATUS_OK;
}

void depth_buffer_fill(depth_buffer_t* depth_buffer, float depth)
{
	// TODO: Look for some sort of blit or fill function 
	const int length = depth_buffer->width * depth_buffer->height;
	float* ptr = depth_buffer->data;

	unsigned int i = length;

	while (i)
	{
		*ptr = depth;
		--i;
		++ptr;
	}
}

void depth_buffer_draw(const depth_buffer_t* source, canvas_t* target, int x_offset, int y_offset)
{
    // TODO: What is this function????
    int* data = target->pixels;

	for (int y = 0; y < source->height; ++y)
	{
		for (int x = 0; x < source->height; ++x)
		{
			float depth = 1-source->data[y * source->width + x];
			int colour = float_rgb_to_int(depth, depth, depth);

			data[(y + y_offset) * target->width + x + x_offset] = colour;
		}
	}
}

void depth_buffer_destroy(depth_buffer_t* depth_buffer)
{
	free(depth_buffer->data);
	depth_buffer->data = 0;

	free(depth_buffer);
	depth_buffer = 0; // TODO: Do we want to do this?
}
