#ifndef DEPTH_BUFFER_H
#define DEPTH_BUFFER_H

// TODO: This is just a canvas that stores float data. Not sure of a better way around this.

#include "core/canvas.h"

#include "common/status.h"

typedef struct
{
	int width, height;
	float* data;

} depth_buffer_t;

status_t depth_buffer_init(depth_buffer_t* depth_buffer, int width, int height);

status_t depth_buffer_rev3_size(depth_buffer_t* depth_buffer, int width, int height);

void depth_buffer_fill(depth_buffer_t* depth_buffer, float depth);

void depth_buffer_draw(const depth_buffer_t* source, canvas_t* target, int x_offset, int y_offset);

void depth_buffer_destroy(depth_buffer_t* depth_buffer);


#endif