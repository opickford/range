#ifndef DEPTH_BUFFER_H
#define DEPTH_BUFFER_H

// TODO: This is just a canvas that stores float data. Not sure of a better way around this.

#include "canvas.h"

#include "common/status.h"

typedef struct
{
	int width, height;
	float* data;

} DepthBuffer;

Status depth_buffer_init(DepthBuffer* depth_buffer, int width, int height);

Status depth_buffer_resize(DepthBuffer* depth_buffer, int width, int height);

void depth_buffer_fill(DepthBuffer* depth_buffer, float depth);

void depth_buffer_draw(const DepthBuffer* source, Canvas* target, int x_offset, int y_offset);

void depth_buffer_destroy(DepthBuffer* depth_buffer);


#endif