#ifndef CANVAS_H
#define CANVAS_H

#include "common/status.h"

#include <chds/vec.h>

#include <stdint.h>

typedef struct
{
	int width, height;
	chds_vec(uint32_t) pixels;

} canvas_t;

status_t canvas_init(canvas_t* canvas, int width, int height);

status_t canvas_rev3_size(canvas_t* canvas, int width, int height);

void canvas_fill(canvas_t* canvas, const unsigned int colour);

void canvas_draw(const canvas_t* source, canvas_t* target, int x_offset, int y_offset);

void canvas_destroy(canvas_t* canvas);

status_t canvas_init_from_bitmap(canvas_t* canvas, const char* file);

status_t canvas_write_to_bmp(const canvas_t* canvas, const char* file);

#endif