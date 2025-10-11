#ifndef TEXTURE_H
#define TEXTURE_H

#include "common/status.h"

#include <chds/vec.h>

#include <stdint.h>

// TODO: Refactor to simply be a canvas?? (e.g. remove texture?)
typedef struct {
	//float* data; // TODO: Should texture be uint32_t like canvas_t again? And we can still access each component?
    chds_vec(float) pixels;
	int width;
	int height;

} texture_t;

status_t texture_load_from_bmp(texture_t* texture, const char* file);

void texture_destroy(texture_t* texture);

#endif