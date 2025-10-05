#ifndef TEXTURE_H
#define TEXTURE_H

#include "common/status.h"

#include <chds/vector.h>

#include <stdint.h>

// TODO: Refactor to simply be a canvas?? (e.g. remove texture?)
typedef struct {
	//float* data; // TODO: Should texture be uint32_t like Canvas again? And we can still access each component?
    Vector(float) pixels;
	int width;
	int height;

} Texture;

Status texture_load_from_bmp(Texture* texture, const char* file);

void texture_destroy(Texture* texture);

#endif