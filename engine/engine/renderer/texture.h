#ifndef TEXTURE_H
#define TEXTURE_H

#include "common/status.h"

#include <stdint.h>

typedef struct {
	float* data; // TODO: Should texture be uint32_t like Canvas again? And we can still access each component?
	int width;
	int height;

} Texture;

Status texture_load_from_bmp(Texture* texture, const char* file);

void texture_destroy(Texture* texture);

#endif