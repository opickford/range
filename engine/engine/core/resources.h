#ifndef RESOURCES_H
#define RESOURCES_H

#include "renderer/texture.h"

#include "common/status.h"
#include "utils/logger.h"

#include <string.h>
#include <stdlib.h>

// TODO: Realistically I think these resources could be global.
// TODO: Or at least figure out a nicer way to pass them around.
#if 1
typedef struct
{
	texture_t* textures;
	int textures_count;

} resources_t;

inline void resources_init(resources_t* resources)
{
	memset(resources, 0, sizeof(resources_t));
}

// TODO: This could return some TextureID that can be used for indexing the texture.
inline status_t resources_load_texture(resources_t* resources, const char* file)
{
	// Get the texture's index.
	int i = resources->textures_count;

	// Make room for the new texture.
	texture_t* textures_temp = realloc(resources->textures, (size_t)(resources->textures_count + 1) * sizeof(texture_t));
	if (!textures_temp)
	{
		log_error("Failed to grow textures array in resources_load_texture because of %s.\n", status_to_str(STATUS_ALLOC_FAILURE));
		return STATUS_ALLOC_FAILURE;
	}
	resources->textures = textures_temp;

	// Try load the texture.
	status_t status = texture_load_from_bmp(&textures_temp[i], file);
	if (STATUS_OK != status)
	{
		log_error("Failed to texture_load_from_bmp in resources_load_texture because of %s.\n", status_to_str(status));
		return status;
	}

	// Successfully loaded the texture.
	++resources->textures_count;
	return STATUS_OK;
}

// TODO: Delete texture?
#else
typedef struct
{
	canvas_t* textures;
	int textures_count;

} resources_t;

inline void resources_init(resources_t* resources)
{
	memset(resources, 0, sizeof(resources_t));
}

inline status_t resources_load_texture(resources_t* resources, const char* file)
{
	// Get the texture's index.
	int i = resources->textures_count;

	// Make room for the new texture.
	canvas_t* textures_temp = realloc(resources->textures, (size_t)(resources->textures_count + 1) * sizeof(canvas_t));
	if (!textures_temp)
	{
		log_error("Failed to grow textures array in resources_load_texture because of %s.\n", status_to_str(STATUS_ALLOC_FAILURE));
		return STATUS_ALLOC_FAILURE;
	}
	resources->textures = textures_temp;

	// Try load the texture.
	status_t status = canvas_init_from_bitmap(&textures_temp[i], file);
	if (STATUS_OK != status)
	{
		log_error("Failed to canvas_init_from_bitmap in resources_load_texture because of %s.\n", status_to_str(status));
		return status;
	}

	// Successfully loaded the texture.
	++resources->textures_count;
	return STATUS_OK;
}

#endif

#endif