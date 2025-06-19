#include "lights.h"

#include "strides.h"

#include "renderer/render_buffers.h"

#include "utils/memory_utils.h"
#include "utils/logger.h"

#include <string.h>
#include <stdlib.h>

/*
void point_light_init(PointLight* light, V3 position, V3 colour, float strength)
{
	const int new_count = point_lights->count + 1;

	// Resize the view positions buffers.
	// TODO: I Would like to use STRIDE_POSITION without importing models.h.
	resize_float_array(&point_lights->world_space_positions, new_count * STRIDE_POSITION);
	resize_float_array(&point_lights->view_space_positions, new_count * STRIDE_POSITION);

	// Copy the lights position.
	v3_write(point_lights->world_space_positions + point_lights->count * STRIDE_POSITION, position);

	// Resize the point lights buffer.
	int old_size = point_lights->count * STRIDE_POINT_LIGHT_ATTRIBUTES;
	resize_float_array(&point_lights->attributes, old_size + STRIDE_POINT_LIGHT_ATTRIBUTES);

	// Copy the point light's attributes across.
	point_lights->attributes[old_size] = colour.x;
	point_lights->attributes[++old_size] = colour.y;
	point_lights->attributes[++old_size] = colour.z;
	point_lights->attributes[++old_size] = strength;

	// TODO: Create shadow maps.
	//		 All temporary for now.
	// TODO: No idea what size would be best.
	const int RES = 256;

	// TODO: Need cubemap for point light.
	// TODO: Potentially a depth buffer the same aspect ratio as the window could
	//		 give better results. Not sure just a thought.

	point_lights->count = new_count;

	// Create a new depth buffer for the new light.
	DepthBuffer* temp = realloc(point_lights->depth_maps, (size_t)point_lights->count * sizeof(DepthBuffer));
	if (!temp)
	{
		log_error("Failed to realloc for point_lights->depth_maps.");
		// TODO: Return status.
		return;
	}
	point_lights->depth_maps = temp;

	depth_buffer_init(&point_lights->depth_maps[point_lights->count - 1], RES, RES);
	depth_buffer_fill(&point_lights->depth_maps[point_lights->count - 1], 1.f);

	// Update the rbs lights count so we know to update the buffers.
	rbs->lights_count = point_lights->count;
}
*/

Status lights_init(Lights* lights)
{
	memset(lights, 0, sizeof(Lights));

    component_list_init(&lights->point_lights);

    return STATUS_OK;
}

void lights_destroy(Lights* lights)
{
	//free(lights->point_lights);
	free(lights->shadow_casting_point_lights);
}

Status point_light_init(PointLight* light)
{
    memset(light, 0, sizeof(PointLight));
    light->strength = 1.f;

    return STATUS_OK;
}

void PointLight_destroy(PointLight* light)
{
}

DEFINE_COMPONENT(PointLight)

