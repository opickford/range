#ifndef LIGHTS_H
#define LIGHTS_H

#include "renderer/render_buffers.h"
#include "renderer/depth_buffer.h"

#include "maths/vector3.h"

/*
A point light can be defined by a:
- Position
- Colour
- Strength

A vertex can reflect a certain amount of colour, this is it's colour,
therefore, if a white light shines on it, it will reflect more red if its
a red colour. 

*/

// TODO: Not all point lights should cast shadows.
 
typedef struct
{
	int count;

	// Data
	// Keep separate for when we convert from world space
	// to view space to avoid loading unncessary data into
	// the cache. Also as we cache the view space positions, 
	// we would not be accessing the world space positions again anyways.
	float* world_space_positions;
	float* attributes; // Strength and colour.
	
	// Cache the point light's view space position.
	float* view_space_positions; 

	DepthBuffer* depth_maps;

} PointLights;


void point_lights_init(PointLights* point_lights);

void point_lights_create(PointLights* point_lights, RenderBuffers* rbs, V3 position, V3 colour, float strength);

// TODO: Destroy.


#endif