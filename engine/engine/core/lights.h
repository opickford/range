#ifndef LIGHTS_H
#define LIGHTS_H

#include "renderer/depth_buffer.h"

#include "maths/vector3.h"

#include "common/status.h"

/*
A point light can be defined by a:
- Position
- Colour
- Strength

A vertex can reflect a certain amount of colour, this is it's colour,
therefore, if a white light shines on it, it will reflect more red if its
a red colour. 

*/

typedef int point_light_id_t;
typedef int shadow_casting_point_light_id_t;

// TODO: This having a position is a bit of an issue if we want to connect it to physics, do we want transform_t here?
typedef struct
{
	v3_t position;
	v3_t colour;
	float strength;

} point_light_t;


typedef struct
{
	v3_t position;
	v3_t colour;
	float strength;

	depth_buffer_t* depth_maps;

} shadow_casting_point_light_t;


// I think in terms of API for example, when creating shadow casting lights,
// we just want the user to call this, and this can create the shadow map?
//status_t point_lights_add(PointLights* point_lights, point_light_t point_light);
//status_t point_lights_add_shadow_casting(PointLights* point_lights, shadow_casting_point_light_t point_light);

status_t point_light_init(point_light_t* light);


#endif