#ifndef LIGHTS_H
#define LIGHTS_H

#include "renderer/depth_buffer.h"

#include "maths/vector3.h"

#include "common/status.h"

// TODO: Should this just be a point_light.h file?

/*
A point light can be defined by a:
- Position
- Colour
- Strength

A vertex can reflect a certain amount of colour, this is it's colour,
therefore, if a white light shines on it, it will reflect more red if its
a red colour. 

*/

typedef struct
{
	V3 position;
	V3 colour;
	float strength;

} PointLight;

typedef struct
{
	V3 position;
	V3 colour;
	float strength;

	DepthBuffer* depth_maps;

} ShadowCastingPointLight;

typedef struct
{
	int num_point_lights;
	PointLight* point_lights;

	int num_shadow_casting_point_lights;
	ShadowCastingPointLight* shadow_casting_point_lights;

} Lights;

Status lights_init(Lights* lights);
void lights_destroy(Lights* lights);

//Status lights_init(PointLights* point_lights);
//void lights_destroy(PointLights* point_lights);

// I think in terms of API for example, when creating shadow casting lights,
// we just want the user to call this, and this can create the shadow map?
//Status point_lights_add(PointLights* point_lights, PointLight point_light);
//Status point_lights_add_shadow_casting(PointLights* point_lights, ShadowCastingPointLight point_light);

//Status point_light_init(PointLight* light);

// TODO: Destroy.


#endif