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

// TODO: Crefactor to component.

typedef int PointLightID;
typedef int ShadowCastingPointLightID;

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


// I think in terms of API for example, when creating shadow casting lights,
// we just want the user to call this, and this can create the shadow map?
//Status point_lights_add(PointLights* point_lights, PointLight point_light);
//Status point_lights_add_shadow_casting(PointLights* point_lights, ShadowCastingPointLight point_light);

Status point_light_init(PointLight* light);


#endif