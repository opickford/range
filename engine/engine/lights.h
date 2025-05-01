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

typedef struct
{
	int num_point_lights;
	PointLight* point_lights;

	int num_shadow_casting_point_lights;
	ShadowCastingPointLight* shadow_casting_point_lights;

	/*
	
	store an array of PointLightIDs to point light indices

	so if we add a new PointLight, expand both arrays by 1, return the length of the array -1.





	
	*/

} Lights;

Status lights_init(Lights* lights);
void lights_destroy(Lights* lights);

PointLightID lights_add_point_light(Lights* lights);
void lights_add_shadow_casting_point_light(Lights* lights);

// I think in terms of API for example, when creating shadow casting lights,
// we just want the user to call this, and this can create the shadow map?
//Status point_lights_add(PointLights* point_lights, PointLight point_light);
//Status point_lights_add_shadow_casting(PointLights* point_lights, ShadowCastingPointLight point_light);

//Status point_light_init(PointLight* light);

// TODO: Destroy.


#endif