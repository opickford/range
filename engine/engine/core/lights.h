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

// TODO: Should probably be split into different fiels.

// TODO: Crefactor to component.

typedef int PointLightID;

typedef struct
{
	V3 position;
	V3 colour;
	float strength;

} PointLight;

Status PointLight_init(PointLight* light);

typedef int ShadowCastingPointLightID;

typedef struct
{
	V3 position;
	V3 colour;
	float strength;

	DepthBuffer* depth_maps;

} ShadowCastingPointLight;

Status ShadowCastingPointLight_init(ShadowCastingPointLight* light);


#endif