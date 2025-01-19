#ifndef SCENE_H
#define SCENE_H

#include "models.h"
#include "lights.h"

#include "maths/vector3.h"

#include "common/status.h"

// A scene is essentially a wrapper for models and lights.
typedef struct
{
	Models models;
	PointLights point_lights;
	V3 ambient_light;

} Scene;

/*
TODO:
We don't want Models here, we want entities. An entity will have a

The models are just for rendering and this should be abstracted from the game really.

Gotta think about how this is going to work, might get a bit more complicated.

I don't need to work on any more engine specific stuff for now. We have fully refactored the 
engine and all that. The more important thing now is getting the rendering working perfectly.

*/



Status scene_init(Scene* scene);
Status scene_destroy(Scene* scene);

#endif