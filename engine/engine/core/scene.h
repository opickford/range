#ifndef SCENE_H
#define SCENE_H

#include "mesh_instance.h"
#include "lights.h"

#include "maths/vector3.h"

#include "common/status.h"

// TODO: This should more be storing metadata of the scene, the ecs contains lights and instances etc.
typedef struct
{
	// TODO: really this should be made up of entities.
	MeshBases mesh_bases;

	V3 ambient_light;
    int bg_colour;

} Scene;

// TODO: refactor to remove the 

Status scene_init(Scene* scene);
Status scene_destroy(Scene* scene);

// MeshBase API Wrappers
MeshBaseID scene_add_mesh_base(Scene* scene);

#endif