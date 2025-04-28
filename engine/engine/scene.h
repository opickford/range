#ifndef SCENE_H
#define SCENE_H

#include "mesh_instance.h"
#include "lights.h"

#include "maths/vector3.h"

#include "common/status.h"

typedef struct
{
	// TODO: really this should be made up of entities.
	MeshBases mesh_bases;
	MeshInstances mesh_instances;

	Lights lights;
	V3 ambient_light;

} Scene;

Status scene_init(Scene* scene);
Status scene_destroy(Scene* scene);

#endif