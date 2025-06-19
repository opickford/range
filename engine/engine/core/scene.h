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

// MeshBase API Wrappers
MeshBaseID scene_add_mesh_base(Scene* scene);

// MeshInstance API Wrappers
MeshInstanceID scene_add_mesh_instance(Scene* scene);

// TODO: These names are so long idk.
Status scene_mesh_instance_set_base(Scene* scene, MeshInstanceID mi_id, MeshBaseID mb_id);
Status scene_mesh_instance_set_albedo(Scene* scene, MeshInstanceID mi_id, V3 albedo);

#endif