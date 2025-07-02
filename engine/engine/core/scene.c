#include "scene.h"

#include "utils/logger.h"

#include <string.h>

Status scene_init(Scene* scene)
{
	memset(scene, 0, sizeof(Scene));

	mesh_bases_init(&scene->mesh_bases);

	scene->ambient_light.x = 0.1f;
	scene->ambient_light.y = 0.1f;
	scene->ambient_light.z = 0.1f;

	return STATUS_OK;
}

Status scene_destroy(Scene* scene)
{
	// TODO: Cleanup
	log_warn("Need to implement scene_destroy.");
	return STATUS_OK;
}

MeshBaseID scene_add_mesh_base(Scene* scene)
{
	return mesh_bases_add(&scene->mesh_bases);
}
