#include "scene.h"

#include "utils/logger.h"

#include <string.h>

status_t scene_init(scene_t* scene)
{
	memset(scene, 0, sizeof(scene_t));

	mesh_bases_init(&scene->mesh_bases);

	scene->ambient_light.x = 0.1f;
	scene->ambient_light.y = 0.1f;
	scene->ambient_light.z = 0.1f;

	return STATUS_OK;
}

status_t scene_destroy(scene_t* scene)
{
	// TODO: Cleanup
	log_warn("Need to implement scene_destroy.");
	return STATUS_OK;
}

mesh_base_id_t scene_add_mesh_base(scene_t* scene)
{
	return mesh_bases_add(&scene->mesh_bases);
}
