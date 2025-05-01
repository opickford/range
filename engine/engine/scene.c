#include "scene.h"

#include "utils/logger.h"

#include <string.h>

Status scene_init(Scene* scene)
{
	memset(scene, 0, sizeof(Scene));

	mesh_bases_init(&scene->mesh_bases);
	mesh_instances_init(&scene->mesh_instances);
	lights_init(&scene->lights);

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

MeshInstanceID scene_add_mesh_instance(Scene* scene)
{
	return mesh_instances_add(&scene->mesh_instances);
}

Status scene_mesh_instance_set_base(Scene* scene, MeshInstanceID mi_id, MeshBaseID mb_id)
{
    // TODO: Helper for this?
    MeshInstance* mi = mesh_instances_get(&scene->mesh_instances, mi_id);
	const MeshBase* mb = &scene->mesh_bases.bases[mb_id];

	return mesh_instance_set_base(mi, mb);
}

Status scene_mesh_instance_set_albedo(Scene* scene, MeshInstanceID mi_id, V3 albedo)
{
    MeshInstance* mi = mesh_instances_get(&scene->mesh_instances, mi_id);
    const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
    mesh_instance_set_albedo(mi, mb, albedo);

    return STATUS_OK;
}
