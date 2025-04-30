#include "mesh_instance.h"

#include "strides.h"

#include "utils/logger.h"

#include <stdlib.h>
#include <string.h>

Status mesh_instance_init(MeshInstance* mi)
{
	memset(mi, 0, sizeof(MeshInstance));

	mi->scale = (V3){ 1,1,1 };
	mi->has_scale_changed = 1; // TODO: Rename to recalc bounding sphere?

	return STATUS_OK;
}

Status mesh_instance_set_base(MeshInstance* mi, const MeshBase* mb)
{
	// Grow the vertex albedos buffer, there should be one albedo
	// per vertex.
	Status status = resize_float_buffer(&mi->vertex_alebdos, mb->num_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR);

	if (STATUS_OK == status)
	{
		mi->mb_id = mb->id;

		// Default albedo to white.
		mesh_instance_set_albedo(mi, mb, (V3) { 1.f, 1.f, 1.f });
	}

	return status;
}

void mesh_instance_set_albedo(MeshInstance* mi, const MeshBase* mb, V3 albedo)
{
	for (int i = 0; i < mb->num_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR; i += STRIDE_COLOUR)
	{
		mi->vertex_alebdos[i] = albedo.x;
		mi->vertex_alebdos[i + 1] = albedo.y;
		mi->vertex_alebdos[i + 2] = albedo.z;
	}
}

void mesh_instance_destroy(MeshInstance* mi)
{
	free(mi->vertex_alebdos);
}

Status mesh_instances_init(MeshInstances* mis)
{
	memset(mis, 0, sizeof(MeshInstances));

	return STATUS_OK;
}

void mesh_instances_destroy(MeshInstances* mis)
{
	free(mis->instances);
}

MeshInstanceID mesh_instances_add(MeshInstances* mis)
{
	// Create a new mesh instance and return a pointer to it.

	// Grow the array of mesh instances.
	const int new_count = mis->count + 1;
	MeshInstance* new_instances = realloc(mis->instances, new_count * sizeof(MeshInstance));

	if (!new_instances)
	{
		log_error("Failed to allocate for new mesh instance.");
		return 0;
	}

	const MeshInstanceID mi_id = mis->count;

	mis->instances = new_instances;
	mis->count = new_count;
	
	// Initialise the new mesh instance.
	MeshInstance* mi = &mis->instances[mi_id];
	mesh_instance_init(mi);

	return mi_id;
}
