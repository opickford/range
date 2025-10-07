#include "mesh_instance.h"

#include "strides.h"

#include "utils/logger.h"

#include <stdlib.h>
#include <string.h>

status_t mesh_instance_init(mesh_instance_t* mi, const mesh_base_t* mb)
{
	memset(mi, 0, sizeof(mesh_instance_t));

	mi->has_scale_changed = 1; // TODO: Rename to recalc bounding sphere?
    mi->texture_id = -1; // Default to untextured.

    status_t status = mesh_instance_set_base(mi, mb);

	return status;
}

status_t mesh_instance_set_base(mesh_instance_t* mi, const mesh_base_t* mb)
{
	// Grow the vertex albedos buffer, there should be one albedo
	// per vertex.
    chds_vec_resize(mi->vertex_alebdos, mb->num_faces * STRIDE_FACE_VERTICES);

    // TODO: Nicer way to check if vector resize succeeded?
    if (chds_vec_capacity(mi->vertex_alebdos) != mb->num_faces * STRIDE_FACE_VERTICES)
    {
        return STATUS_ALLOC_FAILURE;
    }
	
	mi->mb_id = mb->id;

	// Default albedo to white.
	mesh_instance_set_albedo(mi, mb, (v3_t) { 1.f, 1.f, 1.f });
	
	return STATUS_OK;
}

void mesh_instance_set_albedo(mesh_instance_t* mi, const mesh_base_t* mb, v3_t albedo)
{
	for (int i = 0; i < mb->num_faces * STRIDE_FACE_VERTICES; ++i)
	{
		mi->vertex_alebdos[i].x = albedo.x;
		mi->vertex_alebdos[i].y = albedo.y;
		mi->vertex_alebdos[i].z = albedo.z;
	}
}

void mesh_instance_destroy(mesh_instance_t* mi)
{
    if (!mi) return;

    chds_vec_destroy(mi->vertex_alebdos);
}