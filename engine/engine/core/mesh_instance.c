#include "mesh_instance.h"

#include "strides.h"

#include "utils/logger.h"

#include <stdlib.h>
#include <string.h>

Status MeshInstance_init(MeshInstance* mi)
{
	memset(mi, 0, sizeof(MeshInstance));

	mi->scale = (V3){ 1,1,1 };
	mi->has_scale_changed = 1; // TODO: Rename to recalc bounding sphere?
    mi->texture_id = -1; // Default to untextured.

	return STATUS_OK;
}

Status MeshInstance_set_base(MeshInstance* mi, const MeshBase* mb)
{
	// Grow the vertex albedos buffer, there should be one albedo
	// per vertex.

    Vector_resize(mi->vertex_alebdos, mb->num_faces * STRIDE_FACE_VERTICES);

    if (!mi->vertex_alebdos.capacity == mb->num_faces * STRIDE_FACE_VERTICES)
    {
        return STATUS_ALLOC_FAILURE;
    }
	
	mi->mb_id = mb->id;

	// Default albedo to white.
	MeshInstance_set_albedo(mi, mb, (V3) { 1.f, 1.f, 1.f });
	
	return STATUS_OK;
}

void MeshInstance_set_albedo(MeshInstance* mi, const MeshBase* mb, V3 albedo)
{
	for (int i = 0; i < mb->num_faces * STRIDE_FACE_VERTICES; ++i)
	{
		mi->vertex_alebdos.data[i].x = albedo.x;
		mi->vertex_alebdos.data[i].y = albedo.y;
		mi->vertex_alebdos.data[i].z = albedo.z;
	}
}

void MeshInstance_destroy(MeshInstance* mi)
{
    if (!mi) return;

    Vector_destroy(mi->vertex_alebdos);
}