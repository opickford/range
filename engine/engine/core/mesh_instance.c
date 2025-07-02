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

	return STATUS_OK;
}

Status MeshInstance_set_base(MeshInstance* mi, const MeshBase* mb)
{
	// Grow the vertex albedos buffer, there should be one albedo
	// per vertex.
	Status status = resize_float_array(&mi->vertex_alebdos, mb->num_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR);

	if (STATUS_OK == status)
	{
		mi->mb_id = mb->id;

		// Default albedo to white.
		MeshInstance_set_albedo(mi, mb, (V3) { 1.f, 1.f, 1.f });
	}

	return status;
}

void MeshInstance_set_albedo(MeshInstance* mi, const MeshBase* mb, V3 albedo)
{
	for (int i = 0; i < mb->num_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR; i += STRIDE_COLOUR)
	{
		mi->vertex_alebdos[i] = albedo.x;
		mi->vertex_alebdos[i + 1] = albedo.y;
		mi->vertex_alebdos[i + 2] = albedo.z;
	}
}

void MeshInstance_destroy(MeshInstance* mi)
{
	free(mi->vertex_alebdos);
}