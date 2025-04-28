#ifndef MESH_INSTANCE_H
#define MESH_INSTANCE_H

#include "mesh_base.h"

#include "maths/vector3.h"
#include "maths/vector4.h"

#include "common/status.h"

#include <stdint.h>

typedef struct
{
	MeshBase* mb; // TODO: Rename as base?
	int texture_id;

	// Transform
	V3 position;
	V3 rotation;
	V3 scale;

	uint8_t has_scale_changed; // Determines if the bounding sphere needs updating.

	// Per instance data
	// TODO: Do we actually want per vertex albedos? I reckon per face at least
	//		 makes more sense.
	// TODO: How do we alloc for this?
	float* vertex_alebdos; // Aligned with the mesh base positions

	// Offsets into FrameData, these exist here as they are tied
	// to the MeshInstance itself.
	int view_space_positions_offset;
	int view_space_normals_offset;

	int num_front_faces;

} MeshInstance;

typedef struct
{
	int count;
	MeshInstance* instances;
} MeshInstances;

// MeshInstance API
Status mesh_instance_init(MeshInstance* mi);
Status mesh_instance_set_base(MeshInstance* mi, MeshBase* mb);
void mesh_instance_set_albedo(MeshInstance* mi, V3 albedo);
void mesh_instance_destroy(MeshInstance* mi);

// MeshInstances API
Status mesh_instances_init(MeshInstances* mis);
void mesh_instances_destroy(MeshInstances* mis);

MeshInstance* mesh_instances_add(MeshInstances* mis);

#endif