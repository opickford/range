#ifndef MeshInstance_H
#define MeshInstance_H

#include "mesh_base.h"
#include "bounding_sphere.h"

#include "maths/vector3.h"
#include "maths/vector4.h"

#include "common/status.h"

#include "utils/vector.h"

#include <stdint.h>


/*


TODO: Comments about how this is a component now.



*/



// TODO: Do we need this.
typedef int MeshInstanceID; 



typedef struct
{
	MeshBaseID mb_id; // TODO: Rename as base?
	int texture_id;

	uint8_t has_scale_changed; // Determines if the bounding sphere needs updating.

	// Per instance data
	// TODO: Do we actually want per vertex albedos? I reckon per face at least
	//		 makes more sense.
    // TODO: Float or V3???
	Vector(V3) vertex_alebdos;

	// Offsets into FrameData, these exist here as they are tied
	// to the MeshInstance itself.

    // TODO: When refactored to component would we have these? probs fine.
	int view_space_positions_offset;
	int view_space_normals_offset;

	int num_front_faces;

    BoundingSphere view_space_bounding_sphere; // Broad phase frustum culling.

} MeshInstance;

// MeshInstance API
Status MeshInstance_init(MeshInstance* mi, const MeshBase* mb);
void MeshInstance_destroy(MeshInstance* mi);

Status MeshInstance_set_base(MeshInstance* mi, const MeshBase* mb);
void MeshInstance_set_albedo(MeshInstance* mi, const MeshBase* mb, V3 albedo);

#endif