#ifndef mesh_instance_H
#define mesh_instance_H

#include "mesh_base.h"
#include "bounding_sphere.h"

#include "maths/vector3.h"
#include "maths/vector4.h"

#include "common/status.h"

#include <chds/vec.h>

#include <stdint.h>


/*


TODO: Comments about how this is a component now.



*/



// TODO: Do we need this.
typedef int mesh_instance_id_t; 

typedef struct mesh_instance
{
	mesh_base_id_t mb_id; // TODO: Rename as base?
	int texture_id;

    // TODO: How do we ensure that this is set????
	uint8_t has_scale_changed; // Determines if the bounding sphere needs updating.

	// Per instance data
	// TODO: Do we actually want per vertex albedos? I reckon per face at least
	//		 makes more sense.
    // TODO: Float or v3_t???
	chds_vec(v3_t) vertex_alebdos;

    // TODO: A scene could have a base ambient light, but mesh instances should also be able to!
    //       This will let them glow. This could potentially be per vertex but probably not worth.

	// Offsets into frame_data_t, these exist here as they are tied
	// to the mesh_instance_t itself.

    // TODO: When refactored to component would we have these? probs fine.
	int view_space_positions_offset;
	int view_space_normals_offset;

	int num_front_faces;

    bounding_sphere_t view_space_bounding_sphere; // Broad phase frustum culling.

} mesh_instance_t;

// mesh_instance_t API
status_t mesh_instance_init(mesh_instance_t* mi, const mesh_base_t* mb);
void mesh_instance_destroy(mesh_instance_t* mi);

status_t mesh_instance_set_base(mesh_instance_t* mi, const mesh_base_t* mb);
void mesh_instance_set_albedo(mesh_instance_t* mi, const mesh_base_t* mb, v3_t albedo);

#endif