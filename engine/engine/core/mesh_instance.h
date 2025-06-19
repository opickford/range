#ifndef MESH_INSTANCE_H
#define MESH_INSTANCE_H

#include "mesh_base.h"

#include "maths/vector3.h"
#include "maths/vector4.h"

#include "common/status.h"

#include <stdint.h>

/*

A MeshInstance is an instance of a MeshBase. All static data will be read from
the base and the instance stores the non-shared data.

TODO: In the future, a MeshInstance could be an entity component but will have 
      to think about this more in the future.

For O(1) adding/removing, two maps are used.
- id_to_index maps a MeshInstanceID to the instance index, this allows us to 
    actually get the MeshInstance from the ID.
- index_to_id maps an instance index to it's MeshInstanceID, this allows us to 
    keep the instances array tightly packed by swapping the last element with 
    the one we're trying to remove.
*/

typedef int MeshInstanceID;

// TODO: This should be refactored to be a MeshComponent?
/*

Gotta think about it a bit.
*/


typedef struct
{
	MeshBaseID mb_id; // TODO: Rename as base?
	int texture_id;

	// Transform
	V3 position;
	V3 rotation;
	V3 scale;

	uint8_t has_scale_changed; // Determines if the bounding sphere needs updating.

	// Per instance data
	// TODO: Do we actually want per vertex albedos? I reckon per face at least
	//		 makes more sense.
	float* vertex_alebdos;

	// Offsets into FrameData, these exist here as they are tied
	// to the MeshInstance itself.

    // TODO: When refactored to component would we have these? probs fine.
	int view_space_positions_offset;
	int view_space_normals_offset;

	int num_front_faces;

} MeshInstance;

typedef struct
{
	int count;    // The number of valid MeshInstances
    int capacity; // The size of the arrays.

	MeshInstance* instances;

    // Map MeshInstanceIDs to indices in instances.
    int* id_to_index;
    MeshInstanceID* index_to_id;

    // Maintain an array of old ids and reuse when possible.
    int* free_ids;
    int free_ids_capacity;
    int free_ids_count;

} MeshInstances;

// MeshInstance API
Status mesh_instance_init(MeshInstance* mi);
void mesh_instance_destroy(MeshInstance* mi);

Status mesh_instance_set_base(MeshInstance* mi, const MeshBase* mb);
void mesh_instance_set_albedo(MeshInstance* mi, const MeshBase* mb, V3 albedo);

// MeshInstances API
Status mesh_instances_init(MeshInstances* mis);
void mesh_instances_destroy(MeshInstances* mis);

// This pointer is temporary, adding another instance will invalidate this 
// pointer, therefore, store a MeshInstanceID instead and only use the pointer
// for temporary convienience.
MeshInstance* mesh_instances_get(MeshInstances* mis, MeshInstanceID mi_id);

// We cannot store pointers to the MeshInstance because if we add another,
// the pointer could now be invalid after resizing the array.
MeshInstanceID mesh_instances_add(MeshInstances* mis);
void mesh_instances_remove(MeshInstances* mis, MeshInstanceID mi_id);

#endif