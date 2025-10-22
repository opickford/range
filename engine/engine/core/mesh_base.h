#ifndef MESH_BASE_H
#define MESH_BASE_H

#include "maths/vector3.h"
#include "maths/vector2.h"

#include "common/status.h"

#include <chds/vec.h>

#include <stdio.h>

/*

A mesh_base_t stores the static data of a Mesh. 
Currently, it essentially is a struct for a parsed .obj file.

A mesh_instance_t will store the ID of a mesh_base_t and read any static data from
there, therefore, if the base's data was edited, all instances would also be
updated.

It doesn't make sense to iterate over each mesh_base_t, therefore, there is no
need for use of DOD.

It also doesn't really make sense to remove a mesh_base_t, so the mesh_base_id_t will
always remain valid as an index without extra indirection logic.

*/

typedef int mesh_base_id_t;

typedef struct mesh_base
{
	mesh_base_id_t id; // The index of the mesh base in it's mesh_bases_t container.

	int num_faces;
	int num_positions;
	int num_normals;
	int num_uvs;

	chds_vec(int) position_indices;
	chds_vec(int) normal_indices;
    chds_vec(int) uv_indices;

    chds_vec(v3_t) object_space_positions;
	chds_vec(v3_t) object_space_normals;
	chds_vec(v2_t) uvs; // TODO: Specifiy for textures?
	
	v3_t centre;

} mesh_base_t;

// TODO: No longer need this, chds_vec(mesh_base_t)?
typedef struct
{
	mesh_base_t* bases;
	int count;
} mesh_bases_t;

// Helpers
void parse_obj_counts(FILE* file, int* num_positions, int* num_uvs, int* num_normals, int* num_faces);

// mesh_base_t API
status_t mesh_base_init(mesh_base_t* mb);
status_t mesh_base_from_obj(mesh_base_t* mb, const char* filename);
void mesh_base_destroy(mesh_base_t* mb);

// mesh_bases_t API
status_t mesh_bases_init(mesh_bases_t* mbs);
mesh_base_id_t mesh_bases_add(mesh_bases_t* mbs);
void mesh_bases_destroy(mesh_bases_t* mbs);

#endif
 