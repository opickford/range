#ifndef MESH_BASE_H
#define MESH_BASE_H

#include "maths/vector3.h"

#include "common/status.h"

#include "utils/memory_utils.h"

#include <stdio.h>

/*

A MeshBase stores the static data of a Mesh. 
Currently, it essentially is a struct for a parsed .obj file.

A MeshInstance will store the ID of a MeshBase and read any static data from
there, therefore, if the base's data was edited, all instances would also be
updated.

It doesn't make sense to iterate over each MeshBase, therefore, there is no
need for use of DOD.

It also doesn't really make sense to remove a MeshBase, so the MeshBaseID will
always remain valid as an index without extra indirection logic.

*/

typedef int MeshBaseID;

typedef struct
{
	MeshBaseID id; // The index of the mesh base in it's MeshBases container.

	int num_faces;
	int num_positions;
	int num_normals;
	int num_uvs;

	int* position_indices;
	int* normal_indices;
	int* uv_indices;

	float* object_space_positions;
	float* object_space_normals;
	float* uvs; // TODO: Specifiy for textures?
	
	V3 centre;

} MeshBase;

typedef struct
{
	MeshBase* bases;
	int count;
} MeshBases;

// Helpers
void parse_obj_counts(FILE* file, int* num_positions, int* num_uvs, int* num_normals, int* num_faces);

// MeshBase API
Status mesh_base_init(MeshBase* mb);
Status mesh_base_from_obj(MeshBase* mb, const char* filename);
void mesh_base_destroy(MeshBase* mb);

// MeshBases API
Status mesh_bases_init(MeshBases* mbs);
MeshBaseID mesh_bases_add(MeshBases* mbs);
void mesh_bases_destroy(MeshBases* mbs);

#endif
 