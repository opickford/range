#ifndef MESH_BASE_H
#define MESH_BASE_H

#include "maths/vector3.h"

#include "common/status.h"
#include "utils/memory_utils.h"

#include <stdio.h>

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
 