#include "mesh_base.h"

#include "maths/vector4.h"
#include "utils/logger.h"

#include "strides.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// Helpers
void parse_obj_counts(FILE* file, int* num_positions, int* num_uvs, int* num_normals, int* num_faces)
{
	// Find the number of vertex components, faces etc so we can define the memory in one go.
	char line[256];

	*num_positions = 0;
	*num_uvs = 0;
	*num_normals = 0;
	*num_faces = 0;

	while (fgets(line, sizeof(line), file))
	{
		char* token;

		token = strtok(line, " ");

		if (strcmp(token, "v") == 0)
		{
			++*num_positions;
		}
		else if (strcmp(token, "vn") == 0)
		{
			++*num_normals;
		}
		else if (strcmp(token, "vt") == 0)
		{
			++*num_uvs;
		}
		else if (strcmp(token, "f") == 0)
		{
			++*num_faces;
		}
	}
}

// MeshBase API
Status mesh_base_init(MeshBase* mb)
{
	memset(mb, 0, sizeof(MeshBase));
	return STATUS_OK;
}

Status mesh_base_from_obj(MeshBase* mb, const char* filename)
{
	// TODO: Eventually could check the filetype.
	FILE* file = fopen(filename, "r");

	if (NULL == file)
	{
		char* msg = format_str("Failed to open '%c' when loading .obj file.", filename);
		log_error(msg);
		free(msg);

		return STATUS_FILE_FAILURE;
	}

	// Read the sizes that we will need to allocate to accomodate for.
	parse_obj_counts(file, &mb->num_positions, &mb->num_uvs, &mb->num_normals, &mb->num_faces);

	// Allocate the buffers.
	resize_float_array(&mb->object_space_positions, mb->num_positions * STRIDE_POSITION);
	resize_float_array(&mb->object_space_normals, mb->num_normals * STRIDE_NORMAL);
	resize_float_array(&mb->uvs, mb->num_uvs * STRIDE_UV);

	const int num_vertices = mb->num_faces * STRIDE_FACE_VERTICES;
	resize_int_array(&mb->position_indices, num_vertices);
	resize_int_array(&mb->normal_indices, num_vertices);
	resize_int_array(&mb->uv_indices, num_vertices);


	/* TODO: ?
	// Recreate the in/out buffers for clipping if this model base has the most
	// faces yet.
	if (models->max_mb_faces < face_count)
	{
		models->max_mb_faces = face_count;

		// Resize the shared render buffers.
		rbs->mbs_max_faces = models->max_mb_faces;
		render_buffers_resize(rbs);
	}*/

	int positions_offset = 0;
	int normals_offset = 0;
	int uvs_offset = 0;

	int faces_positions_offset = 0;
	int faces_normals_offset = 0;
	int faces_uvs_offset = 0;

	// Move to the start of the file again so we can read it.
	rewind(file);

	// Fill the buffers from the file.
	char line[256];
	while (fgets(line, sizeof(line), file))
	{
		// Split the line into its tokens, we only need 4.
		char* tokens[4] = { "", "", "", "" };
		char* token = strtok(line, " ");

		int i = 0;
		while (token != NULL && i < 4)
		{
			tokens[i++] = token;
			token = strtok(NULL, " ");
		}

		if (strcmp(tokens[0], "v") == 0)
		{
			// Apply the model matrix to put the vertices in world space.
			V4 v = {
				(float)atof(tokens[1]),
				(float)atof(tokens[2]),
				(float)atof(tokens[3]),
				1
			};

			// Store the object space position.
			mb->object_space_positions[positions_offset++] = v.x;
			mb->object_space_positions[positions_offset++] = v.y;
			mb->object_space_positions[positions_offset++] = v.z;
		}

		else if (strcmp(tokens[0], "vn") == 0)
		{
			mb->object_space_normals[normals_offset++] = (float)atof(tokens[1]);
			mb->object_space_normals[normals_offset++] = (float)atof(tokens[2]);
			mb->object_space_normals[normals_offset++] = (float)atof(tokens[3]);
		}

		else if (strcmp(tokens[0], "vt") == 0)
		{
			mb->uvs[uvs_offset++] = (float)atof(tokens[1]);

			// Invert the v component so we don't have to do it per pixel.
			// Not sure if this will cause any confusion if we ever want to edit them.
			mb->uvs[uvs_offset++] = 1.f - (float)atof(tokens[2]);
		}

		else if (strcmp(tokens[0], "f") == 0)
		{
			// A face from the obj file is vertex index, uv index, normal index.
			int face_indices[9] = { 0 };

			// Define a buffer to store the string part of a face.
			char buffer[32] = "\0";

			// Faces must be triangulated
			for (int i = 0; i < 3; ++i)
			{
				// Break the string into per vertex indices.
				const char* vertex_str = tokens[1 + i];

				int char_index = 0;

				// For each vertex component, pos, uv, normal.
				for (int component_index = 0; component_index < 3; ++component_index)
				{
					// Overwrite the previous component.
					int buffer_index = 0;

					// Copy the component str into the buffer until we reach the 
					// delimiter or the end of the string.
					while (vertex_str[char_index] != '/' && vertex_str[char_index] != '\0')
					{
						buffer[buffer_index++] = vertex_str[char_index++];
					}

					// Write the null character to terminate the string.
					buffer[buffer_index] = '\0';

					// If the last character was a delimiter, move past it.
					if (vertex_str[char_index] == '/')
					{
						++char_index;
					}

					// Convert the buffer contents to an index.
					if (buffer_index == 0)
					{
						// 3 components per vertex.
						face_indices[i * 3 + component_index] = -1; // No index was defined.
					}
					else
					{
						face_indices[i * 3 + component_index] = atoi(buffer) - 1; // All indices are 1 based.
					}

					// Check for the end of the string.
					if (vertex_str[char_index] == '\0')
					{
						break;
					}
				}
			}

			mb->position_indices[faces_positions_offset++] = face_indices[0];
			mb->position_indices[faces_positions_offset++] = face_indices[3];
			mb->position_indices[faces_positions_offset++] = face_indices[6];

			mb->normal_indices[faces_normals_offset++] = face_indices[2];
			mb->normal_indices[faces_normals_offset++] = face_indices[5];
			mb->normal_indices[faces_normals_offset++] = face_indices[8];

			mb->uv_indices[faces_uvs_offset++] = face_indices[1];
			mb->uv_indices[faces_uvs_offset++] = face_indices[4];
			mb->uv_indices[faces_uvs_offset++] = face_indices[7];
		}
	}

	// Calculate the centre of the model base by taking the average of all the vertices.
	// After testing using indexed rendering for this, we got the wrong centre, works 
	// correctly with just averaging all the vertices, no matter how many times they're
	// used.
	mb->centre.x = 0;
	mb->centre.y = 0;
	mb->centre.z = 0;

	// Iterate through each vertex of each face in one go.
	for (int i = 0; i < mb->num_positions; ++i)
	{
		const int index = i * STRIDE_POSITION;

		const V3 position = {
			mb->object_space_positions[index],
			mb->object_space_positions[index + 1],
			mb->object_space_positions[index + 2]
		};
		v3_add_eq_v3(&mb->centre, position);
	}

	v3_mul_eq_f(&mb->centre, 1.0f / mb->num_positions);

	// Close the file.
	if (fclose(file) != 0)
	{
		log_error("Failed to close file after loading .obj file.");
		return STATUS_FILE_FAILURE;
	}

	return STATUS_OK;
}

void mesh_base_destroy(MeshBase* mb)
{
	free(mb->object_space_positions);
	free(mb->object_space_normals);
	free(mb->uvs);

	free(mb->position_indices);
	free(mb->normal_indices);
	free(mb->uv_indices);
}

// MeshBases API
Status mesh_bases_init(MeshBases* mbs)
{
	memset(mbs, 0, sizeof(MeshBases));
	return STATUS_OK;
}

MeshBaseID mesh_bases_add(MeshBases* mbs)
{
	// Grow the array of mesh instances.
	const int new_count = mbs->count + 1;
	MeshBase* new_bases = realloc(mbs->bases, new_count * sizeof(MeshBase));

	if (!new_bases)
	{
		log_error("Failed to allocate memory for new mesh base.");
		return STATUS_ALLOC_FAILURE;
	}

	const MeshBaseID mb_id = mbs->count;

	mbs->bases = new_bases;
	mbs->count = new_count;

	// Initialise the new MeshBase.
	MeshBase* mb = &mbs->bases[mb_id];
	mesh_base_init(mb);
	mb->id = mb_id;

	return mb_id;
}

void mesh_bases_destroy(MeshBases* mbs)
{
	free(mbs->bases);
}
