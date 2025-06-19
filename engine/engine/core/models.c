#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS // Hide warnings for strok.
#endif

#include "models.h"

#include "renderer/render_buffers.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"

#include "common/status.h"

#include "utils/logger.h"
#include "utils/str_utils.h"
#include "utils/memory_utils.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void models_init(Models* models)
{
	// TODO: Comments
}

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

Status mb_from_obj(Models* models, RenderBuffers* rbs, const char* filename)
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
	int positions_count, normals_count, uvs_count, face_count;
	parse_obj_counts(file, &positions_count, &uvs_count, &normals_count, &face_count);

	const int mb_index = models->mbs_count;
	const int new_mbs_count = models->mbs_count + 1;

	// Write the number of vertex components in the mesh.
	resize_int_array(&models->mbs_positions_counts, new_mbs_count);
	models->mbs_positions_counts[mb_index] = positions_count;

	resize_int_array(&models->mbs_normals_counts, new_mbs_count);
	models->mbs_normals_counts[mb_index] = normals_count;

	resize_int_array(&models->mbs_faces_counts, new_mbs_count);
	models->mbs_faces_counts[mb_index] = face_count;

	resize_int_array(&models->mbs_uvs_counts, new_mbs_count);
	models->mbs_uvs_counts[mb_index] = uvs_count;

	// Resize the offset buffers.
	resize_int_array(&models->mbs_positions_offsets, new_mbs_count);
	models->mbs_positions_offsets[mb_index] = models->mbs_total_positions;

	resize_int_array(&models->mbs_normals_offsets, new_mbs_count);
	models->mbs_normals_offsets[mb_index] = models->mbs_total_normals;

	resize_int_array(&models->mbs_faces_offsets, new_mbs_count);
	models->mbs_faces_offsets[mb_index] = models->mbs_total_faces;

	resize_int_array(&models->mbs_uvs_offsets, new_mbs_count);
	models->mbs_uvs_offsets[mb_index] = models->mbs_total_uvs;

	// Resize the indexing buffers.
	const int new_total_faces = models->mbs_total_faces + face_count;
	const int new_total_vertices = new_total_faces * STRIDE_FACE_VERTICES;

	resize_int_array(&models->mbs_face_position_indices, new_total_vertices);
	resize_int_array(&models->mbs_face_normal_indices, new_total_vertices);
	resize_int_array(&models->mbs_face_uvs_indices, new_total_vertices);

	// Resize the actual object space data buffers.
	const int new_total_positions = models->mbs_total_positions + positions_count;
	const int new_total_normals = models->mbs_total_normals + normals_count;
	const int new_total_uvs = models->mbs_total_uvs + uvs_count;

	resize_float_array(&models->mbs_object_space_positions, new_total_positions * STRIDE_POSITION);
	resize_float_array(&models->mbs_object_space_normals, new_total_normals * STRIDE_NORMAL);
	resize_float_array(&models->mbs_uvs, new_total_uvs * STRIDE_UV);
	resize_float_array(&models->mbs_object_space_centres, new_mbs_count * STRIDE_POSITION);

	// Recreate the in/out buffers for clipping if this model base has the most
	// faces yet.
	if (models->max_mb_faces < face_count)
	{
		models->max_mb_faces = face_count;

		// Resize the shared render buffers.
		rbs->mbs_max_faces = models->max_mb_faces;
		render_buffers_resize(rbs);
	}

	// Move to the start of the file again so we can read it.
	rewind(file);

	// Define offsets to the start of the new mesh in the array.
	int positions_offset = models->mbs_total_positions * STRIDE_POSITION;
	int normals_offset = models->mbs_total_normals * STRIDE_NORMAL;
	int uvs_offset = models->mbs_total_uvs * STRIDE_UV;
	int faces_positions_offset = models->mbs_total_faces * STRIDE_FACE_VERTICES;
	int faces_normals_offset = models->mbs_total_faces * STRIDE_FACE_VERTICES;
	int faces_uvs_offset = models->mbs_total_faces * STRIDE_FACE_VERTICES;

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
			models->mbs_object_space_positions[positions_offset++] = v.x;
			models->mbs_object_space_positions[positions_offset++] = v.y;
			models->mbs_object_space_positions[positions_offset++] = v.z;
		}

		else if (strcmp(tokens[0], "vn") == 0)
		{
			models->mbs_object_space_normals[normals_offset++] = (float)atof(tokens[1]);
			models->mbs_object_space_normals[normals_offset++] = (float)atof(tokens[2]);
			models->mbs_object_space_normals[normals_offset++] = (float)atof(tokens[3]);
		}

		else if (strcmp(tokens[0], "vt") == 0)
		{	
			models->mbs_uvs[uvs_offset++] = (float)atof(tokens[1]);

			// Invert the v component so we don't have to do it per pixel.
			// Not sure if this will cause any confusion if we ever want to edit them.
			models->mbs_uvs[uvs_offset++] = 1.f - (float)atof(tokens[2]);
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

			models->mbs_face_position_indices[faces_positions_offset++] = face_indices[0];
			models->mbs_face_position_indices[faces_positions_offset++] = face_indices[3];
			models->mbs_face_position_indices[faces_positions_offset++] = face_indices[6];

			models->mbs_face_normal_indices[faces_normals_offset++] = face_indices[2];
			models->mbs_face_normal_indices[faces_normals_offset++] = face_indices[5];
			models->mbs_face_normal_indices[faces_normals_offset++] = face_indices[8];
			
			models->mbs_face_uvs_indices[faces_uvs_offset++] = face_indices[1];
			models->mbs_face_uvs_indices[faces_uvs_offset++] = face_indices[4];
			models->mbs_face_uvs_indices[faces_uvs_offset++] = face_indices[7];
		}
	}

	// Calculate the centre of the model base by taking the average of all the vertices.
	// After testing using indexed rendering for this, we got the wrong centre, works 
	// correctly with just averaging all the vertices, no matter how many times they're
	// used.
	V3 centre = { 0, 0, 0 };

	// Iterate through each vertex of each face in one go.
	int offset = models->mbs_positions_offsets[models->mbs_count];
	for (int i = 0; i < positions_count; ++i)
	{
		int index = (offset + i) * STRIDE_POSITION;

		V3 position = {
			models->mbs_object_space_positions[index],
			models->mbs_object_space_positions[++index],
			models->mbs_object_space_positions[++index]
		};
		v3_add_eq_v3(&centre, position);
	}

	v3_mul_eq_f(&centre, 1.0f / positions_count);

	int index_centre = models->mbs_count * STRIDE_POSITION;
	models->mbs_object_space_centres[index_centre] = centre.x;
	models->mbs_object_space_centres[++index_centre] = centre.y;
	models->mbs_object_space_centres[++index_centre] = centre.z;

	// Close the file.
	if (fclose(file) != 0)
	{
		log_error("Failed to close file after loading .obj file.");
		return STATUS_FILE_FAILURE;
	}

	// Update totals.
	models->mbs_count = new_mbs_count;
	models->mbs_total_faces = new_total_faces;
	models->mbs_total_positions = new_total_positions;
	models->mbs_total_normals = new_total_normals;
	models->mbs_total_uvs = new_total_uvs;

	return STATUS_OK;
}

void mi_create(Models* models, RenderBuffers* rbs, int mb_index, int n)
{
	// TODO: create_mis

	// TODO: Return status
	if (mb_index > models->mbs_count - 1)
	{
		log_error("mb_index out of range.");
		return; // TODO: return status.
	}

	const int new_instances_count = models->mis_count + n;

	// Resize buffers
	resize_int_array(&models->mis_base_ids, new_instances_count);
	resize_int_array(&models->mis_texture_ids, new_instances_count);
	resize_int_array(&models->mis_dirty_bounding_sphere_flags, new_instances_count);
	resize_int_array(&models->mis_intersected_planes, new_instances_count * 7); // 1 for how many planes, 6 for the potential plane indices. 
	resize_int_array(&models->mis_passed_broad_phase_flags, new_instances_count);

	for (int i = models->mis_count; i < new_instances_count; ++i)
	{
		models->mis_base_ids[i] = mb_index;

		// TODO: Textures. Should be a parameter?
		models->mis_texture_ids[i] = -1; // Default to no texture.

		models->mis_dirty_bounding_sphere_flags[i] = 1;
		models->mis_passed_broad_phase_flags[i] = 0;

		// TODO: define 7 as a stride.
		for (int j = i * 7; j < (i + 1) * 7; ++j)
		{
			models->mis_intersected_planes[j] = 0;
		}
	}

	resize_float_array(&models->mis_transforms, new_instances_count * STRIDE_MI_TRANSFORM);
	
	// Make space for the bounding sphere, this will be generated by the next render call.
	resize_float_array(&models->mis_bounding_spheres, new_instances_count * STRIDE_SPHERE);

	// Increase the totals to get the new buffer sizes.
	const int faces_count = models->mbs_faces_counts[mb_index] * n;
	const int old_total_faces = models->mis_total_faces;

	models->mis_total_faces += faces_count;
	models->mis_total_positions += models->mbs_positions_counts[mb_index] * n;
	models->mis_total_normals += models->mbs_normals_counts[mb_index] * n;

	// Initialise colours.
	const int new_total_vertices = models->mis_total_faces * STRIDE_FACE_VERTICES;
	resize_float_array(&models->mis_vertex_colours, new_total_vertices * STRIDE_COLOUR);

	// Initialise all vertex diffuse colours to white so they absorb all light.
	for (int i = old_total_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR; i < new_total_vertices * STRIDE_COLOUR; i += STRIDE_COLOUR)
	{
		models->mis_vertex_colours[i]     = 1.f;
		models->mis_vertex_colours[i + 1] = 1.f;
		models->mis_vertex_colours[i + 2] = 1.f;
	}

	// Resize buffers used for indexed rendering.
	resize_float_array(&models->view_space_positions, models->mis_total_positions * STRIDE_POSITION);
	resize_float_array(&models->view_space_normals, models->mis_total_normals * STRIDE_NORMAL);

	
	// Update the number of instances.
	models->mis_count = new_instances_count;
	
	// Update render buffer counts.
	rbs->instances_count = new_instances_count;
	rbs->total_faces = models->mis_total_faces;
}

void free_models(Models* models)
{
	// Free all mesh buffers.
	free(models->mbs_positions_counts);
	free(models->mbs_normals_counts);
	free(models->mbs_faces_counts);
	free(models->mbs_face_position_indices);
	free(models->mbs_face_normal_indices);
	free(models->mbs_object_space_positions);
	free(models->mbs_object_space_normals);
	free(models->mbs_uvs);

	free(models->mbs_faces_offsets);
	free(models->mbs_positions_offsets);
	free(models->mbs_normals_offsets);

	free(models->mis_base_ids);
	free(models->mis_texture_ids);
	free(models->mis_dirty_bounding_sphere_flags);
	free(models->mis_passed_broad_phase_flags);
	free(models->mis_intersected_planes);

	free(models->mis_vertex_colours);
	free(models->mis_transforms);
	free(models->mis_bounding_spheres);

	free(models->view_space_positions);
	free(models->view_space_normals);
}

void mi_set_transform(Models* models, int mi_index, V3 position, V3 eulers, V3 scale)
{
	const int ti = mi_index * STRIDE_MI_TRANSFORM;

	models->mis_transforms[ti + 0] = position.x;
	models->mis_transforms[ti + 1] = position.y;
	models->mis_transforms[ti + 2] = position.z;

	models->mis_transforms[ti + 3] = eulers.x;
	models->mis_transforms[ti + 4] = eulers.y;
	models->mis_transforms[ti + 5] = eulers.z;

	models->mis_transforms[ti + 6] = scale.x;
	models->mis_transforms[ti + 7] = scale.y;
	models->mis_transforms[ti + 8] = scale.z;
}

#undef _CRT_SECURE_NO_WARNINGS