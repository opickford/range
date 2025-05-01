#ifndef RENDER_BUFFERS_H
#define RENDER_BUFFERS_H

#include "strides.h"

#include "utils/memory_utils.h"

#include "common/status.h"

#include <string.h>
#include <math.h>

/*

Must think about this render buffers struct.

Maybe we start with the render pipeline though, think about the task as transformations on data.

so what data does each stage need.

if data is accessed together at the same time, should probably be structs.




*/


typedef struct
{
	// TODO: Eventually move intermediate buffers out of models to here.
	

	// Counts for helping with resizing.
	int mbs_max_faces;
	int lights_count; // TODO: Shadow casting lights only?
	int total_faces; // TODO: mi prefix or do we abstract that.
	int instances_count; // TODO: Same here ^^

	// This approach also means we don't need separate buffers per scene for 
	// clipping etc.

	// Backface culling buffers. // TODO: Redo comments.
	int* front_faces_counts;		// Number of faces that are visible to the camera.
	float* front_faces;				// An interleaved buffer of {x, y, z, u, v, x, y, z, r, g, b, r, g, b } for each vertex of each front face after backface culling.

	// Clipping buffers.
	float* temp_clipped_faces_in;
	float* temp_clipped_faces_out;
	float* clipped_faces;

	// Light space position buffers.
	float* light_space_positions; // Contains vertex positions in light space.
	float* front_face_light_space_positions;




	// Contains the data for rendering to triangles.
	//float* temp_light_space_positions;

	// Temporary buffer for drawing a triangle.
	float* triangle_vertices;

	

	float* light_space_pos_deltas;

	float* scanline_light_space_positions;
	float* scanline_light_space_pos_deltas;
	float* scanline_light_space_current_pos;
	
} RenderBuffers;

inline Status render_buffers_init(RenderBuffers* rbs)
{
	memset(rbs, 0, sizeof(RenderBuffers));

	return STATUS_OK;
}

inline Status render_buffers_resize(RenderBuffers* rbs)
{
	// TODO: Also this should just be done by a flag so at the start of the render,
	//       the buffers are resized. Or even we check each time.
	// TODO: TEMP: Resizing render buffer for storing light stuff.

	// Backface culling buffers.
	const int STRIDE_FRONT_FACE = STRIDE_BASE_FRONT_FACE + rbs->lights_count * STRIDE_V4 * STRIDE_FACE_VERTICES;
	resize_int_array(&rbs->front_faces_counts, rbs->instances_count);
	resize_float_array(&rbs->front_faces, rbs->total_faces * STRIDE_FRONT_FACE);

	// Clipping buffers.
	// TODO: Should the render buffers define this stride?
	// TODO: Surely this is gonna take too much memory. May have to 
	//		 refactor some of the rendering code.
	// Calculate the maximum number of triangles that one could turn into 
	// after clipping against all enabled planes.
	const int MAX_TRIS_FACTOR = (int)pow(2, 6);
	const int STRIDE_CLIPPED = STRIDE_BASE_CLIPPED_FACE + rbs->lights_count * STRIDE_V4 * STRIDE_FACE_VERTICES;
	const int max_clipped_tris = rbs->mbs_max_faces * MAX_TRIS_FACTOR * STRIDE_CLIPPED;

	resize_float_array(&rbs->temp_clipped_faces_in, max_clipped_tris);
	resize_float_array(&rbs->temp_clipped_faces_out, max_clipped_tris);
	resize_float_array(&rbs->clipped_faces, max_clipped_tris);

	// TODO: CALCULATE THE SIZE OF THE STRIDE PROPERLY?
	Status status = resize_float_array(&rbs->light_space_positions, rbs->total_faces * STRIDE_FACE_VERTICES * rbs->lights_count * STRIDE_V4); 
	resize_float_array(&rbs->front_face_light_space_positions, rbs->total_faces * STRIDE_FACE_VERTICES * rbs->lights_count * STRIDE_V4);

	// Pos (V4), UV (V2), albedo (V3), light (V3)

	// TODO: This will need to be different for when we do UVs.
	//const int vertex_components = (12 + rbs->lights_count * STRIDE_V4) * 4; // 4 vertices to allow for the split.

	// pos(v4), albedo(v3), light(v3)
	const int vertex_components = (10 + rbs->lights_count * STRIDE_V4) * 4; // 4 vertices to allow for the split.
	resize_float_array(&rbs->triangle_vertices, vertex_components);


	// Buffer for the light spaces for both edges.
	resize_float_array(&rbs->scanline_light_space_positions, rbs->lights_count * STRIDE_V4 * 2); // interleaved buffer of lsp0_v0, lsp0_v1, ... 
	resize_float_array(&rbs->light_space_pos_deltas, rbs->lights_count * STRIDE_V4 * 2); // interleaved buffer of dlsp0_dy_v0, dlsp0_dy_v1, ... 
	resize_float_array(&rbs->scanline_light_space_pos_deltas, rbs->lights_count * STRIDE_V4); 
	resize_float_array(&rbs->scanline_light_space_current_pos, rbs->lights_count * STRIDE_V4);
	
	return status;
}


// TODO: Cleanup

#endif