#ifndef FRAME_DATA_H
#define FRAME_DATA_H

#include "mesh_instance.h"
#include "scene.h"

#include "common/status.h"

#include <stdint.h>

/*
A collection of transient buffers used during the render pipeline per frame.
*/

// TODO: Allocate for all these buffers.

typedef struct
{
	// Transform Stage
	float* view_space_positions;
	float* view_space_normals;
	float* view_space_bounding_spheres; // Used for broad phase frustum culling.

	float* point_lights_view_space_positions;

	// Broad Phase Frustum Culling
	int* visible_mi_indices;
	int num_visible_mis;
	uint8_t* intersected_planes;

	// Backface Culling Output
	int* front_face_indices;

	// Lighting
	// TODO: Where should the vertex light output be written to?
	/*
	outputs:
	- total contribution from ambient and non shadow casting point lights - note,
	  this shouldnt be clamped because then we might subtract from it 
	- contribution (albedo * diffuse) from each shadow casting point light.
	 - then just subtract from the total?

	 just forget about shadows for now mayber.
	*/
	
	float* vertex_lighting;





	// Clipping
	float* faces_to_clip;  // Input to clip.
	float* clipped_faces;  // Clipping output.
	int num_clipped_faces; // Number of faces in the output buffer.

    // Intermediate buffers for clipping, alternate between per plane.
    float* temp_clipped_faces0; 
    float* temp_clipped_faces1;


} FrameData;

Status frame_data_init(FrameData* frame_data, Scene* scene);


void frame_data_destroy(FrameData* frame_data);


#endif