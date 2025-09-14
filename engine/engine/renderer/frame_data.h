#ifndef FRAME_DATA_H
#define FRAME_DATA_H

#include "core/mesh_instance.h"
#include "core/scene.h"

#include "maths/bounding_sphere.h"

#include "common/status.h"

#include "utils/vector.h"

#include <cecs/ecs.h>

#include <stdint.h>

// Transient buffers used during the render pipeline per frame.

typedef struct
{
	// Transform Stage
	Vector(V3) view_space_positions;
    Vector(V3) view_space_normals;
    Vector(BoundingSphere) view_space_bounding_spheres; // Used for broad phase frustum culling.

    Vector(V3) point_lights_view_space_positions;

	// Broad Phase Frustum Culling
    Vector(MeshInstance) visible_mis;
	int num_visible_mis;
	Vector(uint8_t) intersected_planes;

	// Backface Culling Output
	Vector(int) front_face_indices;

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


    // TODO: Refactoring this so it's not just float arrays.
	
    Vector(V3) vertex_lighting;

	// Clipping
    Vector(float) faces_to_clip;  // Input to clip.
	Vector(float) clipped_faces;  // Clipping output.
	int num_clipped_faces; // Number of faces in the output buffer.

    // Intermediate buffers for clipping, alternate between per plane.
    Vector(float) temp_clipped_faces0; 
    Vector(float) temp_clipped_faces1;

} FrameData;

Status frame_data_init(
    ECS* ecs, 
    System* render_system, 
    System* lighting_system, 
    FrameData* frame_data, 
    Scene* scene);


void frame_data_destroy(FrameData* frame_data);


#endif