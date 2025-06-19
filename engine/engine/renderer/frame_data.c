#include "frame_data.h"

#include "frustum_culling.h"

#include "core/scene.h"
#include "core/strides.h"

#include "utils/memory_utils.h"
#include "utils/logger.h"

#include <Windows.h> // max - TODO: Somewhere else?

Status frame_data_init(FrameData* frame_data, Scene* scene)
{
	// TODO: A flag to determine if the scene has changed and therefore we 
	//		 need to reinit the frame data buffer.

    // TODO: These will be resizing buffers down as well, do we want that??
    //       probably not, but then we need arenas.

	// Transform Stage
    resize_array(BoundingSphere, frame_data->view_space_bounding_spheres, scene->mesh_instances.count);
    
	int total_positions = 0;
	int total_normals = 0;
	int total_faces = 0;
    int most_faces = 0;
	for (int i = 0; i < scene->mesh_instances.count; ++i)
	{
		const MeshInstance* mi = &scene->mesh_instances.instances[i];
		const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
		total_positions += mb->num_positions;
		total_normals += mb->num_normals;
		total_faces += mb->num_faces;
        most_faces = max(most_faces, mb->num_faces);
	}

	resize_float_array(&frame_data->view_space_positions, total_positions * STRIDE_POSITION);
	
	// TODO: Fails in release?? heap corruption, so potentially from before?
	resize_float_array(&frame_data->view_space_normals, total_normals * STRIDE_NORMAL);

	resize_float_array(&frame_data->point_lights_view_space_positions, scene->lights.point_lights.count * STRIDE_POSITION);

	// Broad Phase Frustum Culling
	resize_int_array(&frame_data->visible_mi_indices, scene->mesh_instances.count);
	frame_data->num_visible_mis = 0;

	// Intersected planes stores the number of intersected planes and then the index of each plane,
	// so give room for the max planes + 1.
	resize_uint8_array(&frame_data->intersected_planes, scene->mesh_instances.count * (MAX_FRUSTUM_PLANES + 1));

	// Backface Culling Output
	resize_int_array(&frame_data->front_face_indices, total_faces);

	// Lighting
	// TODO: TEMP: For now no shadows so only 3 comps?
	resize_float_array(&frame_data->vertex_lighting, total_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR);

	// Clipping
    // x,y,z,r,g,b
	const int components_per_vertex = STRIDE_POSITION + STRIDE_COLOUR; // Temp hardcoded.

    const int max_tris_at_once = most_faces * MAX_CLIPPED_TRIS_FACTOR;

	frame_data->num_clipped_faces = 0;
	resize_float_array(&frame_data->faces_to_clip, components_per_vertex * total_faces * STRIDE_FACE_VERTICES);
	resize_float_array(&frame_data->clipped_faces, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);

    resize_float_array(&frame_data->temp_clipped_faces0, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);
    resize_float_array(&frame_data->temp_clipped_faces1, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);


	return STATUS_OK;
}

void frame_data_destroy(FrameData* frame_data)
{
}
