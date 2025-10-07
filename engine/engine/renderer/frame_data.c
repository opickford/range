#include "frame_data.h"

#include "frustum_culling.h"

#include "core/scene.h"
#include "core/strides.h"
#include "core/components.h"

#include "utils/logger.h"

#include <cecs/ecs.h>

#include <Windows.h> // max - TODO: Somewhere else?

status_t frame_data_init(
    cecs_t* ecs, 
    cecs_view_id_t render_view,
    cecs_view_id_t lighting_view,
    frame_data_t* frame_data, 
    scene_t* scene)
{
    // TODO: A big issue with this view is that it will simply fail if 
    //       we don't have enough memory, honestly no clue how we could
    //       get around this.

	// TODO: A flag to determine if the scene has changed and therefore we 
	//		 need to reinit the frame data buffer.

    // TODO: Do we want to cache this data? Could be useful for debugging.

    int mis_count = 0;
    int total_positions = 0;
    int total_normals = 0;
    int total_faces = 0;
    int most_faces = 0;

    const mesh_base_t* mbs = scene->mesh_bases.bases;

    cecs_view_iter_t it = cecs_view_iter(ecs, render_view);
    while (cecs_view_iter_next(&it))
    {
        mesh_instance_t* mis = cecs_get_column(it, COMPONENT_MESH_INSTANCE);

        for (int i = 0; i < it.num_entities; ++i)
        {
            mesh_instance_t* mi = &mis[i];

            const mesh_base_t* mb = &scene->mesh_bases.bases[mi->mb_id];
            total_positions += mb->num_positions;
            total_normals += mb->num_normals;
            total_faces += mb->num_faces;
            most_faces = max(most_faces, mb->num_faces);

            ++mis_count;

        }
    }

	// transform_t Stage
    chds_vec_reserve(frame_data->view_space_positions, total_positions * STRIDE_POSITION);
	
	// TODO: Fails in release?? heap corruption, so potentially from before?
    chds_vec_reserve(frame_data->view_space_normals, total_normals * STRIDE_NORMAL);

	// Broad Phase Frustum Culling
    chds_vec_reserve(frame_data->visible_mis, mis_count);
	frame_data->num_visible_mis = 0;

	// Intersected planes stores the number of intersected planes and then the index of each plane,
	// so give room for the max planes + 1.
    chds_vec_reserve(frame_data->intersected_planes, mis_count * (MAX_FRUSTUM_PLANES + 1));

	// Backface Culling Output
    chds_vec_reserve(frame_data->front_face_indices, total_faces);

	// Lighting
	// TODO: TEMP: For now no shadows so only 3 comps?
    chds_vec_reserve(frame_data->vertex_lighting, total_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR);


    // TODO: I think I'll have to rethink this entire clipping process once we introduce dynamic lights again
    //       otherwise do we have to allocate space for each light for each mesh instance.

	// Clipping
    // x,y,z,r,g,b,u,v
    // UV are optional and may not be included but obviously we need the space incase.
	const int components_per_vertex = STRIDE_POSITION + STRIDE_COLOUR + STRIDE_UV; // Temp hardcoded.

    // TODO: I really don't like this tbf, not sure the extra memory is worth this, just process a clipped face at a time?
    //       IF! It becomes an issue.
    const int max_tris_at_once = most_faces * MAX_CLIPPED_TRIS_FACTOR;

    // TODO: More comments about this stuff and components per vertex etc.

	frame_data->num_clipped_faces = 0;
    chds_vec_reserve(frame_data->faces_to_clip, components_per_vertex * total_faces * STRIDE_FACE_VERTICES);
    chds_vec_reserve(frame_data->clipped_faces, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);

    chds_vec_reserve(frame_data->temp_clipped_faces0, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);
    chds_vec_reserve(frame_data->temp_clipped_faces1, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);

    // Lighting view
    int num_point_lights = 0;
    cecs_view_iter_t lighting_it = cecs_view_iter(ecs, lighting_view);
    while (cecs_view_iter_next(&lighting_it))
    {
        // TODO: If we introduce other lights this logic needs updating.
        num_point_lights += lighting_it.num_entities;
    }

    chds_vec_reserve(frame_data->point_lights_view_space_positions, num_point_lights * STRIDE_POSITION);

	return STATUS_OK;
}

void frame_data_destroy(frame_data_t* frame_data)
{
}
