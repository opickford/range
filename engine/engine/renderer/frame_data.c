#include "frame_data.h"

#include "frustum_culling.h"

#include "core/scene.h"
#include "core/strides.h"
#include "core/components.h"

#include "utils/logger.h"

#include <cecs/ecs.h>

#include <Windows.h> // max - TODO: Somewhere else?

Status frame_data_init(
    ECS* ecs, 
    System* render_system,
    System* lighting_system,
    FrameData* frame_data, 
    Scene* scene)
{
    // TODO: A big issue with this system is that it will simply fail if 
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

    const MeshBase* mbs = scene->mesh_bases.bases;

    for (int si = 0; si < render_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = render_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        int mis_i = Archetype_find_component_list(archetype, COMPONENT_MeshInstance);
        MeshInstance* mis = archetype->component_lists[mis_i];

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            MeshInstance* mi = &mis[i];

            const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
            total_positions += mb->num_positions;
            total_normals += mb->num_normals;
            total_faces += mb->num_faces;
            most_faces = max(most_faces, mb->num_faces);

            ++mis_count;

        }
    }

	// Transform Stage
    Vector_reserve(frame_data->view_space_positions, total_positions * STRIDE_POSITION);
	
	// TODO: Fails in release?? heap corruption, so potentially from before?
    Vector_reserve(frame_data->view_space_normals, total_normals * STRIDE_NORMAL);

	// Broad Phase Frustum Culling
    Vector_reserve(frame_data->visible_mis, mis_count);
	frame_data->num_visible_mis = 0;

	// Intersected planes stores the number of intersected planes and then the index of each plane,
	// so give room for the max planes + 1.
    Vector_reserve(frame_data->intersected_planes, mis_count * (MAX_FRUSTUM_PLANES + 1));

	// Backface Culling Output
    Vector_reserve(frame_data->front_face_indices, total_faces);

	// Lighting
	// TODO: TEMP: For now no shadows so only 3 comps?
    Vector_reserve(frame_data->vertex_lighting, total_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR);


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
    Vector_reserve(frame_data->faces_to_clip, components_per_vertex * total_faces * STRIDE_FACE_VERTICES);
    Vector_reserve(frame_data->clipped_faces, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);

    Vector_reserve(frame_data->temp_clipped_faces0, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);
    Vector_reserve(frame_data->temp_clipped_faces1, components_per_vertex * max_tris_at_once * STRIDE_FACE_VERTICES);

    // Lighting system
    // TODO: I feel like a lot of this could be done with a define or a function
    int num_point_lights = 0;
    for (int si = 0; si < lighting_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = lighting_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        // TODO: If we introduce other lights this logic needs updating.
        num_point_lights += archetype->entity_count;
    }

    Vector_reserve(frame_data->point_lights_view_space_positions, num_point_lights * STRIDE_POSITION);

	return STATUS_OK;
}

void frame_data_destroy(FrameData* frame_data)
{
}
