#include "light_mapping.h"

#include <engine/engine.h>

#include <engine/maths/vector_maths.h>

#define TEXELS_PER_METER 16

void do_light_mapping(Engine* engine)
{
	Scene scene = engine->scenes[engine->current_scene_id];
	Models models = scene.models;
	log_info("Generating lightmaps for %d mis.", scene.models.mis_count);
	
	// Convert all model space data to world space.
	M4 identity;
	m4_identity(identity);

	// By setting the view matrix to the identity, we're keeping all the 
	// data in world space.
	model_to_view_space(&models, identity);

	float* wsps = models.view_space_positions;
	
	int positions_offset = 0;

	// TODO: model_to_view_space should write out the positions offsets so
	//		 we can easily read the vsps of a given mi.
	for (int mi = 0; mi < models.mis_count; ++mi)
	{
		int mb_index = models.mis_base_ids[mi];

		float total_surface_area = 0.f;

		// For each position in the mi.
		for (int j = 0; j < models.mbs_faces_counts[mi]; ++j)
		{
			int face_index = (models.mbs_faces_offsets[mi] + j) * STRIDE_FACE_VERTICES;

			const int index_v0 = models.mbs_face_position_indices[face_index] + positions_offset;
			const int index_v1 = models.mbs_face_position_indices[face_index + 1] + positions_offset;
			const int index_v2 = models.mbs_face_position_indices[face_index + 2] + positions_offset;

			const int index_parts_v0 = index_v0 * STRIDE_POSITION;
			const int index_parts_v1 = index_v1 * STRIDE_POSITION;
			const int index_parts_v2 = index_v2 * STRIDE_POSITION;

			const V3 v0 = v3_read(wsps + index_parts_v0);
			const V3 v1 = v3_read(wsps + index_parts_v1);
			const V3 v2 = v3_read(wsps + index_parts_v2);

			
			// Calculate area of the triangle.
			// The magintude of the cross product equals the area of the parallelogram 
			// spanned by u and v. We can ignore the direction here.
			V3 u = v3_sub_v3(v1, v0);
			V3 v = v3_sub_v3(v2, v0);

			V3 dir = cross(u, v);
			
			// This is twice the size, but we can half later.
			//surface_area += ;
			float surface_area = sqrtf(size(dir) * 0.5f);
			int lightmap_size = (int)(surface_area * TEXELS_PER_METER);

			normalise(&dir);

			Canvas lightmap;
			canvas_init(&lightmap, lightmap_size, lightmap_size);
			

			log_info("surface_area: %f", surface_area);
			log_info("lightmap_size: %d", lightmap_size);
			
			// Calculate center vertex.
			V3 v3 = v3_mul_f(v3_add_v3(v0, v3_add_v3(v1, v2)), 1 / 3.f);

			V3 inv_dir = v3_mul_f(dir, -1);
			
			// TODO: How much do we need to move back to see the face? 
			//	     For now just try the surface area who knows...
			
			V3 from = v3_add_v3(v3, v3_mul_f(inv_dir, surface_area));

			printf("inv_dir %s\n", v3_to_str(inv_dir));

			printf("pos: test: %s\n", v3_to_str(v3_mul_f(from, -1)));

			M4 view_matrix;
			//m4_identity(view_matrix);


			V3 pos = v3_mul_f(from, -1);

			look_at(pos, inv_dir, view_matrix);

			printf("m4: %s\n", m4_to_str(view_matrix));

			// TODO: Find the correct one, this just ignores z,
			//		 could be fine.
			M4 othrographic = {
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 1
			};

			
			M4 combined;
			m4_mul_m4(othrographic, view_matrix, combined);
			
			
			V4 v0_proj, v1_proj, v2_proj;
			m4_mul_v4(combined, v3_to_v4(v0, 1.f), &v0_proj);
			m4_mul_v4(combined, v3_to_v4(v1, 1.f), &v1_proj);
			m4_mul_v4(combined, v3_to_v4(v2, 1.f), &v2_proj);


			printf("v0: %s\n", v4_to_str(v0_proj));
			printf("v1: %s\n", v4_to_str(v1_proj));
			printf("v2: %s\n", v4_to_str(v2_proj));
			

			
			
			/*
			TODO: My plan is to project each face onto it's own mini map and just combine them.
				  If performance is super bad I can probably just import from blender, but I dont
				  see why this would work too badly.


			Potential method of projection is:

			set view matrix look at dir in the opposite direction of the face, so we're looking directly at it.

			use orthographic projection matrix 
			
			TODO: ^^^

			Figure out the orthographic projection or find an alternative.

			TODO: Look at chatgpt most recent chat. LOOKS REALLY USEFUL.

			Thinking about it, the graphics im really aiming for are unturned. Most of the lighting is just ambient
			and looks fine, so I could definitely get away with that, it would be nice to have this though and i still
			do want shadows so static will have to be the way forward with that, because then if its night i can have
			cool little red/orangey lights for good atmosphere. Thinking about it I really don't need dynamic lights
			and don't care about them much unless we wanted a sun moving, but i think instead of that just lerp ambient lighting.


			*/





			printf("\n");
		}


		// Round the lightmap size to 4 TODO: is this necessary?
		//lightmap_size += lightmap_size % 4;

		

		// UV Unwrapping seems to be extremely complex

		// TODO: I'm going to do a very naive method of rendering each face as it's own mini lightmap and combining
		//		 all of these into 1 larger one per mesh.

		/*
		Surely UV unwrapping is as simple as finding faces with similar normals and writing them next to each ot
		*/







		positions_offset += models.mbs_positions_counts[mi];
	}




}