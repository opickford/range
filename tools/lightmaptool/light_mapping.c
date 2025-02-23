#include "light_mapping.h"

// TODO: Use <> or "" ?
#include <engine/engine.h>

// TODO: The utils/common are a little confusing maybe?
#include <engine/common/colour.h>
#include <engine/utils/common.h>
#include <engine/renderer/draw_2d.h>

#include <engine/maths/vector_maths.h>

// TODO: Surely we need to clamp the size yes?
// TODO: Using such a small area may cause issues because the gaps between the faces 
//		 on the lightmap may overlap? How can I sort this?
#define TEXELS_PER_METER 16

/*
For now unwrapping the uvs may be too complex. Eventually I would like to probably implement 
some unwrapping myself as it seems quite interesting, or at least doing a naive method.

For now, just create unwrap uvs with 'lightmap pack' in blender, delete the inital uv map
then export and take the vts.

TODO: Not sure how I want to import these per model. - we should only ever
	  need two textures per model instance, so maybe just have a separate uv
	  channel that isn't optional. We can just load them separately?
*/



void do_light_mapping(Engine* engine)
{

	// TODO: Before all of this, should create shadow maps for
	//		 each light, should make them high resolution and
	//		 even can perform blurring when sampling them, must
	//		 remember that performance is not an issue here.
	//		 On this note, for point lights, literally just convert
	//		 it to 6 ones facing in the proper directions, again, not
	//		 an issue here.



	// Hardcoded cube lightmap uvs.
	// The UVs are wrote in blender in order for each face for each vertex, so I think here
	// we can just maintain this.

	float temp_uvs[] = {
		0.667668, 0.331331,
		0.997998, 0.001001,
		0.667668, 0.001001,
		0.998999, 0.002002,
		0.668669, 0.332332,
		0.998999, 0.332332,
		0.001001, 0.331331,
		0.331331, 0.001001,
		0.001001, 0.001001,
		0.332332, 0.335335,
		0.002002, 0.665666,
		0.332332, 0.665666,
		0.001001, 0.664665,
		0.331331, 0.334334,
		0.001001, 0.334334,
		0.332332, 0.002002,
		0.002002, 0.332332,
		0.332332, 0.332332,
		0.331331, 0.667668,
		0.001001, 0.667668,
		0.001001, 0.997998,
		0.002002, 0.998999,
		0.332332, 0.998999,
		0.332332, 0.668669,
		0.335335, 0.332332,
		0.665666, 0.332332,
		0.665666, 0.002002,
		0.664665, 0.334334,
		0.334334, 0.334334,
		0.334334, 0.664665,
		0.664665, 0.001001,
		0.334334, 0.001001,
		0.334334, 0.331331,
		0.335335, 0.665666,
		0.665666, 0.665666,
		0.665666, 0.335335 };

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
	int normals_offset = 0;

	// TODO: model_to_view_space should write out the positions offsets so
	//		 we can easily read the vsps of a given mi.
	for (int mi = 0; mi < models.mis_count; ++mi)
	{
		int mb_index = models.mis_base_ids[mi];

		float surface_area = 0.f;

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
			// The magnitude of the cross product equals the area of the parallelogram 
			// spanned by u and v. We can ignore the direction here.
			V3 u = v3_sub_v3(v1, v0);
			V3 v = v3_sub_v3(v2, v0);

			// This is twice the size, but we can half later.
			surface_area += size(cross(u, v));
		}

		surface_area *= 0.5f;


		int lightmap_side_len = (int)(sqrtf(surface_area) * TEXELS_PER_METER);

		lightmap_side_len += lightmap_side_len % 4;

		Canvas lightmap;
		canvas_init(&lightmap, lightmap_side_len, lightmap_side_len);

		log_info("surface_area: %f", surface_area);
		log_info("lightmap_side_len: %d", lightmap_side_len);

		// Remember, view matrix is indentity here, so these are world space.
		const float* wsns = models.view_space_normals;


		// TODO: Now must draw each face of the triangle, but actually render in uv space given the uv
		//		 from blender.
		// For each world space position in the mi.

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

			/* TODO: Normals...
			const int index_n0 = models.mbs_face_normal_indices[face_index] + normals_offset;
			const int index_n1 = models.mbs_face_normal_indices[face_index + 1] + normals_offset;
			const int index_n2 = models.mbs_face_normal_indices[face_index + 2] + normals_offset;

			int index_parts_n0 = index_n0 * STRIDE_NORMAL;
			int index_parts_n1 = index_n1 * STRIDE_NORMAL;
			int index_parts_n2 = index_n2 * STRIDE_NORMAL;

			const V3 n0 = v3_read(wsns + index_parts_n0);
			const V3 n1 = v3_read(wsns + index_parts_n1);
			const V3 n2 = v3_read(wsns + index_parts_n2);
			*/

			// TODO: Gotta think of how to attach the second channel of uvs to
			//		 my models. 
			//		 For now just dealing with the 1 cube so can just access the 
			//		 hardcoded ones.

			V2 lm_uv0 = v2_read(temp_uvs + j * 6);
			V2 lm_uv1 = v2_read(temp_uvs + j * 6 + 2);
			V2 lm_uv2 = v2_read(temp_uvs + j * 6 + 4);

			V2 lm_uv0_pixels = v2_mul_f(lm_uv0, lightmap_side_len);
			V2 lm_uv1_pixels = v2_mul_f(lm_uv1, lightmap_side_len);
			V2 lm_uv2_pixels = v2_mul_f(lm_uv2, lightmap_side_len);

			// Draw the triangle with uv to draw to the triangle.

			// TODO: For testing purposes, just render the triangles so hopefully
			//		 we get the same image as in blender.

			const int min_x = (int)(min(lm_uv0_pixels.x, min(lm_uv1_pixels.x, lm_uv2_pixels.x)));
			const int min_y = (int)(min(lm_uv0_pixels.y, min(lm_uv1_pixels.y, lm_uv2_pixels.y)));
			const int max_x = (int)(max(lm_uv0_pixels.x, max(lm_uv1_pixels.x, lm_uv2_pixels.x)));
			const int max_y = (int)(max(lm_uv0_pixels.y, max(lm_uv1_pixels.y, lm_uv2_pixels.y)));

			printf("%d %d %d %d\n", min_x, min_y, max_x, max_y);

			// TODO: Rasterise in UV space, LERP world space coordinates.

			// Calculate the 'area' of the triangle, this is actually the double
			// signed area of the triangle, but it's fine as it's only used to normalise
			// the areas of the sub triangles for barycentric coordinates.
			float abc_area = double_signed_area(lm_uv0, lm_uv1, lm_uv2);

			printf("total_area: %f\n", abc_area);

			int colour = float_rgb_to_int(random_float(), random_float(), random_float());

			/*
			draw_line(&lightmap, lm_uv0_pixels.x, lm_uv0_pixels.y, lm_uv1_pixels.x, lm_uv1_pixels.y, colour);
			draw_line(&lightmap, lm_uv0_pixels.x, lm_uv0_pixels.y, lm_uv2_pixels.x, lm_uv2_pixels.y, colour);
			draw_line(&lightmap, lm_uv1_pixels.x, lm_uv1_pixels.y, lm_uv2_pixels.x, lm_uv2_pixels.y, colour);
			*/
			// TODO: Outlines are correct, something going wrong with filling.
			//	     it seems like the triangles are overwritting each other, 

			// TODO: Use barycentric coordinates for triangle rasterization
			//		 to make the lerp code more simple, as again, perf is
			//		 not an issue or priority here.
			for (int y = min_y; y < max_y; ++y)
			{
				for (int x = min_x; x < max_x; ++x)
				{
					V2 P = { (float)x / lightmap_side_len, (float)y / lightmap_side_len };

					// U = area BCP / total area
					float u = double_signed_area(lm_uv1, lm_uv2, P) / abc_area;

					// V = area CAP / total area
					float v = double_signed_area(lm_uv2, lm_uv0, P) / abc_area;
					
					// u + v + w = 1.f if P in triangle.
					float w = 1.f - u - v;

					
					if (u < 0.f || u > 1.0f || v < 0.f || v > 1.0f || w < 0.f || w > 1.0f)
						continue;

					lightmap.pixels[y * lightmap.width + x] = colour;				
				}
			}

			// TODO: Implement phong shading (per pixel with ambient,
			//		 diffuse and specular).
		}

		// TODO: Is the behaviour where I can just keep writing intended?
		Status status = canvas_write_to_bitmap(&lightmap, "C:\\Users\\olive\\source\\repos\\range\\res\\lightmaps\\out.bmp");
		if (STATUS_OK != status)
		{
			log_error("Failed to canvas_write_to_bitmap");
		}

		positions_offset += models.mbs_positions_counts[mi];
		normals_offset += models.mbs_normals_counts[mi];
	}
}