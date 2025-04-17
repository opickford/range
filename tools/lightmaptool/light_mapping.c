#include "light_mapping.h"

#include <engine/engine.h>

// TODO: The utils/common are a little confusing maybe?
#include <engine/common/colour.h>
#include <engine/utils/common.h>
#include <engine/renderer/draw_2d.h>

#include <engine/maths/vector_maths.h>

// TODO: Surely we need to clamp the size yes?
// TODO: Using such a small area may cause issues because the gaps between the faces 
//		 on the lightmap may overlap? How can I sort this?
//#define TEXELS_PER_METER 16
#define TEXELS_PER_METER 128

/*
I could even implement ray tracing for the static lighting one day.
*/

// TODO: It would be nice to have an editor for everything... one day you may.

/*
TODO: UV Unwrapping.

- I thought I was getting there but actually we must implement some sort of unwrapping, even if it's terrible.

Must generate unique UVs for lightmapping



- All three vertices must be on a plane. This plane is then the 2D texture. Planar Projection?

Steps:
- Find normal of triangle.


*/


void do_light_mapping(Engine* engine)
{
	/*
	Steps:

	- Make shadow maps.

	- Find surface area of model

	- Create lightmap

	- Calculate uvs for model

	- Pack uvs into lightmap

	- Actually render to lightmap
	
	
	
	*/

	/*
	UV Unwrapping:

	Steps:
	- Either split mesh up into flatter parts
	  OR
	  Just per face. - i think this is fine for now.

	- Map 3D vertices to 2D UVs.

	- Pack UVs into atlas.
	
	
	
	*/




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


	float* gen_uvs = malloc(engine->scenes[0].models.mbs_faces_counts[0] * 6 * sizeof(float));
	if (!gen_uvs)
	{
		log_error("failed to alloc for gen_uvs.\n");
	}
	

	// TODO: Need to look into how this works, surely this has just made a copy of the entire scene.
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
		const int mb_index = models.mis_base_ids[mi];

		// SECTION: Lightmap canvas creation.
		
		// Calculate the whole surface area of the model so we can create a 
		// lightmap that is appropriately sized, we will have to clamp this at 
		// some point. TODO: Clamp size.
		float surface_area = 0.f;

		// For each position in the mi.
		for (int j = 0; j < models.mbs_faces_counts[mi]; ++j)
		{
			const int face_index = (models.mbs_faces_offsets[mi] + j) * STRIDE_FACE_VERTICES;

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
			const V3 u = v3_sub_v3(v1, v0);
			const V3 v = v3_sub_v3(v2, v0);

			// This is twice the size, but we can half later.
			surface_area += size(cross(u, v));
		}

		surface_area *= 0.5f;

		int lightmap_side_len = (int)(sqrtf(surface_area) * TEXELS_PER_METER);
			
		// Make the size a multiple of 4 so the bitmap is aligned properly.
		lightmap_side_len += lightmap_side_len % 4;

		Canvas lightmap;
		canvas_init(&lightmap, lightmap_side_len, lightmap_side_len);

		log_info("surface_area: %f", surface_area);
		log_info("lightmap_side_len: %d", lightmap_side_len);

		// END

		// Remember, view matrix is indentity here, so these are world space.
		const float* wsns = models.view_space_normals;


		typedef struct 
		{
			V2 v0, v1, v2;
		} UVTriangle;

		UVTriangle* uv_tris = malloc(models.mbs_faces_counts[mi] * sizeof(UVTriangle));

		// SECTION: UV Unwrapping
		for (int j = 0; j < models.mbs_faces_counts[mi]; ++j)
		{
			// Unpack vertices, again, TODO: Shouldn't have to do all of this.
			// precalculating will save performance and effort. Better to read the indices!
			int face_index = (models.mbs_faces_offsets[mi] + j) * STRIDE_FACE_VERTICES;

			const int index_v0 = models.mbs_face_position_indices[face_index] + positions_offset;
			const int index_v1 = models.mbs_face_position_indices[face_index + 1] + positions_offset;
			const int index_v2 = models.mbs_face_position_indices[face_index + 2] + positions_offset;

			const int index_parts_v0 = index_v0 * STRIDE_POSITION;
			const int index_parts_v1 = index_v1 * STRIDE_POSITION;
			const int index_parts_v2 = index_v2 * STRIDE_POSITION;

			V3 v0 = v3_read(wsps + index_parts_v0);
			V3 v1 = v3_read(wsps + index_parts_v1);
			V3 v2 = v3_read(wsps + index_parts_v2);
			
			// Find the normal of the face and the plane.
			const V3 normal = tri_normal(v0, v1, v2);

			// We now want to find two perpendicular vectors to form a 2D
			// coordinate system on the plane, the tangent and bitangent.

			// There are infinite number of tangents and bitangents right?
			// So for now, just pick randomly essentially. Although this 
			// might cause me problems.

			// Calculate a vector along the face to be the tangent.
			const V3 tangent = normalised(v3_sub_v3(v1, v0));
			
			// Cross product between normal and tangent should give bitangent.
			const V3 bitangent = normalised(cross(normal, tangent));

			// TODO: Is normalising the correct behaviour? Defintiely not. as we then lose the scale, we must normalise based on what? 
			//		 the TEXELS_PER_METER?

			V2 P_uv0 = {
				0, 0
			};

			V3 AB = v3_sub_v3(v1, v0);
			V2 P_uv1 = {
				dot(AB, tangent),
				dot(AB, bitangent)
			};


			V3 AC = v3_sub_v3(v2, v0);
			V2 P_uv2 = {
				dot(AC, tangent),
				dot(AC, bitangent)
			};

			printf("ab: %f\n", v2_size(v2_sub_v2(P_uv1, P_uv0)));
			printf("bc: %f\n", v2_size(v2_sub_v2(P_uv2, P_uv1)));
			printf("ca: %f\n", v2_size(v2_sub_v2(P_uv0, P_uv2)));

			//v2_normalise(&P_uv0);
			//v2_normalise(&P_uv1);
			//v2_normalise(&P_uv2);


			float scale = (float)lightmap_side_len / TEXELS_PER_METER;

			P_uv0.x /= scale;
			P_uv0.y /= scale;

			P_uv1.x /= scale;
			P_uv1.y /= scale;

			P_uv2.x /= scale;
			P_uv2.y /= scale;



			printf("P_uv0: %s\n", v2_to_str(P_uv0));
			printf("P_uv1: %s\n", v2_to_str(P_uv1));
			printf("P_uv2: %s\n", v2_to_str(P_uv2));


			/*
			int uv_index = j * 6;

			gen_uvs[uv_index] = P_uv0.x;
			gen_uvs[++uv_index] = P_uv0.y;
			gen_uvs[++uv_index] = P_uv1.x;
			gen_uvs[++uv_index] = P_uv1.y;
			gen_uvs[++uv_index] = P_uv2.x;
			gen_uvs[++uv_index] = P_uv2.y;
			*/

			uv_tris[j] = (UVTriangle){ P_uv0, P_uv1, P_uv2 };

			/*
			TODO:

			Gotta think about how we scale these uvs. 

			We get the UVs, but they should be scaled relative to the world space size of the light map in metres?

			*/


		}

		// SECTION: UV Packing
		/*
		TODO: I'm not so sure about how to do this to be honest.

		A simple way is to pack the largest triangles first into the map.
		but need to think about this stuff properly as probably not that easys

		I think the idea might be to essentially try and line up triangle edges.
		and rotate the triangles until they fit with the desired gap between?
		



		Pseudocode:
		sorted = sort triangles from largest to smallest

		# we now want to pair the triangles, obviously, if we take the largest, we can use it 
		# to create a square, then the second largest will fit inside of that square. and keep doing this.
		# NOTE: the issue would then be a large change in size, where multiple triangles could fit inside 
		#		an existing cell. - but for now this is fine i think.


		


		for tri in tris_sorted_largest_to_smallest:
			search for 

			assign a square/cell for the triangle to go into





		*/

		/*
		# TODO: Look at this seems better
		https://gamedev.stackexchange.com/questions/137656/triangle-texture-packing-problem

		Idea seems to be:

		- place biggest in middle 
		- next biggest line up largest edges
		- repeat for all triangles

		
		
		*/
		for (int j = 0; j < models.mbs_faces_counts[mi]; ++j)
		{
			// TODO: 
			UVTriangle tri = uv_tris[j];




		}














		// END

		// SECTION: Baking lightmaps

		

		// TODO: Now must draw each face of the triangle, but actually render in uv space given the uv
		//		 from blender.
		// For each world space position in the mi.

		for (int j = 0; j < models.mbs_faces_counts[mi]; ++j)
		{
			// TODO: This code each time I access a face is annoying, the indices will help,
			//		 should think about that properly.
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

			/*
			V2 lm_uv0 = v2_read(temp_uvs + j * 6);
			V2 lm_uv1 = v2_read(temp_uvs + j * 6 + 2);
			V2 lm_uv2 = v2_read(temp_uvs + j * 6 + 4);*/


			V2 lm_uv0 = v2_read(gen_uvs + j * 6);
			V2 lm_uv1 = v2_read(gen_uvs + j * 6 + 2);
			V2 lm_uv2 = v2_read(gen_uvs + j * 6 + 4);

			/*
			V2 lm_uv0 = { 0.581111, 0.581111 };
			V2 lm_uv1 = { 0.382149, 0.617851 };
			V2 lm_uv2 = { 0.418889, 0.581111 };
			*/
 
			printf("P_uv0: %s\n", v2_to_str(lm_uv0));
			printf("P_uv1: %s\n", v2_to_str(lm_uv1));
			printf("P_uv2: %s\n", v2_to_str(lm_uv2));
				 
				 

			V2 lm_uv0_pixels = v2_mul_f(lm_uv0, lightmap_side_len);
			V2 lm_uv1_pixels = v2_mul_f(lm_uv1, lightmap_side_len);
			V2 lm_uv2_pixels = v2_mul_f(lm_uv2, lightmap_side_len);

			// Draw the triangle with uv to draw to the triangle.
			const int min_x = (int)(min(lm_uv0_pixels.x, min(lm_uv1_pixels.x, lm_uv2_pixels.x)));
			const int min_y = (int)(min(lm_uv0_pixels.y, min(lm_uv1_pixels.y, lm_uv2_pixels.y)));
			const int max_x = (int)(max(lm_uv0_pixels.x, max(lm_uv1_pixels.x, lm_uv2_pixels.x)));
			const int max_y = (int)(max(lm_uv0_pixels.y, max(lm_uv1_pixels.y, lm_uv2_pixels.y)));

			// TODO: Rasterise in UV space, LERP world space coordinates.

			// Calculate the 'area' of the triangle, this is actually the double
			// signed area of the triangle, but it's fine as it's only used to normalise
			// the areas of the sub triangles for barycentric coordinates.
			float abc_area = double_signed_area(lm_uv0, lm_uv1, lm_uv2);

			printf("total_area: %f\n", abc_area);

			int colour = float_rgb_to_int(random_float(), random_float(), random_float());

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

			// Save the lightmap.
			const char* fmt = "C:\\Users\\olive\\source\\repos\\range\\res\\lightmaps\\out_%d.bmp";

			int length = snprintf(NULL, 0, fmt, j);
			char* out = malloc(length + 1);
			snprintf(out, length + 1, fmt, j);

			Status status = canvas_write_to_bitmap(&lightmap, out);
			if (STATUS_OK != status)
			{
				log_error("Failed to canvas_write_to_bitmap");
			}

			free(out);

			canvas_fill(&lightmap, 0);
			
		}
		/*
		// Save the lightmap.
		//char* out = "C:\\Users\\olive\\source\\repos\\range\\res\\lightmaps\\out.bmp";
		char* out = "C:\\Users\\olive\\source\\repos\\range\\res\\lightmaps\\";

		int length = snprintf(NULL, 0, "out_%d", j);
		char* num_str = malloc(length + 1);
		snprintf(num_str, length + 1, "out_%d", j);
		
			
		snprintf(str, size, "%d", x);


		strcat_s(out, strlen(out) + strlen(num_str)
		Status status = canvas_write_to_bitmap(&lightmap, );
		if (STATUS_OK != status)
		{
			log_error("Failed to canvas_write_to_bitmap");
		}

		free(num_str);
		*/
		// END

		// TODO: Again instead of this, should write these out to a buffer.
		positions_offset += models.mbs_positions_counts[mi];
		normals_offset += models.mbs_normals_counts[mi];
	}


	// TODO: Should lightmaps be packed together?
}