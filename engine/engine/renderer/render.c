#include "render.h"

#include "renderer.h"
#include "render_target.h"
#include "draw_2d.h"
#include "frame_data.h"
#include "frustum_culling.h"

#include "core/strides.h" // TODO: Surely this isn't a 'core' thing.
#include "core/globals.h"
#include "core/resources.h"
#include "core/components.h"
#include "core/transform.h"

#include "common/colour.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"
#include "maths/vector_maths.h"
#include "maths/utils.h"

#include "utils/timer.h"
#include "utils/common.h"

#include <stdio.h>
#include <string.h>


// TODO: Consider global render target.
// TODO: Think about my DOD. Should use structs rather than float* because will compile to the same thing.
// TODO: Something important would be to make the buffer accesses easier.

void debug_draw_point_lights(
    Canvas* canvas, 
    const ECS* ecs,
    const System* lighting_system,
    const FrameData* frame_data, 
    const RenderSettings* settings
    )
{
    const V3* vsps = frame_data->point_lights_view_space_positions.data;
    int vsps_offset = 0;

    // Debug draw point light icons as rects.
    for (int si = 0; si < lighting_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = lighting_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        int pls_i = Archetype_find_component_list(archetype, COMPONENT_PointLight);
        PointLight* pls = archetype->component_lists[pls_i];

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            V4 p = v3_to_v4(vsps[vsps_offset], 1.f);
            ++vsps_offset;

            // Only draw if depth is visibile in clip space.
            if (p.z > -settings->near_plane)
            {
                continue;
            }

            V4 projected;
            project(canvas, settings->projection_matrix, p, &projected);

            V3 colour_v3 = pls[i].colour;

            int colour = float_rgb_to_int(colour_v3.x, colour_v3.y, colour_v3.z);

            // Scale the radius so it's at a maximum of 10.
            const float radius = 10.f * (-settings->near_plane / p.z); // Square radius nice.

            int x0 = (int)(projected.x - radius);
            int x1 = (int)(projected.x + radius);

            int y0 = (int)(projected.y - radius);
            int y1 = (int)(projected.y + radius);

            draw_rect(canvas, x0, y0, x1, y1, colour);
        }
    }
}

void debug_draw_view_space_point(Canvas* canvas, const RenderSettings* settings, V3 point, int colour)
{
    // Convert from world space to screen space.
    V4 vsp = v3_to_v4(point, 1.f);

    // Don't draw points behind the camera.
    if (vsp.z > -settings->near_plane)
    {
        return;
    }

    V4 ssp;
    project(canvas, settings->projection_matrix, vsp, &ssp);

    int n = 2;
    int y0 = (int)(ssp.y - n);
    int y1 = (int)(ssp.y + n);
    int x0 = (int)(ssp.x - n);
    int x1 = (int)(ssp.x + n);

    draw_rect(canvas, x0, y0, x1, y1, colour);
}

void debug_draw_normals(Canvas* canvas, const FrameData* frame_data, const RenderSettings* settings, const Scene* scene)
{
    /*
    const MeshInstance* mis = scene->mesh_instances.instances;
    const MeshBase* mbs = scene->mesh_bases.bases;

    // TODO: Will need to calculate the actual stride of a vertex with the number of light space positions
    //		 too. Maybe to simplify this for all of my functions, should calculate this and keep it in the renderer.

    const int* visible_mi_indices = frame_data->visible_mi_indices;
    const int num_visible_mis = frame_data->num_visible_mis;

    const float* vsps = frame_data->view_space_positions;
    const float* vsns = frame_data->view_space_normals;
    const int* front_face_indices = frame_data->front_face_indices;

    int front_faces_offset = 0;

    for (int i = 0; i < num_visible_mis; ++i)
    {
        const MeshInstance* mi = &mis[visible_mi_indices[i]];
        const MeshBase* mb = &mbs[mi->mb_id];

        const int* position_indices = mb->position_indices;
        const int* normal_indices = mb->normal_indices;
        
        const int vsp_offset = mi->view_space_positions_offset;
        const int vsn_offset = mi->view_space_normals_offset;

        for (int j = 0; j < mi->num_front_faces; ++j)
        {
            const int face_index = front_face_indices[front_faces_offset + j] * STRIDE_FACE_VERTICES;
            for (int k = 0; k < STRIDE_FACE_VERTICES; ++k)
            {
                const int p_index = vsp_offset + position_indices[face_index + k] * STRIDE_POSITION;
                const int n_index = vsn_offset + normal_indices[face_index + k] * STRIDE_NORMAL;

                // TODO: Would just storing as a V3 be better?
                const V3 p = v3_read(vsps + p_index);
                const V3 n = v3_read(vsns + n_index);

                const float length = 2.5f;

                V3 end = v3_add_v3(p, v3_mul_f(n, length));

                V4 start_v4 = v3_to_v4(p, 1.f);
                V4 end_v4 = v3_to_v4(end, 1.f);

                V4 ss_start, ss_end;
                project(canvas, settings->projection_matrix, start_v4, &ss_start);
                project(canvas, settings->projection_matrix, end_v4, &ss_end);

                draw_line(canvas, (int)ss_start.x, (int)ss_start.y, (int)ss_end.x, (int)ss_end.y, COLOUR_LIME);
               
            }
        }
        front_faces_offset += mi->num_front_faces;
    }*/
}

/*
void debug_draw_bounding_spheres(Canvas* canvas, const RenderSettings* settings, const Models* models, const M4 view_matrix)
{
	// TODO: This doesn't really work.

	// TODO: For a function like this, I should be able to do debug_draw_bounding_sphere and pass in the mi index.

	int colour = int_rgb_to_int(1, 0, 0);

	for (int i = 0; i < models->mis_count; ++i)
	{
		int sphere_index = i * STRIDE_SPHERE;

		V3 view_centre_v3 = {
			models->mis_bounding_spheres[sphere_index],
			models->mis_bounding_spheres[sphere_index + 1],
			models->mis_bounding_spheres[sphere_index + 2]
		};

		debug_draw_view_space_point(canvas, settings, view_centre_v3, COLOUR_LIME);

		V4 view_centre = v3_to_v4(view_centre_v3, 1.f);

		V4 view_centre = m4_mul_v4(view_matrix, world_centre_v4);

		if (view_centre.z > -settings->near_plane)
		{
			continue;
		}

		float radius = models->mis_bounding_spheres[3];

		V4 world_bottom = world_centre_v4;
		world_bottom.y -= radius;

		V4 view_bottom = m4_mul_v4(view_matrix, world_bottom);

		V4 world_top = world_centre_v4;
		world_top.y += radius;

		V4 view_top = m4_mul_v4(view_matrix, world_top);

		V4 pc = project(canvas, settings->projection_matrix, view_centre);
		V4 pt = project(canvas, settings->projection_matrix, view_top);
		V4 pb = project(canvas, settings->projection_matrix, view_bottom);
		
		float pr = fabsf(pb.y - pt.y) / 2.f;
		
		draw_circle(canvas, (int)pc.x, (int)pc.y, (int)pr, colour);
	}
}*/

void debug_draw_world_space_point(Canvas* canvas, const RenderSettings* settings, V3 point, const M4 view_matrix, int colour)
{
	// Convert from world space to screen space.
	V4 wsp = v3_to_v4(point, 1.f);

	V4 vsp;
	m4_mul_v4(view_matrix, wsp, &vsp);
	
	// Don't draw points behind the camera.
	if (vsp.z > -settings->near_plane) 
	{
		return;
	}

	V4 ssp;
	project(canvas, settings->projection_matrix, vsp, &ssp);

	// TODO: Could be a draw 2d rect function.
	int n = 2;
	int y0 = (int)(ssp.y - n);
	int y1 = (int)(ssp.y + n);
	int x0 = (int)(ssp.x - n);
	int x1 = (int)(ssp.x + n);

	draw_rect(canvas, x0, y0, x1, y1, colour);
}

void debug_draw_world_space_line(Canvas* canvas, const RenderSettings* settings, const M4 view_matrix, V3 v0, V3 v1, V3 colour)
{
	V4 ws_v0 = v3_to_v4(v0, 1.f);
	V4 ws_v1 = v3_to_v4(v1, 1.f);

	V4 vs_v0, vs_v1;
	m4_mul_v4(view_matrix, ws_v0, &vs_v0);
	m4_mul_v4(view_matrix, ws_v1, &vs_v1);

	// Don't draw if behind the camera.
	if (vs_v0.z > -settings->near_plane && vs_v1.z > -settings->near_plane)
	{
		return;
	}

	// Clip the points to the near plane so we don't draw anything behind.
	// I think this works okay.
	if (vs_v0.z > -settings->near_plane)
	{
		float dist = vs_v0.z + settings->near_plane;
		
		float dxdz = (vs_v1.x - vs_v0.x) / (vs_v1.z - vs_v0.z);
		float dydz = (vs_v1.y - vs_v0.y) / (vs_v1.z - vs_v0.z);

		vs_v0.x += dxdz * dist;
		vs_v0.y += dydz * dist;
		vs_v0.z = -settings->near_plane;
	}
	else if (vs_v1.z > -settings->near_plane)
	{
		float dist = vs_v1.z + settings->near_plane;

		float dxdz = (vs_v0.x - vs_v1.x) / (vs_v0.z - vs_v1.z);
		float dydz = (vs_v0.y - vs_v1.y) / (vs_v0.z - vs_v1.z);

		vs_v1.x += dxdz * dist;
		vs_v1.y += dydz * dist;
		vs_v1.z = -settings->near_plane;
	}

	V4 ss_v0, ss_v1; 
	project(canvas, settings->projection_matrix, vs_v0, &ss_v0);
	project(canvas, settings->projection_matrix, vs_v1, &ss_v1);

	const int colour_int = float_rgb_to_int(colour.x, colour.y, colour.z);

	draw_line(canvas, (int)ss_v0.x, (int)ss_v0.y, (int)ss_v1.x, (int)ss_v1.y, colour_int);
}


void draw_scanline(RenderTarget* rt,
	int x0, int x1,
	int y,
	float z0, float z1,
	float w0, float w1,
	V3 start_colour, V3 end_colour)
{
	// TODO: Globals could be used for the render target pixels and depth buffer to make faster? Maybe? Would need to profile idk.


    if (y < 0 || y >= rt->canvas.height) log_error("invalid y %d \n", y);
    if (x0 < 0 || x0 >= rt->canvas.width) log_error("invalid x0 %d \n",x0);
    if (x1 < 0 || x1 >= rt->canvas.width) log_error("invalid x1 %d \n",x1);

	if (x0 == x1) return;

	// Precalculate deltas.
	const unsigned int dx = x1 - x0;
	float inv_dx = 1.f / dx;

	float w_step = (w1 - w0) * inv_dx;

	// Offset x by the given y.
	int row_offset = rt->canvas.width * y;

	int start_x = x0 + row_offset;
	int end_x = x1 + row_offset;

	// Render the scanline
	unsigned int* pixels = rt->canvas.pixels.data + start_x;
	float* depth_buffer = rt->depth_buffer + start_x;

	float inv_w = w0;
	float z = z0;
	float z_step = (z1 - z0) * inv_dx;

	
	V3 colour_step = v3_mul_f(v3_sub_v3(end_colour, start_colour), inv_dx);

	
	V3 colour = start_colour;

	for (unsigned int i = 0; i < dx; ++i)
	{
		// Depth test, only draw closer values.
		if (*depth_buffer > z)
		{
			// Recover w
			const float w = 1.0f / inv_w;

			// Calculate the colour of the vertex.
			
            *pixels = float_rgb_to_int(colour.x * w, colour.y * w, colour.z * w);
			
			*depth_buffer = z;			
		}

		// Move to the next pixel
		++pixels;
		++depth_buffer;

		// Step per pixel values.
		z += z_step;
		inv_w += w_step;

		v3_add_eq_v3(&colour, colour_step);
	}
}

void draw_flat_bottom_triangle(RenderTarget* rt, float* vc0, float* vc1, float* vc2)
{
	// Sort the flat vertices left to right.
	if (vc1[0] > vc2[0])
	{  
        SWAP(float*, vc1, vc2);
	}

	const V4 v0 = v4_read(vc0);
	const V4 v1 = v4_read(vc1);
	const V4 v2 = v4_read(vc2);

	const V3 c0 = v3_read(vc0 + 4);
	const V3 c1 = v3_read(vc1 + 4);
	const V3 c2 = v3_read(vc2 + 4);

	float inv_dy = 1 / (v2.y - v0.y);

	float dxdy0 = (v1.x - v0.x) * inv_dy;
	float dxdy1 = (v2.x - v0.x) * inv_dy;

	float dzdy0 = (v1.z - v0.z) * inv_dy;
	float dzdy1 = (v2.z - v0.z) * inv_dy;

	float dwdy0 = (v1.w - v0.w) * inv_dy;
	float dwdy1 = (v2.w - v0.w) * inv_dy;

	int start_y = (int)(ceil(v0.y - 0.5f));
	int end_y = (int)(ceil(v2.y - 0.5f));

	V3 dcdy0 = v3_mul_f(v3_sub_v3(c1, c0), inv_dy);
	V3 dcdy1 = v3_mul_f(v3_sub_v3(c2, c0), inv_dy);

	for (int y = start_y; y < end_y; ++y) {
		// Must lerp for the vertex attributes otherwise the accuracy is poor.
		// TODO: Would be nice to not have to actually lerp but step instead.
		float a = (y + 0.5f - v0.y);

		// Calculate the start and ends of the scanline
		float x0 = v0.x + dxdy0 * a;
		float x1 = v0.x + dxdy1 * a;

		float z0 = v0.z + dzdy0 * a;
		float z1 = v0.z + dzdy1 * a;

		float start_w = v0.w + dwdy0 * a;
		float end_w = v0.w + dwdy1 * a;

		int start_x = (int)((x0 - 0.5f));
		int end_x = (int)((x1 - 0.5f));

		V3 start_colour = v3_add_v3(c0, v3_mul_f(dcdy0, a));
		V3 end_colour = v3_add_v3(c0, v3_mul_f(dcdy1, a));
		
		draw_scanline(rt, start_x, end_x, y, z0, z1, start_w, end_w, start_colour, end_colour);
	}
}

void draw_flat_top_triangle(RenderTarget* rt, float* vc0, float* vc1, float* vc2)
{
	// Sort the flat vertices left to right.
	if (vc0[0] > vc1[0])
	{	
        SWAP(float*, vc0, vc1);
	}

	V4 v0 = v4_read(vc0);
	V4 v1 = v4_read(vc1);
	V4 v2 = v4_read(vc2);

	const V3 c0 = v3_read(vc0 + 4);
	const V3 c1 = v3_read(vc1 + 4);
	const V3 c2 = v3_read(vc2 + 4);

	float inv_dy = 1 / (v2.y - v0.y);

	float dxdy0 = (v2.x - v0.x) * inv_dy;
	float dxdy1 = (v2.x - v1.x) * inv_dy;

	float dzdy0 = (v2.z - v0.z) * inv_dy;
	float dzdy1 = (v2.z - v1.z) * inv_dy;

	float dwdy0 = (v2.w - v0.w) * inv_dy;
	float dwdy1 = (v2.w - v1.w) * inv_dy;

	int start_y = (int)(ceil(v0.y - 0.5f));
	int end_y = (int)(ceil(v2.y - 0.5f));

	V3 dcdy0 = v3_mul_f(v3_sub_v3(c2, c0), inv_dy);
	V3 dcdy1 = v3_mul_f(v3_sub_v3(c2, c1), inv_dy);
	
	for (int y = start_y; y < end_y; ++y) {
		// Must lerp for the vertex attributes to get them accurately.
		// TODO: Would be nice to find a way to step not lerp.
		//		 - not sure why stepping wouldn't work....

		// TODO: Step and not lerp.

		float a = (y + 0.5f - v0.y);

		float x0 = v0.x + dxdy0 * a;
		float x1 = v1.x + dxdy1 * a;

		float z0 = v0.z + dzdy0 * a;
		float z1 = v1.z + dzdy1 * a;

		float start_w = v0.w + dwdy0 * a;
		float end_w = v1.w + dwdy1 * a;

		int start_x = (int)((x0 - 0.5f));
		int end_x = (int)((x1 - 0.5f));

		V3 start_colour = v3_add_v3(c0, v3_mul_f(dcdy0, a));
		V3 end_colour = v3_add_v3(c1, v3_mul_f(dcdy1, a));

		draw_scanline(rt, start_x, end_x, y, z0, z1, start_w, end_w, start_colour, end_colour);
	}
}

void draw_triangle(RenderTarget* rt, float* vc0, float* vc1, float* vc2)
{
	// vc = vertex components

	// Sort vertices in ascending order.
	if (vc0[1] > vc1[1])
	{ 
        SWAP(float*, vc0, vc1);
	}
	if (vc0[1] > vc2[1])
	{
        SWAP(float*, vc0, vc2);
	}
	if (vc1[1] > vc2[1])
	{
        SWAP(float*, vc1, vc2);
	}
	
	// Handle if the triangle is already flat.
	if (vc0[1] == vc1[1])
	{
		draw_flat_top_triangle(rt, vc0, vc1, vc2);
		return;
	}

	if (vc1[1] == vc2[1])
	{
		draw_flat_bottom_triangle(rt, vc0, vc1, vc2);
		return;
	}
	
	// The triangle isn't flat, so split it into two flat triangles.
	float t = (vc1[1] - vc0[1]) / (vc2[1] - vc0[1]);
    
    // TODO: THis 7 should be defined as like BASE_VERTEX_COMPONENTS, x,y,z,w,r,g,b
    static float vc3[7] = { 0 };

	// Lerp for v3. 
	vc3[0] = vc0[0] + (vc2[0] - vc0[0]) * t;
	vc3[1] = vc1[1];
	vc3[2] = vc0[2] + (vc2[2] - vc0[2]) * t;
	vc3[3] = vc0[3] + (vc2[3] - vc0[3]) * t;

	
    // TODO: This should be defined somewhere!!!!!!!!!

	// Lerp the components / attributes?.
	// TODO: Could be nice to have this as a stride?
	// albedo, diffuse, light space pos * count
	const int STRIDE_COMPONENTS = STRIDE_COLOUR;
	//const int STRIDE_COMPONENTS = STRIDE_COLOUR + STRIDE_COLOUR + lights_count * STRIDE_V4;
	
	// Lerp for each component.
	for (int i = 4; i < 4 + STRIDE_COMPONENTS; ++i)
	{	
		vc3[i] = vc0[i] + (vc2[i] - vc0[i]) * t;
	}

	draw_flat_top_triangle(rt, vc1, vc3, vc2);
	draw_flat_bottom_triangle(rt, vc0, vc1, vc3);
}

void draw_textured_scanline(RenderTarget* rt, int x0, int x1, int y, float z0, float z1, float w0, float w1, const V3 c0, const V3 c1, const V2 uv0, const V2 uv1, const Texture* texture)
{
	// TODO: Refactor function args.
	if (x0 == x1) return;

	// Precalculate deltas.
	const unsigned int dx = x1 - x0;
	float inv_dx = 1.f / dx;

	float w_step = (w1 - w0) * inv_dx;

	// Offset x by the given y.
	int row_offset = rt->canvas.width * y;

	int start_x = x0 + row_offset;
	int end_x = x1 + row_offset;

	V3 c_step = v3_mul_f(v3_sub_v3(c1, c0), inv_dx);
	V3 c = c0;

	v3_mul_eq_f(&c_step, 255);
	v3_mul_eq_f(&c, 255);

	V2 uv = uv0;

	V2 uv_step = 
	{
		(uv1.x - uv0.x) * inv_dx,
		(uv1.y - uv0.y) * inv_dx
	};


	// Render the scanline
	const float* texture_data = texture->pixels.data;

	// Write out like this to avoid packing the int back together.
	uint8_t* rgbas = (uint8_t*)(rt->canvas.pixels.data + start_x);

	float* depth_buffer = rt->depth_buffer + start_x;

	float inv_w = w0;
	float z = z0;
	float z_step = (z1 - z0) * inv_dx;

	const int texture_width = texture->width;
	
	for (unsigned int i = 0; i < dx; ++i)
	{
		// Depth test, only draw closer values.
		if (*depth_buffer > z)
		{
			// Recover w
			const float w = 1.0f / inv_w;

			// The texture width/height scaling is applied
			// before drawing the triangle to reduce calculations.
			int cols = (int)((uv.x * w));
			int rows = (int)((uv.y * w));

			int n = (rows * texture_width + cols) * 3; // Texture is split into float r,g,b components.
			
			float r = texture_data[n] * c.x * w;
			float g = texture_data[n + 1] * c.y * w;
			float b = texture_data[n + 2] * c.z * w;
			
            //float r = texture_data[n];
            //float g = texture_data[n + 1];
            //float b = texture_data[n + 2];

            //printf("%f %f %f\n", r, g, b);

			*rgbas = b;
			*(rgbas + 1) = g;
			*(rgbas + 2) = r;
			//*(rgbas + 3) = 0; // It's not actually necessary to set this.
			
			*depth_buffer = z;
		}

		// Move to the next pixel
		rgbas += 4;
		++depth_buffer;

		// Step per pixel values.
		z += z_step;
		inv_w += w_step;

		uv.x += uv_step.x;
		uv.y += uv_step.y;

		v3_add_eq_v3(&c, c_step);
	}
}

void draw_textured_flat_bottom_triangle(RenderTarget* rt, float* vc0, float* vc1, float* vc2, const Texture* texture)
{
    // Sort the flat vertices left to right.
    if (vc1[0] > vc2[0])
    {
        SWAP(float*, vc1, vc2);
    }

    const V4 v0 = v4_read(vc0);
    const V4 v1 = v4_read(vc1);
    const V4 v2 = v4_read(vc2);

    const V3 c0 = v3_read(vc0 + 4);
    const V3 c1 = v3_read(vc1 + 4);
    const V3 c2 = v3_read(vc2 + 4);

    const V2 uv0 = v2_read(vc0 + 7);
    const V2 uv1 = v2_read(vc1 + 7);
    const V2 uv2 = v2_read(vc2 + 7);

    float inv_dy = 1 / (v2.y - v0.y);

    float dxdy0 = (v1.x - v0.x) * inv_dy;
    float dxdy1 = (v2.x - v0.x) * inv_dy;

    float dzdy0 = (v1.z - v0.z) * inv_dy;
    float dzdy1 = (v2.z - v0.z) * inv_dy;

    float dwdy0 = (v1.w - v0.w) * inv_dy;
    float dwdy1 = (v2.w - v0.w) * inv_dy;

    int start_y = (int)(ceil(v0.y - 0.5f));
    int end_y = (int)(ceil(v2.y - 0.5f));

    V3 dcdy0 = v3_mul_f(v3_sub_v3(c1, c0), inv_dy);
    V3 dcdy1 = v3_mul_f(v3_sub_v3(c2, c0), inv_dy);

    V2 duvdy0 = v2_mul_f(v2_sub_v2(uv1, uv0), inv_dy);
    V2 duvdy1 = v2_mul_f(v2_sub_v2(uv2, uv0), inv_dy);

    for (int y = start_y; y < end_y; ++y) {
        // Must lerp for the vertex attributes otherwise the accuracy is poor.
        // TODO: Would be nice to not have to actually lerp but step instead.
        float a = (y + 0.5f - v0.y);

        // Calculate the start and ends of the scanline
        float x0 = v0.x + dxdy0 * a;
        float x1 = v0.x + dxdy1 * a;

        float z0 = v0.z + dzdy0 * a;
        float z1 = v0.z + dzdy1 * a;

        float start_w = v0.w + dwdy0 * a;
        float end_w = v0.w + dwdy1 * a;

        int start_x = (int)((x0 - 0.5f));
        int end_x = (int)((x1 - 0.5f));

        V3 start_colour = v3_add_v3(c0, v3_mul_f(dcdy0, a));
        V3 end_colour = v3_add_v3(c0, v3_mul_f(dcdy1, a));

        V2 start_uv = v2_add_v2(uv0, v2_mul_f(duvdy0, a));
        V2 end_uv = v2_add_v2(uv0, v2_mul_f(duvdy1, a));

        draw_textured_scanline(rt, start_x, end_x, y, z0, z1, start_w, end_w, start_colour, end_colour, start_uv, end_uv, texture);
    }
}

void draw_textured_flat_top_triangle(RenderTarget* rt, float* vc0, float* vc1, float* vc2, const Texture* texture)
{
    // Sort the flat vertices left to right.
    if (vc0[0] > vc1[0])
    {
        SWAP(float*, vc0, vc1);
    }

    V4 v0 = v4_read(vc0);
    V4 v1 = v4_read(vc1);
    V4 v2 = v4_read(vc2);

    const V3 c0 = v3_read(vc0 + 4);
    const V3 c1 = v3_read(vc1 + 4);
    const V3 c2 = v3_read(vc2 + 4);

    const V2 uv0 = v2_read(vc0 + 7);
    const V2 uv1 = v2_read(vc1 + 7);
    const V2 uv2 = v2_read(vc2 + 7);

    float inv_dy = 1 / (v2.y - v0.y);

    float dxdy0 = (v2.x - v0.x) * inv_dy;
    float dxdy1 = (v2.x - v1.x) * inv_dy;

    float dzdy0 = (v2.z - v0.z) * inv_dy;
    float dzdy1 = (v2.z - v1.z) * inv_dy;

    float dwdy0 = (v2.w - v0.w) * inv_dy;
    float dwdy1 = (v2.w - v1.w) * inv_dy;

    int start_y = (int)(ceil(v0.y - 0.5f));
    int end_y = (int)(ceil(v2.y - 0.5f));

    V3 dcdy0 = v3_mul_f(v3_sub_v3(c2, c0), inv_dy);
    V3 dcdy1 = v3_mul_f(v3_sub_v3(c2, c1), inv_dy);

    V2 duvdy0 = v2_mul_f(v2_sub_v2(uv2, uv0), inv_dy);
    V2 duvdy1 = v2_mul_f(v2_sub_v2(uv2, uv1), inv_dy);

    for (int y = start_y; y < end_y; ++y) {
        // Must lerp for the vertex attributes to get them accurately.
        // TODO: Would be nice to find a way to step not lerp.
        //		 - not sure why stepping wouldn't work....

        // TODO: Step and not lerp.

        float a = (y + 0.5f - v0.y);

        float x0 = v0.x + dxdy0 * a;
        float x1 = v1.x + dxdy1 * a;

        float z0 = v0.z + dzdy0 * a;
        float z1 = v1.z + dzdy1 * a;

        float start_w = v0.w + dwdy0 * a;
        float end_w = v1.w + dwdy1 * a;

        int start_x = (int)((x0 - 0.5f));
        int end_x = (int)((x1 - 0.5f));

        V3 start_colour = v3_add_v3(c0, v3_mul_f(dcdy0, a));
        V3 end_colour = v3_add_v3(c1, v3_mul_f(dcdy1, a));

        V2 start_uv = v2_add_v2(uv0, v2_mul_f(duvdy0, a));
        V2 end_uv = v2_add_v2(uv1, v2_mul_f(duvdy1, a));

        draw_textured_scanline(rt, start_x, end_x, y, z0, z1, start_w, end_w, start_colour, end_colour, start_uv, end_uv, texture);
    }
}

void draw_textured_triangle(RenderTarget* rt, float* vc0, float* vc1, float* vc2, const Texture* texture)
{
    // TODO: why is it called vc0 (vertex copmonents) not just vertex.

	// TODO: All this stuff needs changing, however, important to remember that each texel actually defines the albedo
	//		 at that texel, the colour we're passing in here is really the diffuse contribution. I think I should refactor my code to focus on 
	//		 textured surfaces really. Obviously we're going to be introducing static lighting and light maps too.
	// 
	// 
    
    // TODO: NOTE, THIS is almost identical to the draw_triangle......


    // Sort vertices in ascending order.
    if (vc0[1] > vc1[1])
    {
        SWAP(float*, vc0, vc1);
    }
    if (vc0[1] > vc2[1])
    {
        SWAP(float*, vc0, vc2);
    }
    if (vc1[1] > vc2[1])
    {
        SWAP(float*, vc1, vc2);
    }

    // Handle if the triangle is already flat.
    if (vc0[1] == vc1[1])
    {
        draw_textured_flat_top_triangle(rt, vc0, vc1, vc2, texture);
        return;
    }

    if (vc1[1] == vc2[1])
    {
        draw_textured_flat_bottom_triangle(rt, vc0, vc1, vc2, texture);
        return;
    }

    // The triangle isn't flat, so split it into two flat triangles.
    float t = (vc1[1] - vc0[1]) / (vc2[1] - vc0[1]);

    // TODO: THis 9 should be defined as like TEXTURED_VERTEX_COMPONENTS, x,y,z,w,r,g,b,u,v
    static float vc3[9] = { 0 };

    // Lerp for v3. 
    vc3[0] = vc0[0] + (vc2[0] - vc0[0]) * t;
    vc3[1] = vc1[1];
    vc3[2] = vc0[2] + (vc2[2] - vc0[2]) * t;
    vc3[3] = vc0[3] + (vc2[3] - vc0[3]) * t;


    // TODO: This should be defined somewhere!!!!!!!!!

    // Lerp the components / attributes?.
    const int STRIDE_COMPONENTS = STRIDE_COLOUR + STRIDE_UV;

    // Lerp for each component.
    for (int i = 4; i < 4 + STRIDE_COMPONENTS; ++i)
    {
        vc3[i] = vc0[i] + (vc2[i] - vc0[i]) * t;
    }

    draw_textured_flat_top_triangle(rt, vc1, vc3, vc2, texture);
    draw_textured_flat_bottom_triangle(rt, vc0, vc1, vc3, texture);
}

void draw_depth_scanline(DepthBuffer* db, int x0, int x1, int y, float z0, float z1)
{
	// TODO: TEMP: Whilst no clipping.
	if (y < 0) return;
	if (y >= db->height) return;
	if (x0 == x1) return;
	if (x0 >= db->width) return;
	if (x1 < 0) return;

	if (x0 < 0) x0 = 0;
	if (x1 >= db->width) x1 = db->width - 1;



	// Precalculate deltas.
	const unsigned int dx = x1 - x0;
	float inv_dx = 1.f / dx;

	// Offset x by the given y.
	int row_offset = db->width * y;

	int start_x = x0 + row_offset;
	int end_x = x1 + row_offset;

	// Render the scanline
	float* depth_buffer = db->data + start_x;

	float z = z0;
	float z_step = (z1 - z0) * inv_dx;

	for (unsigned int i = 0; i < dx; ++i)
	{
		// Depth test, only draw closer values.
		//*depth_buffer = min(*depth_buffer, z);

		
		// TODO: profile
		if (*depth_buffer > z)
		{
			// TODO: min function for better speed?
			*depth_buffer = z;
		}

		// Move to the next pixel
		++depth_buffer;

		// Step per pixel values.
		z += z_step;
	}
}

void draw_depth_flat_bottom_triangle(DepthBuffer* db, V4 v0, V4 v1, V4 v2)
{
	// Sort the flat vertices left to right.
	if (v1.x > v2.x)
	{
        SWAP(V4, v1, v2);
	}

	float inv_dy = 1 / (v2.y - v0.y);

	float dxdy0 = (v1.x - v0.x) * inv_dy;
	float dxdy1 = (v2.x - v0.x) * inv_dy;

	float dzdy0 = (v1.z - v0.z) * inv_dy;
	float dzdy1 = (v2.z - v0.z) * inv_dy;

	float dwdy0 = (v1.w - v0.w) * inv_dy;
	float dwdy1 = (v2.w - v0.w) * inv_dy;

	int yStart = (int)(ceil(v0.y - 0.5f));
	int yEnd = (int)(ceil(v2.y - 0.5f));

	for (int y = yStart; y < yEnd; ++y) {
		// Must lerp for the vertex attributes otherwise the accuracy is poor.
		// TODO: Would be nice to not have to actually lerp but step instead.
		float a = (y + 0.5f - v0.y);

		// Calculate the start and ends of the scanline
		float x0 = v0.x + dxdy0 * a;
		float x1 = v0.x + dxdy1 * a;

		float z0 = v0.z + dzdy0 * a;
		float z1 = v0.z + dzdy1 * a;

		float wStart = v0.w + dwdy0 * a;
		float wEnd = v0.w + dwdy1 * a;


		int xStart = (int)(ceilf(x0 - 0.5f));
		int xEnd = (int)(ceilf(x1 - 0.5f));

		draw_depth_scanline(db, xStart, xEnd, y, z0, z1);
	}
}

void draw_depth_flat_top_triangle(DepthBuffer* db, V4 v0, V4 v1, V4 v2)
{
	// Sort the flat vertices left to right.
	if (v0.x > v1.x)
	{
        SWAP(V4, v0, v1);
	}

	float inv_dy = 1 / (v2.y - v0.y);

	float dxdy0 = (v2.x - v0.x) * inv_dy;
	float dxdy1 = (v2.x - v1.x) * inv_dy;

	float dzdy0 = (v2.z - v0.z) * inv_dy;
	float dzdy1 = (v2.z - v1.z) * inv_dy;

	float dwdy0 = (v2.w - v0.w) * inv_dy;
	float dwdy1 = (v2.w - v1.w) * inv_dy;

	int yStart = (int)(ceil(v0.y - 0.5f));
	int yEnd = (int)(ceil(v2.y - 0.5f));

	for (int y = yStart; y < yEnd; ++y) {
		// Must lerp for the vertex attributes to get them accurately.
		// TODO: Would be nice to find a way to step not lerp.
		float a = ((float)y + 0.5f - v0.y);

		float x0 = v0.x + dxdy0 * a;
		float x1 = v1.x + dxdy1 * a;

		float z0 = v0.z + dzdy0 * a;
		float z1 = v1.z + dzdy1 * a;

		float wStart = v0.w + dwdy0 * a;
		float wEnd = v1.w + dwdy1 * a;

		int xStart = (int)(ceil(x0 - 0.5f));
		int xEnd = (int)(ceil(x1 - 0.5f));
		draw_depth_scanline(db, xStart, xEnd, y, z0, z1);
	}
}

void draw_depth_triangle(DepthBuffer* db, V4 v0, V4 v1, V4 v2)
{
	// TODO: I don't think we need w, so should make these V3.

	// Sort vertices in ascending order.
	if (v0.y > v1.y)
	{
        SWAP(V4, v0, v1);
	}
	if (v0.y > v2.y)
	{
        SWAP(V4, v0, v2);
	}
	if (v1.y > v2.y)
	{
        SWAP(V4, v1, v2);
	}

	// Handle if the triangle is already flat.
	if (v0.y == v1.y)
	{
		draw_depth_flat_top_triangle(db, v0, v1, v2);
		return;
	}

	if (v1.y == v2.y)
	{
		draw_depth_flat_bottom_triangle(db, v0, v1, v2);
		return;
	}

	// The triangle isn't flat, so split it into two flat triangles.

	// Linear interpolate for v3.
	float t = (v1.y - v0.y) / (v2.y - v0.y);

	V4 v3 = {
		v0.x + (v2.x - v0.x) * t,
		v1.y,
		v0.z + (v2.z - v0.z) * t,
		v0.w + (v2.w - v0.w) * t
	};

	draw_depth_flat_top_triangle(db, v1, v3, v2);
	draw_depth_flat_bottom_triangle(db, v0, v1, v3);
}

float calculate_diffuse_factor(const V3 v, const V3 n, const V3 light_pos, float a, float b)
{
	// TODO: Comments, check maths etc.

	// calculate the direction of the light to the vertex
	V3 light_dir = v3_sub_v3(light_pos, v);


	float light_distance = size(light_dir);

	v3_mul_eq_f(&light_dir, 1.f / light_distance);
	
    // Calculate how much the vertex is lit
	float diffuse_factor = max(0.0f, dot(light_dir, n));

	// TODO: Just hardcode the attenuation factors here? Not sure we will need to change them.

	float attenuation = 1.0f / (1.0f + (a * light_distance) + (b * light_distance * light_distance));
	float dp = diffuse_factor * attenuation;

	// TODO: What is the name for this after attentuation is applied to the 
	// diffsue factor?
	return dp;
}

void project(const Canvas* canvas, const M4 projection_matrix, V4 v, V4* out)
{
	// TODO: Rename view space to screen space?

	// Opengl uses a right handed coordinate system, camera looks down the -z axis,
	// however, NDC space is left handed, from -1 to 1 in all axis. 
	// Therefore, the perspective projection matrix copies and inverts the 
	// initial depth z, to w' in v_projected.

	// Apply the perspective projection matrix to project
	// the 3D coordinates into 2D.
	V4 v_projected;
	m4_mul_v4(projection_matrix, v, &v_projected);

	// Perform perspective divide to bring to NDC space.
	// NDC space is a left handed coordinate system from -1 to 1 in all axis.
	const float inv_w = 1.0f / v_projected.w; // Precalculate the perspective divide.

	v_projected.x *= inv_w;
	v_projected.y *= inv_w;
	v_projected.z *= inv_w;

	// Convert from NDC space to screen space.
	// Convert from [-1:1] to [0:1], then scale to the screen dimensions.

	out->x = (v_projected.x + 1) * 0.5f * canvas->width;
	out->y = (-v_projected.y + 1) * 0.5f * canvas->height;

	// Projecting depth z results in z' which encodes a nonlinear transformation
	// of the depth, just like with x' and y'. So use this to depth test for 
	// more accurate results closer to the camera. Also, this means we only need
	// to recover the depth from w' if the depth test passes, saving a divison
	// per pixel.
	out->z = (v_projected.z + 1) * 0.5f; // Offset from [-1:1] to [0:1]

	// Save inv of w' for perspective correct interpolation. This allows us to lerp
	// between vertex components. 
	out->w = inv_w;
}

void model_to_view_space(ECS* ecs, System* render_system, FrameData* frame_data, Scene* scene, const M4 view_matrix)
{
    // TODO: Should this function really be defined as transform stage as 
    //		 we also do the bounding sphere stuff.....
    const MeshBase* mbs = scene->mesh_bases.bases;

    int vsps_offset = 0;
    int vsns_offset = 0;

    for (int si = 0; si < render_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = render_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        int mis_i = Archetype_find_component_list(archetype, COMPONENT_MeshInstance);
        MeshInstance* mis = archetype->component_lists[mis_i];

        int transforms_i = Archetype_find_component_list(archetype, COMPONENT_Transform);
        Transform* transforms = archetype->component_lists[transforms_i];

        V3* vsps = frame_data->view_space_positions.data;
        BoundingSphere* view_space_bounding_spheres = frame_data->view_space_bounding_spheres.data;

        for (int i = 0; i < archetype->entity_count; ++i)
        {            
            MeshInstance* mi = &mis[i];
            Transform transform = transforms[i];

            // Save the offset to the start of the view space positions for the 
            // mesh instance.
            mi->view_space_positions_offset = vsps_offset;

            // Calculate model matrix.
            // TODO: rotation or direction or eulers?
            M4 model_matrix;
            m4_model_matrix(transform.position, transform.rotation, transform.scale, model_matrix);

            M4 model_view_matrix;
            m4_mul_m4(view_matrix, model_matrix, model_view_matrix);

            const MeshBase* mb = &mbs[mi->mb_id];

            // Update the centre of the bounding sphere.
            V4 view_space_centre;
            m4_mul_v4(
                model_view_matrix,
                v3_to_v4(mb->centre, 1.f),
                &view_space_centre);

            // Write out the bounding sphere.
            BoundingSphere* bs = &view_space_bounding_spheres[i];
            bs->centre = v4_xyz(view_space_centre);

            // Read the mesh base's object space positions and convert them to 
            // view space for the mesh instance.

            for (int j = 0; j < mb->num_positions; ++j)
            {
                const V4 osp = v3_to_v4(mb->object_space_positions.data[j], 1.f);
                V4 vsp;
                m4_mul_v4(model_view_matrix, osp, &vsp);

                vsps[vsps_offset].x = vsp.x;
                vsps[vsps_offset].y = vsp.y;
                vsps[vsps_offset].z = vsp.z;

                ++vsps_offset;
            }
        }

        // Convert object space normals to view space.
        V3* vsns = frame_data->view_space_normals.data;
        
        for (int i = 0; i < archetype->entity_count; ++i)
        {
            MeshInstance* mi = &mis[i];
            Transform transform = transforms[i];

            // Save the offset to the start of the view space positions for the 
            // mesh instance.
            mi->view_space_normals_offset = vsns_offset;

            // Calculate normal matrix.
            // TODO: rotation or direction or eulers?
            M4 normal_matrix;
            m4_normal_matrix(transform.rotation, transform.scale, normal_matrix);

            M4 view_normal_matrix;
            m4_mul_m4(view_matrix, normal_matrix, view_normal_matrix);

            const MeshBase* mb = &mbs[mi->mb_id];

            for (int j = 0; j < mb->num_normals; ++j)
            {
                const V4 osn = v3_to_v4(mb->object_space_normals.data[j], 0.f);

                V4 vsn;
                m4_mul_v4(view_normal_matrix, osn, &vsn);

                V3 normal = normalised(v4_xyz(vsn));

                vsns[vsns_offset].x = normal.x;
                vsns[vsns_offset].y = normal.y;
                vsns[vsns_offset].z = normal.z;

                ++vsns_offset;
            }
        }

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            MeshInstance* mi = &mis[i];

            // Update the radius of each mesh instance bounding sphere.
            const V3* vsps = frame_data->view_space_positions.data;
            BoundingSphere* view_space_bounding_spheres = frame_data->view_space_bounding_spheres.data;

            // Only update the bounding sphere if the scale is changed, otherwise
            // we don't need to update it.
            if (mi->has_scale_changed)
            {
                mi->has_scale_changed = 0;

                const MeshBase* mb = &mbs[mi->mb_id];

                // Read bounding sphere (x,y,z,radius).
                BoundingSphere* bs = &view_space_bounding_spheres[i];
                V3 view_space_centre = bs->centre;

                // Calculate the new radius of the mi's bounding sphere.
                float radius_squared = -1;

                const int end = mi->view_space_positions_offset + mb->num_positions;
                for (int j = mi->view_space_positions_offset; j < end; ++j)
                {
                    V3 v = vsps[j];
                    V3 between = v3_sub_v3(v, view_space_centre);

                    radius_squared = max(size_squared(between), radius_squared);
                }

                // Save the radius (4th component of bounding sphere).
                bs->radius = sqrtf(radius_squared);
            }
        }
    }	
}

void lights_world_to_view_space(ECS* ecs, System* lighting_system, FrameData* frame_data, const Scene* scene, const M4 view_matrix)
{
	// This could be made more efficient by having an array of input world
	// space positions, however, we will never be able to support enough lights
	// to make this worthwhile (i think).


	V3* view_space_positions = frame_data->point_lights_view_space_positions.data;

    int vsps_offset = 0;

    for (int si = 0; si < lighting_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = lighting_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        int pls_i = Archetype_find_component_list(archetype, COMPONENT_PointLight);
        PointLight* pls = archetype->component_lists[pls_i];

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            V4 v_view_space;
            m4_mul_v4(view_matrix, v3_to_v4(pls[i].position, 1.f),
                &v_view_space);

            // There is no need to save the w component as it is always 1 until 
            // after projection.
            view_space_positions[vsps_offset].x = v_view_space.x;
            view_space_positions[vsps_offset].y = v_view_space.y;
            view_space_positions[vsps_offset].z = v_view_space.z;

            ++vsps_offset;
        }
    }
}

void broad_phase_frustum_culling(ECS* ecs, System* render_system, FrameData* frame_data, Scene* scene, const ViewFrustum* view_frustum)
{
	// Performs broad phase frustum culling on the models, writes out the planes
	// that need to be clipped against.
	
    MeshInstance* visible_mis = frame_data->visible_mis.data;
	int visible_mis_count = 0;
	uint8_t* intersected_planes = frame_data->intersected_planes.data;
	int intersected_planes_out_index = 0;

	const BoundingSphere* bounding_spheres = frame_data->view_space_bounding_spheres.data;

	const int num_planes = view_frustum->planes_count;
	const Plane* planes = view_frustum->planes;

    BoundingSphere* view_space_bounding_spheres = frame_data->view_space_bounding_spheres.data;

    for (int si = 0; si < render_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = render_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        int mis_i = Archetype_find_component_list(archetype, COMPONENT_MeshInstance);
        MeshInstance* mis = archetype->component_lists[mis_i];
    
        for (int i = 0; i < archetype->entity_count; ++i)
        {
            MeshInstance* mi = &mis[i];

            const BoundingSphere bs = bounding_spheres[i];

            // Store what planes need clipping against.
            int clip_against_plane[MAX_FRUSTUM_PLANES] = { 0 };
            int num_intersected_planes = 0;

            // Broad phase bounding sphere test.
            for (int j = 0; j < num_planes; ++j)
            {
                const float dist = signed_distance(&planes[j], bs.centre);
                if (dist < -bs.radius)
                {
                    // Completely outside the plane, therefore, no need to check against the others.
                    num_intersected_planes = -1; // -1 here means the mi is not visible.
                    break;
                }
                else if (dist < bs.radius)
                {
                    // Mark that we need to clip against this plane.
                    clip_against_plane[num_intersected_planes] = j;
                    ++num_intersected_planes;
                }
            }

            // Mark whether the mi passed the broad phase and store the 
            // intersection data if it passed.
            if (-1 != num_intersected_planes)
            {
                // Copy the mi to the visible mis array. TODO: Is this faster than using pointers??
                visible_mis[visible_mis_count++] = *mi;

                // Write out for narrow phase to use.
                // In format: num_planes_intersecting, plane_index_0, plane_index_1, ...
                intersected_planes[intersected_planes_out_index++] = num_intersected_planes;



                for (int j = 0; j < num_intersected_planes; ++j)
                {
                    intersected_planes[intersected_planes_out_index++] = clip_against_plane[j];
                }
            }
        }
    }

	frame_data->num_visible_mis = visible_mis_count;
}

void cull_backfaces(ECS* ecs, System* render_system, FrameData* frame_data, Scene* scene)
{
	const MeshBase* mbs = scene->mesh_bases.bases;

	// Determines what faces are facing the camera and prepares
	// the vertex data for the lighting calculations.
	const V3* vsps = frame_data->view_space_positions.data;
	MeshInstance* visible_mis = frame_data->visible_mis.data;
	const int num_visible_mis = frame_data->num_visible_mis;

	int* front_face_indices = frame_data->front_face_indices.data;
	int front_face_out = 0;


    

	for (int i = 0; i < num_visible_mis; ++i)
	{
		int front_faces = 0;

		MeshInstance* mi = &visible_mis[i];
		const MeshBase* mb = &mbs[mi->mb_id];

		const int* position_indices = mb->position_indices.data;
		const int vsp_offset = mi->view_space_positions_offset;

		for (int j = 0; j < mb->num_faces; ++j)
		{
			const int face_index = j * STRIDE_FACE_VERTICES;

			const int p0_index = vsp_offset + position_indices[face_index];
			const int p1_index = vsp_offset + position_indices[face_index + 1];
			const int p2_index = vsp_offset + position_indices[face_index + 2];
			
			// TODO: Would just storing as a V3 be better?
			if (is_front_face(vsps[p0_index], vsps[p1_index], vsps[p2_index]))
			{
				// Store the index of the face for culling later.
				front_face_indices[front_face_out++] = j;

				// TODO: Is it better to prep the data for the light stage here?
				//	     I don't really see why it would be tbf.

				++front_faces;
			}
		}

		// Update the number of front faces visible in the mi.
		// TODO: Gotta figure out whether it makes more sense to keep this
		//		 in the mi or in a separate frame data buffer, not sure.
		//		 Having it in frame data would mean that the next values
		//		 get prefetched which is dumb.
		mi->num_front_faces = front_faces;
	}
}

void light_front_faces(
    ECS* ecs, 
    System* render_system, 
    System* lighting_system, 
    FrameData* frame_data, 
    Scene* scene, 
    const V3 ambient)
{
	/*
	
	For this we need position, normal and abledo.

	However, what about with my shadows.

	For all none shadow casting lights, we can just calculate the total 
	contribution to the vertex.

	But for a shadow casting light....

	we must calculate the contribution to the vertex by each light, then per
	pixel we determine if that pixel is in shadow by that light and we must 
	add the light from the shadow.

	tbf my old approach is pretty solid.

	so output:
	- total vertex colour from non-shadow casting lights, including ambient.
	- a contribution from each shadow casting light
		when rendering the pixel, if it's not in shadow, we will add the light's
		contribution to the total.

	the outputs should be written to a transient framedata buffer i think, and
	then accessed from the array? Or how do we write them in a nicer way, not
	sure.

	but for the next step we need every component. i think though that we 
	just access them if visible, using a face/vertex index from the backface
	culling step.

	for this function do we want to iterate through faces that have passed
	the broad phase? Obviously we would have indirection when processing the 
	faces, but we don't have the branch, i believe the indirection is better.

	
	i think we have to calculate light contribution per vertex for the non-
	shadowing lights and the shadowing lights in the same step? Or not? Not
	sure what is better.
	*/

	// Input
	//const Lights* lights = &scene->lights;
    //const ComponentList* point_lights_list = &scene->lights.point_lights;
    //const PointLight* point_lights = (PointLight*)point_lights_list->elements;
	const MeshBase* mbs = scene->mesh_bases.bases;

	const int* front_face_indices = frame_data->front_face_indices.data;
	const V3* vsps = frame_data->view_space_positions.data;
	const V3* vsns = frame_data->view_space_normals.data;
	const V3* point_light_vsps = frame_data->point_lights_view_space_positions.data;

	const MeshInstance* mis = frame_data->visible_mis.data;
	const int num_visible_mis = frame_data->num_visible_mis;
	
	// Output
	int out_index = 0;

    // TODO: float* or V3*
	V3* vertex_lighting = frame_data->vertex_lighting.data;

	int front_faces_offset = 0;

	for (int i = 0; i < num_visible_mis; ++i)
	{
		const MeshInstance* mi = &mis[i];
		const MeshBase* mb = &mbs[mi->mb_id];

		const V3* vertex_albedos = mi->vertex_alebdos.data;

		const int vsp_offset = mi->view_space_positions_offset;
		const int vsn_offset = mi->view_space_normals_offset;

		for (int j = 0; j < mi->num_front_faces; ++j)
		{
			const int face_index = front_face_indices[j + front_faces_offset] * STRIDE_FACE_VERTICES;

			for (int vi = 0; vi < STRIDE_FACE_VERTICES; ++vi)
			{
				const V3 point = vsps[vsp_offset + mb->position_indices.data[face_index + vi]];
				const V3 normal = vsns[vsn_offset + mb->normal_indices.data[face_index + vi]];

                // Albedos are stored per vertex. 
                // TODO: How do these line up when we cull faces????? They dont..... do they? could cause issues when setting albedo for specific vertices....

                const V3 albedo = vertex_albedos[j * STRIDE_FACE_VERTICES + vi];
				V3 diffuse = { 0 };
                
                int pl_vsps_offset = 0;
                for (int si = 0; si < lighting_system->num_archetypes; ++si)
                {
                    const ArchetypeID archetype_id = lighting_system->archetype_ids[si];
                    Archetype* archetype = &ecs->archetypes[archetype_id];

                    int pls_i = Archetype_find_component_list(archetype, COMPONENT_PointLight);
                    PointLight* pls = archetype->component_lists[pls_i];

                    for (int i = 0; i < archetype->entity_count; ++i)
                    {
                        const V3 pl_vsp = point_light_vsps[pl_vsps_offset];
                        ++pl_vsps_offset;

                        const PointLight* pl = &pls[i];

                        const float a = 0.1f / pl->strength;
                        const float b = 0.01f / pl->strength;

                        const float df = calculate_diffuse_factor(point, normal, pl_vsp, a, b);

                        const V3 colour = v3_mul_f(pl->colour, df);

                        v3_add_eq_v3(&diffuse, colour);
                    }
                }


                /*
				for (int pl_index = 0; pl_index < point_lights_list->count; ++pl_index)
				{
					const PointLight* pl = &point_lights[pl_index];

					const V3 pl_vsp = v3_read(point_light_vsps + pl_index * STRIDE_POSITION);

					const float a = 0.1f / pl->strength;
					const float b = 0.01f / pl->strength;

					const float df = calculate_diffuse_factor(point, normal, pl_vsp, a, b);
                    
					const V3 colour = v3_mul_f(pl->colour, df);

					v3_add_eq_v3(&diffuse, colour);
				}*/

				// TODO: shadow casting lights etc.

				// Clamp diffuse contribution to a valid range 0-1. 
                V3 light = v3_mul_v3(albedo, v3_add_v3(diffuse, ambient));


				// TODO: Do we clamp here? How does it work with shadows..
                
                //i think the idea with the shadows this time will be to add
                //their colour contribution if not in shadow.
                
                
				// Clamp to max, should never be negative.
				light.x = min(1.f, light.x);
				light.y = min(1.f, light.y);
				light.z = min(1.f, light.z);


				vertex_lighting[out_index].x = light.x;
				vertex_lighting[out_index].y = light.y;
				vertex_lighting[out_index++].z = light.z;

                // TODO: TEMP
                //break;
			}
            // TODO: TEMP
            //break;

		}

		front_faces_offset += mi->num_front_faces;
	}
}

void prepare_for_clipping(ECS* ecs, System* render_system, FrameData* frame_data, Scene* scene)
{
    // TODO: Comments on how the data is packed together.

	// Input
	const MeshInstance* mis = frame_data->visible_mis.data;
	const MeshBase* mbs = scene->mesh_bases.bases;
	
	const int* front_face_indices = frame_data->front_face_indices.data;

	const int num_visible_mis = frame_data->num_visible_mis;

	const V3* vertex_lighting = frame_data->vertex_lighting.data;
    int vertex_lighting_in = 0;

	// Output
	float* faces_to_clip = frame_data->faces_to_clip.data;
	int out_index = 0;

	int front_faces_offset = 0;

	for (int i = 0; i < num_visible_mis; ++i)
	{
		const MeshInstance* mi = &mis[i];
		const MeshBase* mb = &mbs[mi->mb_id];

		const int* position_indices = mb->position_indices.data;
		const int vsp_offset = mi->view_space_positions_offset;
	    const V3* vsps = frame_data->view_space_positions.data + vsp_offset;

		const int* uv_indices = mb->uv_indices.data;
		const V2* uvs = mb->uvs.data;

        // Only write out UVs if the mesh is textured.
        if (mi->texture_id != -1)
        {
            for (int j = 0; j < mi->num_front_faces; ++j)
            {
                const int face_index = front_face_indices[j + front_faces_offset] * STRIDE_FACE_VERTICES;

                const V3 p0 = vsps[position_indices[face_index]];
                const V3 p1 = vsps[position_indices[face_index + 1]];
                const V3 p2 = vsps[position_indices[face_index + 2]];

                // TODO: Lighting is not the nicest name for this.
                const V3 lighting0 = vertex_lighting[vertex_lighting_in];
                const V3 lighting1 = vertex_lighting[vertex_lighting_in + 1];
                const V3 lighting2 = vertex_lighting[vertex_lighting_in + 2];
                vertex_lighting_in += STRIDE_FACE_VERTICES;

                const V2 uv0 = uvs[uv_indices[face_index]];
                const V2 uv1 = uvs[uv_indices[face_index + 1]];
                const V2 uv2 = uvs[uv_indices[face_index + 2]];

                // Write out the vertices.
#define WRITE_VERTEX(p, l, uv) \
                do { \
                    faces_to_clip[out_index++] = (p).x; \
                    faces_to_clip[out_index++] = (p).y; \
                    faces_to_clip[out_index++] = (p).z; \
                    faces_to_clip[out_index++] = (l).x; \
                    faces_to_clip[out_index++] = (l).y; \
                    faces_to_clip[out_index++] = (l).z; \
                    faces_to_clip[out_index++] = (uv).x; \
                    faces_to_clip[out_index++] = (uv).y; \
                } while (0)

                WRITE_VERTEX(p0, lighting0, uv0);
                WRITE_VERTEX(p1, lighting1, uv1);
                WRITE_VERTEX(p2, lighting2, uv2);
#undef WRITE_VERTEX
            }
        }
        else
        {
            for (int j = 0; j < mi->num_front_faces; ++j)
            {
                const int face_index = front_face_indices[j + front_faces_offset] * STRIDE_FACE_VERTICES;

                const V3 p0 = vsps[position_indices[face_index]];
                const V3 p1 = vsps[position_indices[face_index + 1]];
                const V3 p2 = vsps[position_indices[face_index + 2]];

                // TODO: Lighting is not the nicest name for this.
                const V3 lighting0 = vertex_lighting[vertex_lighting_in];
                const V3 lighting1 = vertex_lighting[vertex_lighting_in + 1];
                const V3 lighting2 = vertex_lighting[vertex_lighting_in + 2];
                vertex_lighting_in += STRIDE_FACE_VERTICES;

                // Write out the vertices.
#define WRITE_VERTEX(p, l) \
                do { \
                    faces_to_clip[out_index++] = (p).x; \
                    faces_to_clip[out_index++] = (p).y; \
                    faces_to_clip[out_index++] = (p).z; \
                    faces_to_clip[out_index++] = (l).x; \
                    faces_to_clip[out_index++] = (l).y; \
                    faces_to_clip[out_index++] = (l).z; \
                } while (0)

                WRITE_VERTEX(p0, lighting0);
                WRITE_VERTEX(p1, lighting1);
                WRITE_VERTEX(p2, lighting2);
#undef WRITE_VERTEX
            }
        }

        // Processed the entire mi so offset to the next one.
        front_faces_offset += mi->num_front_faces;
    }	
}

void clip_project_and_draw(
    System* render_system,
	Renderer* renderer,
	RenderTarget* rt,
	FrameData* frame_data,
	Scene* scene,
    const Resources* resources)
{
	// Input
	const MeshInstance* mis = frame_data->visible_mis.data;
	const MeshBase* mbs = scene->mesh_bases.bases;

	const int num_visible_mis = frame_data->num_visible_mis;

	const uint8_t* intersected_planes = frame_data->intersected_planes.data;
	int intersected_planes_index = 0;

	// Output
	float* faces_to_clip = frame_data->faces_to_clip.data;
	float* clipped_faces = frame_data->clipped_faces.data;
	int out_index = 0;

	int face_offset = 0;

	for (int i = 0; i < num_visible_mis; ++i)
	{
		const MeshInstance* mi = &mis[i];
		const MeshBase* mb = &mbs[mi->mb_id];

		const uint8_t num_planes_to_clip_against = intersected_planes[intersected_planes_index++];
        
        // TODO: Should be defined somewhere.
        // If textured, include UVs.
        const int VERTEX_COMPONENTS = mi->texture_id == -1 ? 6 : 8; // x,y,z,r,g,b + u,v
        const int FACE_COMPONENTS = VERTEX_COMPONENTS * STRIDE_FACE_VERTICES;
        
        // TODO: At this point we know if we need UVs or not.

		// Entire mesh is visible so just copy the vertices over.
		if (num_planes_to_clip_against == 0)
		{			
			memcpy(clipped_faces, 
                   faces_to_clip + face_offset, 
                   (size_t)mi->num_front_faces * FACE_COMPONENTS * sizeof(float));

            // Write out the number of faces to draw.
			frame_data->num_clipped_faces = mi->num_front_faces;
		}
		else
        {
            // Calculate the number of components to lerp per vertex.
            const int COMPS_TO_LERP = VERTEX_COMPONENTS - STRIDE_POSITION;

			// Partially inside so must clip the vertices against the planes.

			// Initially read from the front_faces buffer.
			float* temp_clipped_faces_in = frame_data->faces_to_clip.data;
			float* temp_clipped_faces_out = frame_data->temp_clipped_faces1.data;

			// Store the index to write out to, needs to be defined here so we can
			// update the clipped_faces_index after writing to the clipped_faces buffer.
			int index_out = 0;

			// After each plane, we will have a different number of faces to clip again.
			// Initially set this to the number of front faces.
			int num_faces_to_process = mi->num_front_faces;

			// This is needed as an offset into the front_faces buffer for the first plane.
			int clipped_faces_offset = face_offset;

			// The logic for setting in/out buffers depends on how many planes we actually
			// render against, not the plane index.
			int num_planes_clipped_against = 0;

			for (int j = 0; j < num_planes_to_clip_against; ++j)
			{
				const Plane* plane = &renderer->settings.view_frustum.planes[intersected_planes[intersected_planes_index++]];
				
				// Reset the index to write out to.
				index_out = 0;

				// Store how many triangles were wrote to the out buffer.				
				int temp_visible_faces_count = 0;

				// After the first plane we want to read from the in buffer.
                // Now we're not reading from the 
				if (num_planes_clipped_against == 1)
				{
					// Initially we used the in front faces buffer, after the first iteration 
					// we have wrote to the out buffer, so that can now be our in buffer.
					temp_clipped_faces_out = frame_data->temp_clipped_faces0.data;

					// Now we want to read from the start of the in buffer, 
					// not the offset into the front faces buffer.
					clipped_faces_offset = 0;
				}

				// If we're processing the last plane, write out to the clipped_faces buffer.
				if (num_planes_clipped_against == num_planes_to_clip_against - 1)
				{
					// On the last plane, we want to write out to the clipped faces.
					temp_clipped_faces_out = frame_data->clipped_faces.data;
				}

				// For faces in mesh.
				for (int j = 0; j < num_faces_to_process; ++j)
				{
					int index_face = clipped_faces_offset + j * FACE_COMPONENTS;

					int num_inside_points = 0;
					int num_outside_points = 0;

					int inside_points_indices[3] = { 0 };
					int outside_points_indices[3] = { 0 };

                    // TODO: Do we want to be dealing with indices here? Could just have the positions from the 
                    //       step before right?
					const int index_v0 = index_face;
					const int index_v1 = index_face + VERTEX_COMPONENTS;
					const int index_v2 = index_face + VERTEX_COMPONENTS + VERTEX_COMPONENTS;

					const float d0 = signed_distance(plane, v3_read(temp_clipped_faces_in + index_v0));
					const float d1 = signed_distance(plane, v3_read(temp_clipped_faces_in + index_v1));
					const float d2 = signed_distance(plane, v3_read(temp_clipped_faces_in + index_v2));

					// Determine what points are inside and outside the plane.
					if (d0 >= 0)
					{
						inside_points_indices[num_inside_points++] = index_v0;
					}
					else
					{
						outside_points_indices[num_outside_points++] = index_v0;
					}
					if (d1 >= 0)
					{
						inside_points_indices[num_inside_points++] = index_v1;
					}
					else
					{
						outside_points_indices[num_outside_points++] = index_v1;
					}
					if (d2 >= 0)
					{
						inside_points_indices[num_inside_points++] = index_v2;
					}
					else
					{
						outside_points_indices[num_outside_points++] = index_v2;
					}

					if (num_inside_points == 3)
					{
						// The whole triangle is inside the plane, so copy the face.
						memcpy(temp_clipped_faces_out + index_out, temp_clipped_faces_in + index_face, FACE_COMPONENTS * sizeof(float));

						index_out += FACE_COMPONENTS;
						++temp_visible_faces_count;
					}
					else if (num_inside_points == 1 && num_outside_points == 2)
					{
						// Form a new triangle with the plane edge.
						const int index_iv0 = inside_points_indices[0];
						const int index_ov0 = outside_points_indices[0];
						const int index_ov1 = outside_points_indices[1];

                        const float* iv0 = temp_clipped_faces_in + index_iv0;
                        const float* ov0 = temp_clipped_faces_in + index_ov0;
                        const float* ov1 = temp_clipped_faces_in + index_ov1;

						const V3 ip0 = v3_read(iv0);
						const V3 op0 = v3_read(ov0);
						const V3 op1 = v3_read(ov1);

						// Copy the inside vertex.
						memcpy(temp_clipped_faces_out + index_out, iv0, VERTEX_COMPONENTS * sizeof(float));
						index_out += VERTEX_COMPONENTS;

						// Lerp for the first new vertex.
						V3 p0;
						float t = line_intersect_plane(ip0, op0, plane, &p0);

						temp_clipped_faces_out[index_out++] = p0.x;
						temp_clipped_faces_out[index_out++] = p0.y;
						temp_clipped_faces_out[index_out++] = p0.z;

						// Lerp the vertex components straight into the out buffer.
                        const float* from = iv0 + STRIDE_POSITION;
                        const float* to = ov0 + STRIDE_POSITION;
						for (int k = 0; k < COMPS_TO_LERP; ++k)
						{
							temp_clipped_faces_out[index_out++] = lerp(from[k], to[k], t);
						}

						// Lerp for the second vertex.
						V3 p1;
						t = line_intersect_plane(ip0, op1, plane, &p1);

						temp_clipped_faces_out[index_out++] = p1.x;
						temp_clipped_faces_out[index_out++] = p1.y;
						temp_clipped_faces_out[index_out++] = p1.z;

                        // Lerp components.
                        to = ov1 + STRIDE_POSITION;
						for (int k = 0; k < COMPS_TO_LERP; ++k)
						{
							temp_clipped_faces_out[index_out++] = lerp(from[k], to[k], t);
						}
						
						++temp_visible_faces_count;
					}
					else if (num_inside_points == 2 && num_outside_points == 1)
					{
						// Form two new triangles with the plane edge.
						const int index_iv0 = inside_points_indices[0];
						const int index_iv1 = inside_points_indices[1];
						const int index_ov0 = outside_points_indices[0];

                        const float* iv0 = temp_clipped_faces_in + index_iv0;
                        const float* iv1 = temp_clipped_faces_in + index_iv1;
                        const float* ov0 = temp_clipped_faces_in + index_ov0;

						const V3 ip0 = v3_read(iv0);
						const V3 ip1 = v3_read(iv1);
						const V3 op0 = v3_read(ov0);

						// Copy the first inside vertex.
						memcpy(temp_clipped_faces_out + index_out, iv0, VERTEX_COMPONENTS * sizeof(float));
						index_out += VERTEX_COMPONENTS;

						// Copy the second inside vertex.
						memcpy(temp_clipped_faces_out + index_out, iv1, VERTEX_COMPONENTS * sizeof(float));
						index_out += VERTEX_COMPONENTS;

						// Lerp for the first new vertex.
						V3 p0;
						float t = line_intersect_plane(ip0, op0, plane, &p0);

                        // We will reuse this vertex for the second triangle, so copy it's index.
						int new_v0_index = index_out;

						temp_clipped_faces_out[index_out++] = p0.x;
						temp_clipped_faces_out[index_out++] = p0.y;
						temp_clipped_faces_out[index_out++] = p0.z;

                        const float* from = iv0 + STRIDE_POSITION;
                        const float* to = ov0 + STRIDE_POSITION;
						for (int k = 0; k < COMPS_TO_LERP; ++k)
						{
							temp_clipped_faces_out[index_out++] = lerp(from[k], to[k], t);
						}
						
						++temp_visible_faces_count;

						// First triangle done.
						
						// Copy the initial new vertex for the second triangle.
						for (int k = new_v0_index; k < new_v0_index + VERTEX_COMPONENTS; ++k)
						{
							temp_clipped_faces_out[index_out++] = temp_clipped_faces_out[k];
						}

						// Lerp for the second new point.
						V3 p1;
						t = line_intersect_plane(ip1, op0, plane, &p1);
						
						temp_clipped_faces_out[index_out++] = p1.x;
						temp_clipped_faces_out[index_out++] = p1.y;
						temp_clipped_faces_out[index_out++] = p1.z;

                        from = iv1 + STRIDE_POSITION;
						for (int k = 0; k < COMPS_TO_LERP; ++k)
						{
							temp_clipped_faces_out[index_out++] = lerp(from[k], to[k], t);
						}

						// Copy over the inside vertex for this triangle.
						memcpy(temp_clipped_faces_out + index_out, iv1, VERTEX_COMPONENTS * sizeof(float));
						index_out += VERTEX_COMPONENTS;
						
						++temp_visible_faces_count;
					}
				}

				// Update how many faces are visible after being clipped.
				num_faces_to_process = temp_visible_faces_count;

				// Swap the in and out buffers.
                SWAP(float*, temp_clipped_faces_in, temp_clipped_faces_out);

				// Increment to show we actually clipped a plane.
				++num_planes_clipped_against;
			}			

            frame_data->num_clipped_faces = num_faces_to_process;
		}

        // Draw the clipped face
        if (frame_data->num_clipped_faces > 0)
        {
            if (mi->texture_id == -1)
            {
                project_and_draw_clipped(renderer, rt, frame_data, mi);
            }
            else
            {
                project_and_draw_clipped_textured(renderer, rt, frame_data, mi,
                    &resources->textures[mi->texture_id]);
            }
        }

		face_offset += mi->num_front_faces * FACE_COMPONENTS;
	}
}

void project_and_draw_clipped(
	Renderer* renderer,
	RenderTarget* rt,
	FrameData* frame_data,
	const MeshInstance* mi
	)
{
    // TODO: I would at least like to project everything at once if we're going to do it like this right??
    //       should be much better for cache.

    // TODO: VERTEX_COMPONENTS should be defined somewhere as it's used in multiple functions.

	// TODO: Hardcoded for now.
	// TODO: WITH NO UVS FOR NOW
#define VERTEX_COMPONENTS 7 //x,y,z,w,r,g,b
	const int INPUT_VERTEX_COMPONENTS = 6;

    // Input
	const float* clipped_faces = frame_data->clipped_faces.data;
	const int num_clipped_faces = frame_data->num_clipped_faces;

	for (int i = 0; i < num_clipped_faces; ++i)
	{
		const int clipped_face_index = i * INPUT_VERTEX_COMPONENTS * STRIDE_FACE_VERTICES;

		const float* face = clipped_faces + clipped_face_index;

		const float* v0 = face;
		const float* v1 = face + INPUT_VERTEX_COMPONENTS;
		const float* v2 = face + INPUT_VERTEX_COMPONENTS * 2;

		const V4 p0 = v3_read_to_v4(v0, 1.f);
		const V4 p1 = v3_read_to_v4(v1, 1.f);
		const V4 p2 = v3_read_to_v4(v2, 1.f);

		V4 ssp0, ssp1, ssp2;
		project(&rt->canvas, renderer->settings.projection_matrix, p0, &ssp0);
		project(&rt->canvas, renderer->settings.projection_matrix, p1, &ssp1);
		project(&rt->canvas, renderer->settings.projection_matrix, p2, &ssp2);

		static float vc0[VERTEX_COMPONENTS] = { 0 };
		v4_write(vc0, ssp0);
		for (int j = 4; j < VERTEX_COMPONENTS; ++j)
		{
			// -1 because incoming x,y,z -> out x,y,z,w
			vc0[j] = v0[j-1] * ssp0.w;
		}

        static float vc1[VERTEX_COMPONENTS] = { 0 };
		v4_write(vc1, ssp1);
		for (int j = 4; j < VERTEX_COMPONENTS; ++j)
		{
			vc1[j] = v1[j-1] * ssp1.w;
		}

        static float vc2[VERTEX_COMPONENTS] = { 0 };
		v4_write(vc2, ssp2);
		for (int j = 4; j < VERTEX_COMPONENTS; ++j)
		{
			vc2[j] = v2[j-1] * ssp2.w;
		}

		draw_triangle(rt, vc0, vc1, vc2);
	}

#undef VERTEX_COMPONENTS
}

void project_and_draw_clipped_textured(
    Renderer* renderer,
    RenderTarget* rt,
    FrameData* frame_data,
    const MeshInstance* mi,
    const Texture* texture
)
{
    // TODO: Almost identical to project_and_draw_clipped..... maybe not worth
    //       the tiny performance gain.

    // TODO: I would at least like to project everything at once if we're going to do it like this right??
    //       should be much better for cache.

    // TODO: VERTEX_COMPONENTS should be defined somewhere as it's used in multiple functions.

    // TODO: Hardcoded for now.
#define VERTEX_COMPONENTS 9 //x,y,z,w,r,g,b
    const int INPUT_VERTEX_COMPONENTS = 8;

    // Input
    const float* clipped_faces = frame_data->clipped_faces.data;
    const int num_clipped_faces = frame_data->num_clipped_faces;

    for (int i = 0; i < num_clipped_faces; ++i)
    {
        const int clipped_face_index = i * INPUT_VERTEX_COMPONENTS * STRIDE_FACE_VERTICES;

        const float* face = clipped_faces + clipped_face_index;

        const float* v0 = face;
        const float* v1 = face + INPUT_VERTEX_COMPONENTS;
        const float* v2 = face + INPUT_VERTEX_COMPONENTS * 2;

        const V4 p0 = v3_read_to_v4(v0, 1.f);
        const V4 p1 = v3_read_to_v4(v1, 1.f);
        const V4 p2 = v3_read_to_v4(v2, 1.f);

        V4 ssp0, ssp1, ssp2;
        project(&rt->canvas, renderer->settings.projection_matrix, p0, &ssp0);
        project(&rt->canvas, renderer->settings.projection_matrix, p1, &ssp1);
        project(&rt->canvas, renderer->settings.projection_matrix, p2, &ssp2);

        static float vc0[VERTEX_COMPONENTS] = { 0 };
        v4_write(vc0, ssp0);
        for (int j = 4; j < VERTEX_COMPONENTS; ++j)
        {
            // -1 because incoming x,y,z -> out x,y,z,w
            vc0[j] = v0[j - 1] * ssp0.w;
        }

        static float vc1[VERTEX_COMPONENTS] = { 0 };
        v4_write(vc1, ssp1);
        for (int j = 4; j < VERTEX_COMPONENTS; ++j)
        {
            vc1[j] = v1[j - 1] * ssp1.w;
        }

        static float vc2[VERTEX_COMPONENTS] = { 0 };
        v4_write(vc2, ssp2);
        for (int j = 4; j < VERTEX_COMPONENTS; ++j)
        {
            vc2[j] = v2[j - 1] * ssp2.w;
        }

        // TODO: Would be nice to do this elsewhere
        // Scale the uvs to the size of the texture to avoid doing it 
        // per pixel. 
        // TODO: This could be done via a MeshInstance set texture
        // function?
        int w = texture->width - 1;
        int h = texture->height - 1;

        vc0[7] *= w;
        vc0[8] *= h;
        vc1[7] *= w;
        vc1[8] *= h;
        vc2[7] *= w;
        vc2[8] *= h;

        draw_textured_triangle(rt, vc0, vc1, vc2, texture);
    }
#undef VERTEX_COMPONENTS
}

void render(
    ECS* ecs,
    System* render_system,
    System* lighting_system,
	Renderer* renderer,
	Scene* scene,
	const Resources* resources,
	const M4 view_matrix)
{
    // This render function is using the render and lighting 'systems' honestly
    // they're more like views at this point.

    // TODO: I realy don't like passing the ecs and systems, feels very messy.
    //       Not sure though, doesn't matter loads. 
	
	frame_data_init(ecs, render_system, lighting_system, &renderer->frame_data, scene);

	// Convert positions and normals from object space to view space.
	// Also update mesh instance's bounding spheres.
	// TODO: Rename object space? or model space??/
	model_to_view_space(ecs, render_system, &renderer->frame_data, 
		scene, view_matrix);

	// Convert light positions from world space to view space.
	lights_world_to_view_space(ecs, lighting_system, &renderer->frame_data,
		scene, view_matrix);

	broad_phase_frustum_culling(ecs, render_system, &renderer->frame_data, scene, &renderer->settings.view_frustum);

	cull_backfaces(ecs, render_system, &renderer->frame_data, scene);

    // TODO: Need to profile actually to find the issues before doing all this stupid stuff.

    // TODO: Honestly a lot of this code just feels awful to read.

	light_front_faces(ecs, render_system, lighting_system, &renderer->frame_data, scene, scene->ambient_light);

	prepare_for_clipping(ecs, render_system, &renderer->frame_data, scene);

	clip_project_and_draw(render_system, renderer, &renderer->target, &renderer->frame_data, scene, resources);

    // DEBUGGING
    debug_draw_point_lights(
        &renderer->target.canvas, 
        ecs,
        lighting_system,
        &renderer->frame_data,
        &renderer->settings
    );

   // debug_draw_normals(&renderer->target.canvas, &renderer->frame_data, &renderer->settings, scene);
}
/*
void update_depth_maps(Renderer* renderer, const Scene* scene)
{

	// TODO: To avoid recalculating shadows for each mesh. We should only do it for ones that have moved.
	//		 This should save us a lot of computation I believe. But also means we need a static and dynamic
	//		 depth buffer. Something like this anyways too tired to think of now. 
	// 
	//		 Need to look into this properly.
	//		 But for example, for the main map, that's not moving, we might have a single directional (or other) light as the 
	//		 sun, we could generate a depth map from all the static geometry, we would have to draw dynamic to a separate. But
	//		 this would give us proper shadows and essentially just make it a texture lookup.

	// TODO: At some point we definitely want to be able to render a directional light for the sun/moon, 
	//		 then the environment can have a static shadow map. Then the dynamic stuff can be renderered to a separate map.
	
	// TODO: Rename shadow maps?

	// TODO: Instead of 6 different textures for a point light, have one big/long cubemap?

	const Models* models = &scene->models;
	const PointLights* pls = &scene->point_lights;

	for (int i = 0; i < pls->count; ++i)
	{
		DepthBuffer* depth_map = &pls->depth_maps[i];
		depth_buffer_fill(depth_map, 1.f);

		int pos_i = i * STRIDE_POSITION;

		V3 pos = v3_read(pls->world_space_positions + pos_i);
		
		// TODO: TEMP, hardcoded.
		V3 dir = { 0, 0, -1 };
		
		// Create MV matrix for light.
		M4 view;
		look_at(v3_mul_f(pos, -1.f), v3_mul_f(dir, -1.f), view);

		// 90 degrees will give us a face of the cube map we want.
		float fov = 90.f;
		float aspect_ratio = depth_map->width / (float)depth_map->height;

		// TODO: TEMP: Hardcoded settings

		float near_plane = 1.f;
		float far_plane = 100.f; // TODO: Defined from strength of point light? with attenuation taken into account?

		M4 proj;
		m4_projection(fov, aspect_ratio, near_plane, far_plane, proj);
		
		const int mis_count = models->mis_count;

		const float* mis_transforms = models->mis_transforms;

		const int* mis_base_ids = models->mis_base_ids;

		const int* mbs_positions_counts = models->mbs_positions_counts;
		const int* mbs_positions_offsets = models->mbs_positions_offsets;
		const float* object_space_positions = models->mbs_object_space_positions;

		
		// TODO: Rename vars.

		int positions_offset = 0;

		int light_out_index = 0;

		for (int j = 0; j < mis_count; ++j)
		{
			// Convert the model base object space positions to world space
			// for the current model instance.
			const int mb_index = mis_base_ids[j];
			
			// Calculate the new model matrix from the mi's transform.
			int transform_index = j * STRIDE_MI_TRANSFORM;
			
			M4 model_matrix;
			m4_model_matrix(
				v3_read(mis_transforms + transform_index), 
				v3_read(mis_transforms + transform_index + 3), 
				v3_read(mis_transforms + transform_index + 6), 
				model_matrix
			);

			M4 model_view;
			m4_mul_m4(view, model_matrix, model_view);

			for (int k = 0; k < models->mbs_faces_counts[mb_index]; ++k)
			{
				const int face_index = (models->mbs_faces_offsets[mb_index] + k) * STRIDE_FACE_VERTICES;

				// Get the indices to the first component of each vertex position.
				const int index_v0 = models->mbs_face_position_indices[face_index] + mbs_positions_offsets[mb_index];
				const int index_v1 = models->mbs_face_position_indices[face_index + 1] + mbs_positions_offsets[mb_index];
				const int index_v2 = models->mbs_face_position_indices[face_index + 2] + mbs_positions_offsets[mb_index];

				const int index_parts_v0 = index_v0 * STRIDE_POSITION;
				const int index_parts_v1 = index_v1 * STRIDE_POSITION;
				const int index_parts_v2 = index_v2 * STRIDE_POSITION;

				// Get the vertices from the face indices.
				V4 osp0 = v3_read_to_v4(object_space_positions + index_parts_v0, 1.f);
				V4 osp1 = v3_read_to_v4(object_space_positions + index_parts_v1, 1.f);
				V4 osp2 = v3_read_to_v4(object_space_positions + index_parts_v2, 1.f);

				// TODO: For this, refactor to do model view transformation first per vertex
				//		 so we're not doing it multiple times with the indexed rendering. 
				//		
				// TODO: This should all work the same as the normal rendering really, frustum
				//		 culling and clipping etc.

				V4 vsp0, vsp1, vsp2;
				m4_mul_v4(model_view, osp0, &vsp0);
				m4_mul_v4(model_view, osp1, &vsp1);
				m4_mul_v4(model_view, osp2, &vsp2);

				V3 vsp0_v3 = v4_xyz(vsp0);
				V3 vsp1_v3 = v4_xyz(vsp1);
				V3 vsp2_v3 = v4_xyz(vsp2);

				V3 face_normal = normalised(cross(v3_sub_v3(vsp1_v3, vsp0_v3), v3_sub_v3(vsp2_v3, vsp0_v3)));

				// Only fill depth map from back faces, need the normal so doing this manually.
				if (dot(vsp0_v3, face_normal) <= 0)
				{
					// TODO: What do we write out now???
					continue;
				}

				// The perspective projection transforms the coordinates into clip space, before the perpsective divide.
				// which just converts the homogeneous coordinates to cartesian ones. 
				V4 clip0, clip1, clip2;
				m4_mul_v4(proj, vsp0, &clip0);
				m4_mul_v4(proj, vsp1, &clip1);
				m4_mul_v4(proj, vsp2, &clip2);

				// Save the light space positions for each vertex.
				const int index_lsp_parts_v0 = (models->mbs_face_position_indices[face_index] + positions_offset) * 4;
				const int index_lsp_parts_v1 = (models->mbs_face_position_indices[face_index + 1] + positions_offset) * 4;
				const int index_lsp_parts_v2 = (models->mbs_face_position_indices[face_index + 2] + positions_offset) * 4;

				// Perspective-correct intertpolation of values that have undergone perspective divide dont work?
				// TODO: Comments and understand all this a bit more.
				v4_write(renderer->buffers.light_space_positions + index_lsp_parts_v0, clip0);
				v4_write(renderer->buffers.light_space_positions + index_lsp_parts_v1, clip1);
				v4_write(renderer->buffers.light_space_positions + index_lsp_parts_v2, clip2);

				// Perform perspective divide to convert Clip Space to NDC space (-1:1 for x,y,z).
				const float inv_w0 = 1.0f / clip0.w;
				const float inv_w1 = 1.0f / clip1.w;
				const float inv_w2 = 1.0f / clip2.w;
				
				V4 ndc0 = {
					clip0.x * inv_w0,
					clip0.y * inv_w0,
					clip0.z * inv_w0,
					inv_w0
				};

				V4 ndc1 = {
					clip1.x * inv_w1,
					clip1.y * inv_w1,
					clip1.z * inv_w1,
					inv_w1
				};

				V4 ndc2 = {
					clip2.x * inv_w2,
					clip2.y * inv_w2,
					clip2.z * inv_w2,
					inv_w2
				};

				// Convert NDC to screen space by first converting to 0-1 in all axis.
				V4 ssp0 = {
					(ndc0.x + 1) * 0.5f * depth_map->width,
					(-ndc0.y + 1) * 0.5f * depth_map->height,
					(ndc0.z + 1) * 0.5f,
					inv_w0
				};

				V4 ssp1 = {
					(ndc1.x + 1) * 0.5f * depth_map->width,
					(-ndc1.y + 1) * 0.5f * depth_map->height,
					(ndc1.z + 1) * 0.5f,
					inv_w1
				};

				V4 ssp2 = {
					(ndc2.x + 1) * 0.5f * depth_map->width,
					(-ndc2.y + 1) * 0.5f * depth_map->height,
					(ndc2.z + 1) * 0.5f,
					inv_w2
				};
				
				// Apply slope scaled depth bias to fix shadow acne and peter panning.
				
				// Don't allow the cos angle to be negative, the bias should push the shadow away from the light.
				float cos_theta = fabsf(dot(face_normal, v3_mul_f(dir, -1.f)));

				// TODO: This will have to be changed for each scene i think.
				const float constant_bias = 0.00001f;

				// Clamp the cos_theta to a value near 0 so we don't divide by 0.
				cos_theta = max(cos_theta, 0.00001f);
				
				// We want a large bias when the light dir and surface dir are perpendicular
				// because shadow acne is most common there.
				float slope_bias = constant_bias * sqrtf(1.f - cos_theta * cos_theta) / cos_theta;
				ssp0.z += slope_bias;
				ssp1.z += slope_bias;
				ssp2.z += slope_bias;


				// Draw to the depth buffer.
				draw_depth_triangle(depth_map, ssp0, ssp1, ssp2);
			}
		
			positions_offset += models->mbs_positions_counts[mb_index];
		}
	}
}
*/