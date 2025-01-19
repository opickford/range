#include "frustum_culling.h"

#include "maths/utils.h"

#include <string.h>

void view_frustum_init(ViewFrustum* view_frustum, float near_dist, float far_dist, float fov, float aspect_ratio)
{
	// Reset the struct.
	memset(view_frustum, 0, sizeof(ViewFrustum));

	// Calculate the dimensions of the near and far planes.
	float double_tanf_half_fov = 2.f * tanf(radians(fov) / 2.f);

	float near_height = double_tanf_half_fov * near_dist;
	float near_width = near_height * aspect_ratio;

	float far_height = double_tanf_half_fov * far_dist;
	float far_width = far_height * aspect_ratio;

	// Define the forward direction.
	// TODO: Should these be defined globally somewhere? or at least in a function? like v3_world_up?
	// TODO: Also these are really world up as well. Not specific to view space. 
	// TODO: Should be defined in some engine_globals.h maybe.
	V3 view_forward = { 0, 0, -1.f };
	V3 view_up = { 0, 1.f, 0 };
	V3 view_right = { 1.f, 0, 0 };

	// Define offsets for the near and far planes.
	V3 near_centre = v3_mul_f(view_forward, near_dist);
	V3 near_top_offset = v3_mul_f(view_up, near_height * 0.5f);
	V3 near_right_offset = v3_mul_f(view_right, near_width * 0.5f);

	V3 far_centre = v3_mul_f(view_forward, far_dist);
	V3 far_top_offset = v3_mul_f(view_up, far_height * 0.5f);
	V3 far_right_offset = v3_mul_f(view_right, far_width * 0.5f);

	// To calculate each plane normal, we can define 3 points on each plane
	// and use the cross product of the edges. So define the four corners
	// of each plane.
	V3 near_top_left = v3_sub_v3(v3_add_v3(near_centre, near_top_offset), near_right_offset);
	V3 near_top_right = v3_add_v3(v3_add_v3(near_centre, near_top_offset), near_right_offset);
	V3 near_bottom_left = v3_sub_v3(v3_sub_v3(near_centre, near_top_offset), near_right_offset);
	V3 near_bottom_right = v3_add_v3(v3_sub_v3(near_centre, near_top_offset), near_right_offset);

	V3 far_top_left = v3_sub_v3(v3_add_v3(far_centre, far_top_offset), far_right_offset);
	V3 far_top_right = v3_add_v3(v3_add_v3(far_centre, far_top_offset), far_right_offset);
	//V3 far_bottom_left = v3_sub_v3(v3_sub_v3(far_centre, far_top_offset), far_right_offset);
	//V3 far_bottom_right = v3_add_v3(v3_sub_v3(far_centre, far_top_offset), far_right_offset);

	// Near and far are trivial to define.
	Plane near = {
		.point = near_centre,
		.normal = { 0, 0, -1.f}
	};

	Plane far = {
		.point = far_centre,
		.normal = { 0, 0, 1.f}
	};

	// Define the left/right planes, opposite x direction.

	// TODO: Function for calculating the plane normal.
	Plane right = { 
		.point = near_top_right, 
		.normal = normalised(cross(v3_sub_v3(near_top_right, near_bottom_right), v3_sub_v3(far_top_right, near_bottom_right)))
	};
	
	// Define the left plane.
	Plane left = { 
		.point = near_top_left,
		.normal = right.normal
	};
	left.normal.x *= -1;

	// Define the top/bottom planes, opposite y direction.
	Plane top = {
		.point = near_top_left,
		.normal = normalised(cross(v3_sub_v3(far_top_right, far_top_left), v3_sub_v3(near_top_left, far_top_left)))
	};

	Plane bottom = {
		.point = near_bottom_left,
		.normal = top.normal
	};
	bottom.normal.y *= -1;
	
	// As we're not doing screen space clipping, all are necessary except far,
	// however, if we're rendering stuff far away, far is a great optimisation.
	// Also, there is almost no real cost to enabling far even if we're not 
	// rendering anything far away because it's just a broad phase check that
	// will fail.
	view_frustum->planes[view_frustum->planes_count++] = near;
	view_frustum->planes[view_frustum->planes_count++] = far;
	view_frustum->planes[view_frustum->planes_count++] = right;
	view_frustum->planes[view_frustum->planes_count++] = left;
	view_frustum->planes[view_frustum->planes_count++] = top;
	view_frustum->planes[view_frustum->planes_count++] = bottom;
}