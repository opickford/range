#include "plane.h"

#include "vector3.h"

#include <stdio.h>

// TODO: Doesn't need to be passed by pointer but i think old code using it 
//       has ptrs already.
float signed_distance(const plane_t* plane, const v3_t point)
{
	float d = -dot(plane->normal, plane->point);
	return dot(point, plane->normal) + d;
}

float line_intersect_plane(const v3_t v0, const v3_t v1, const plane_t* plane, v3_t* out)
{
	// This uses the fact that a plane can be expressed as the set of points p for which Dot((p - p0), n) = 0
	v3_t ray = v3_sub_v3(v1, v0);

	float normalDotRay = dot(plane->normal, ray);

    // Check if planes are parallel, this should not happen as we're using this
    // for the triangle plane intersection which means the triangle edge does
    // not intersect the plane...
	if (normalDotRay == 0)
	{
		// TODO: Removing this could improve performance a little.
		printf("normal_dot_ray == 0. Should not happen\n");
		return 0;
	}

	v3_t v0_to_plane_point = v3_sub_v3(v0, plane->point);

	// Calculate the time of intersection.
	float t = -(dot(plane->normal, v0_to_plane_point)) / normalDotRay;
	
	// Interpolate for the point of intersection.
	*out = v3_mul_f(ray, t);
	v3_add_eq_v3(out, v0);

	return t;
}

plane_t plane_from_points(const v3_t v0, const v3_t v1, const v3_t v2)
{
    return (plane_t) {
        .normal = v3_normalised(cross(v3_sub_v3(v1, v0), v3_sub_v3(v2, v0))),
        .point = v1
    };
}