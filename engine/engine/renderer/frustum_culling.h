#ifndef FRUSTUM_CULLING_H
#define FRUSTUM_CULLING_H

#include "maths/plane.h"

#define MAX_FRUSTUM_PLANES 6

// TODO: Apparently this is the max 1 triangle can become, 
//       figure out or test for myself.
#define MAX_CLIPPED_TRIS_FACTOR 18

typedef struct
{
	plane_t planes[MAX_FRUSTUM_PLANES];
	int planes_count;

} view_frustum_t;

void view_frustum_init(view_frustum_t* view_frustum, float near_dist, float far_dist, float fov, float aspect_ratio);


#endif