#ifndef FRUSTUM_CULLING_H
#define FRUSTUM_CULLING_H

#include "maths/plane.h"

#define MAX_FRUSTUM_PLANES 6

typedef struct
{
	Plane planes[MAX_FRUSTUM_PLANES];
	int planes_count;

} ViewFrustum;

void view_frustum_init(ViewFrustum* view_frustum, float near_dist, float far_dist, float fov, float aspect_ratio);


#endif