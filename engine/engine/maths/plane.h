#ifndef PLANE_H
#define PLANE_H

#include "vector3.h"

typedef struct
{
	v3_t normal;
	v3_t point;
} plane_t;


float signed_distance(const plane_t* plane, const v3_t point);

float line_intersect_plane(const v3_t v0, const v3_t v1, const plane_t* plane, v3_t* out);

plane_t plane_from_points(const v3_t v0, const v3_t v1, const v3_t v2);

#endif