#ifndef PLANE_H
#define PLANE_H

#include "vector3.h"

typedef struct
{
	V3 normal;
	V3 point;
} Plane;


float signed_distance(const Plane* plane, const V3 point);

float line_intersect_plane(const V3 v0, const V3 v1, const Plane* plane, V3* out);

#endif