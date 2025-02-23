#ifndef VECTOR_MATHS_H
#define VECTOR_MATHS_H

#include "vector3.h"
#include "vector4.h"

// TODO: Per function comments.


// Avoid circular dependencies by having shared functions here.
inline V3 v4_xyz(V4 v)
{
	return (V3) { .x = v.x, v.y, v.z };
}

inline V4 v3_to_v4(V3 in, float w)
{
	return (V4) { in.x, in.y, in.z, w };
}

// Combine these to potentially minimise function calls.
inline V4 v3_read_to_v4(const float* in, float w)
{
	return (V4) { in[0], in[1], in[2], w };
}

inline float double_signed_area(V2 a, V2 b, V2 c)
{
	// Returns double the signed area of a triangle a,b,c.
	// Uses the shoelace formula with the abs and half omitted.
	return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

#endif