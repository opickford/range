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

#endif