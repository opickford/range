#ifndef VECTOR_MATHS_H
#define VECTOR_MATHS_H

#include "vector3.h"
#include "vector4.h"

// TODO: Per function comments.
// Avoid circular dependencies by having shared functions here.

inline v3_t v4_xyz(v4_t v)
{
	return (v3_t) { .x = v.x, v.y, v.z };
}

inline v4_t v3_to_v4(v3_t in, float w)
{
	return (v4_t) { in.x, in.y, in.z, w };
}

// Combine these to potentially minimise function calls.
inline v4_t v3_read_to_v4(const float* in, float w)
{
	return (v4_t) { in[0], in[1], in[2], w };
}

#endif