#ifndef VECTOR_MATHS_H
#define VECTOR_MATHS_H

#include "vector3.h"
#include "vector4.h"

#include "utils/common.h"

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

inline uint8_t point_in_triangle(v3_t p, v3_t a, v3_t b, v3_t c)
{
    // FROM: https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html

    // Translate ABC so P is at the origin, then check if world 
    // origin is contained in ABC.
    v3_sub_eq_v3(&a, p);
    v3_sub_eq_v3(&b, p);
    v3_sub_eq_v3(&c, p);

    // Calculate normals of triangles using two vertices and the 
    // origin (P).
    v3_t u = cross(b, c); // PBC
    v3_t v = cross(c, a); // PCA
    v3_t w = cross(a, b); // PAB

    // If all normals face the same direction, then ABC contains P.
    if (dot(u, v) < 0.f) return 0;
    if (dot(u, w) < 0.f) return 0;

    return 1;
}

#endif