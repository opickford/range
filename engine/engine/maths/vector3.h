#ifndef VECTOR3_H
#define VECTOR3_H

#include "utils/str_utils.h"

#include <math.h>

// TODO: Top of file comments explaining stuff properly?

/*
Originally VXs were defined as float[x], however, after profiling, there was no speed 
difference, most likely because they compile to the same thing.

Also, we pass these by value which should result in better performance for a small struct,
see: https://austinmorlan.com/posts/pass_by_value_vs_pointer/#:~:text=When%20I%20was%20in%20college,of%20that%20on%20the%20stack

Returning the v3_t seems to have no performance impact and makes the code much more readable.

No need to overcomplicate this stuff, I can always do it differently in a per-pixel loop if
profiling shows this is an issue.
*/

typedef struct v3
{
	float x, y, z;

} v3_t;

inline v3_t cross(v3_t v0, v3_t v1)
{
	return (v3_t) { v0.y * v1.z - v0.z * v1.y, v0.z * v1.x - v0.x * v1.z, v0.x * v1.y - v0.y * v1.x };
}

// TODO: Rename v3_size etc?
inline float v3_size(v3_t v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline float v3_size_sqrd(v3_t v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline void v3_normalise(v3_t* v)
{	
	const float inv_size = 1.f / sqrtf(v->x * v->x + v->y * v->y + v->z * v->z);
	v->x *= inv_size;
	v->y *= inv_size;
	v->z *= inv_size;
}

// TODO: v3_ prefix?
inline v3_t v3_normalised(v3_t v)
{
	const float inv_size = 1.f / v3_size(v);

	return (v3_t) { v.x * inv_size, v.y * inv_size, v.z * inv_size, };
}

inline float dot(v3_t v0, v3_t v1)
{
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

inline void v3_mul_eq_v3(v3_t* v0, v3_t v1)
{
	v0->x *= v1.x;
	v0->y *= v1.y;
	v0->z *= v1.z;
}

inline v3_t v3_mul_v3(v3_t v0, v3_t v1)
{
	return (v3_t) { v0.x * v1.x, v0.y * v1.y, v0.z * v1.z };
}

inline void v3_mul_eq_f(v3_t* v, float f)
{
	v->x *= f;
	v->y *= f;
	v->z *= f;
}

inline v3_t v3_mul_f(v3_t v, float f)
{
	return (v3_t) { v.x * f, v.y * f, v.z * f };
}

inline void v3_add_eq_f(v3_t* v, float f)
{
	v->x += f;
	v->y += f;
	v->z += f;
}

inline void v3_add_eq_v3(v3_t* v0, v3_t v1)
{
	v0->x += v1.x;
	v0->y += v1.y;
	v0->z += v1.z;
}

inline v3_t v3_add_v3(v3_t v0, v3_t v1)
{
	return (v3_t) { v0.x + v1.x, v0.y + v1.y, v0.z + v1.z };
}

inline void v3_sub_eq_f(v3_t* v, float f)
{
	v->x -= f;
	v->y -= f;
	v->z -= f;
}

inline void v3_sub_eq_v3(v3_t* v0, v3_t v1)
{
	v0->x -= v1.x;
	v0->y -= v1.y;
	v0->z -= v1.z;
}

inline v3_t v3_sub_v3(v3_t v0, v3_t v1)
{
	return (v3_t) { v0.x - v1.x, v0.y - v1.y, v0.z - v1.z };
}

inline v3_t v3_uniform(float n)
{
    return (v3_t) { n, n, n };
}

inline v3_t v3_inv(v3_t v)
{
    return (v3_t) { 1.f / v.x, 1.f / v.y, 1.f / v.z };
}

inline v3_t v3_lerp(v3_t v0, v3_t v1, float t)
{
    return v3_add_v3(v0, v3_mul_f(v3_sub_v3(v1, v0), t));
}

/*
inline void v3_write(float* out, v3_t v)
{
	out[0] = v.x;
	out[1] = v.y;
	out[2] = v.z;
}
*/



inline v3_t v3_read(const float* in)
{
	return (v3_t) {in[0], in[1], in[2] };
}


inline char* v3_to_str(v3_t v)
{
	return format_str("%f %f %f", v.x, v.y, v.z);
}

// TODO: This shouldn't be in vector3.h?
inline int is_front_face(v3_t v0, v3_t v1, v3_t v2)
{
    // Camera is at origin as coordinates are given in view space, so 
    // view vector from tri to cam is V = v0 - P or just -v0
    // 
    // dot(A,B) = |A||B|cos(theta), we only care about the sign, so
    // just becomes dot(A,B) = cos(theta).
    // 
    // Dot product of face normal and V gives cos(theta) between
    // vectors, cos(90) = 0, cos(<90) -> positive. Therefore, 
    // Front facing IF dot(-v0, normal) > 0
    // To avoid the extra inversion of v0, convert to:
    // Front facing IF dot(v0, normal) <= 0
	return dot(v0, cross(v3_sub_v3(v1, v0), v3_sub_v3(v2, v0))) <= 0;
}

#endif