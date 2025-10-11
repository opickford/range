#ifndef VECTOR4_H
#define VECTOR4_H

#include "utils/str_utils.h"

// TODO: Comments

typedef struct
{
	float x, y, z, w;

} v4_t;

inline void v4_mul_eq_v4(v4_t* v0, v4_t v1)
{
	v0->x *= v1.x;
	v0->y *= v1.y;
	v0->z *= v1.z;
	v0->w *= v1.w;
}

inline void v4_mul_eq_f(v4_t* v, float f)
{
	v->x *= f;
	v->y *= f;
	v->z *= f;
	v->w *= f;
}

inline v4_t v4_mul_f(v4_t v, float f)
{
	return (v4_t) { v.x * f, v.y * f, v.z * f, v.w * f };
}

inline void v4_add_eq_f(v4_t* v, float f)
{
	v->x += f;
	v->y += f;
	v->z += f;
	v->w += f;
}

inline v4_t v4_add_f(v4_t v, float f)
{
	return (v4_t) { v.x + f, v.y + f, v.z + f, v.w + f };
}

inline void v4_add_eq_v4(v4_t* v0, v4_t v1)
{
	v0->x += v1.x;
	v0->y += v1.y;
	v0->z += v1.z;
	v0->w += v1.w;
}

inline void v4_sub_eq_f(v4_t* v, float f)
{
	v->x -= f;
	v->y -= f;
	v->z -= f;
	v->w -= f;
}

inline void v4_sub_eq_v4(v4_t* v0, v4_t v1)
{
	v0->x -= v1.x;
	v0->y -= v1.y;
	v0->z -= v1.z;
	v0->w -= v1.w;
}

inline v4_t v4_sub_v4(v4_t v0, v4_t v1)
{
	return (v4_t) { v0.x - v1.x, v0.y - v1.y, v0.z - v1.z, v0.w - v1.w };
}

// TODO: If it turns out these are taking a noticeable amount of CPU,
//		 convert to #define.
inline void v4_write(float* out, v4_t v)
{
	out[0] = v.x;
	out[1] = v.y;
	out[2] = v.z;
	out[3] = v.w;
}

// TODO: Feels a little messy
inline void v4_write_xyz(float* out, v4_t v)
{
	out[0] = v.x;
	out[1] = v.y;
	out[2] = v.z;
}

inline v4_t v4_read(const float* in)
{
	return (v4_t) { in[0], in[1], in[2], in[3] };
}

inline char* v4_to_str(v4_t v)
{
	return format_str("%f %f %f %f", v.x, v.y, v.z, v.w);
}

#endif