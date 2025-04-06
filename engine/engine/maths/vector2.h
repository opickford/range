#ifndef VECTOR2_H
#define VECTOR2_H

#include "utils/str_utils.h"

typedef struct
{
	float x, y;

} V2;

inline void v2_swap(V2* v0, V2* v1)
{
	V2 temp = *v0;
	*v0 = *v1;
	*v1 = temp;
}

inline V2 v2_read(const float* in)
{
	return (V2) { in[0], in[1] };
}

inline V2 v2_mul_f(V2 v, float f)
{
	return (V2) { v.x * f, v.y * f };
}

inline char* v2_to_str(V2 v)
{
	return format_str("%f %f", v.x, v.y);
}

inline float v2_size(V2 v)
{
	return sqrtf(v.x * v.x + v.y * v.y);
}

inline void v2_normalise(V2* v)
{
	const float size = v2_size(*v);
	if (size == 0) return;

	const float inv_size = 1.f / size;
	v->x *= inv_size;
	v->y *= inv_size;
}

inline V2 v2_sub_v2(V2 v0, V2 v1)
{
	return (V2) { v0.x - v1.x, v0.y - v1.y };
}

#endif