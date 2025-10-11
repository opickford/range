#ifndef VECTOR2_H
#define VECTOR2_H

typedef struct
{
	float x, y;

} v2_t;

inline v2_t v2_add_v2(v2_t v0, v2_t v1)
{
    return (v2_t) { v0.x + v1.x, v0.y + v1.y };
}

inline v2_t v2_sub_v2(v2_t v0, v2_t v1)
{
    return (v2_t) { v0.x - v1.x, v0.y - v1.y };
}

inline v2_t v2_mul_f(v2_t v, float f)
{
    return (v2_t) { v.x * f, v.y * f };
}

inline v2_t v2_read(const float* in)
{
	return (v2_t) { in[0], in[1] };
}

#endif