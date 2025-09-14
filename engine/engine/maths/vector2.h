#ifndef VECTOR2_H
#define VECTOR2_H

typedef struct
{
	float x, y;

} V2;

inline V2 v2_add_v2(V2 v0, V2 v1)
{
    return (V2) { v0.x + v1.x, v0.y + v1.y };
}

inline V2 v2_sub_v2(V2 v0, V2 v1)
{
    return (V2) { v0.x - v1.x, v0.y - v1.y };
}

inline V2 v2_mul_f(V2 v, float f)
{
    return (V2) { v.x * f, v.y * f };
}

inline V2 v2_read(const float* in)
{
	return (V2) { in[0], in[1] };
}

#endif