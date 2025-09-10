#ifndef VECTOR2_H
#define VECTOR2_H

typedef struct
{
	float x, y;

} V2;

inline void v2_swap(V2* v0, V2* v1)
{
    // TODO: Generic function for this?
	V2 temp = *v0;
	*v0 = *v1;
	*v1 = temp;
}

/*
inline V2 v2_read(const float* in)
{
	return (V2) { in[0], in[1] };
}*/

#endif