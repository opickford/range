#ifndef MATRIX4_H
#define MATRIX4_H

#include "vector4.h"
#include "vector3.h"

#include "utils/str_utils.h"

/*
- Memory Layout:
We're using a column major memory layout for better cache usage
when multiplying m4 and v4 as columns are stored contingous memory.
After testing, there actually seemed to be no real noticeable 
difference in performance though.

- Multiplication
We post-multiply (matrix first) a m4 and v4, this means that the v4
is a column vector (4row x 1col) matrix. This is because the cols
in the first matrix much match the rows in the second.


// TODO: I want to stop using column major storage as both opengl and directx do
		 not, therefore it's just confusing and weird. Right?.... 

		 Not sure it even matters. Can be a refactor in future if necessary.

*/
typedef float m4_t[16];

// TODO: Rename _out?
void m4_mul_m4(const m4_t m0, const m4_t m1, m4_t out);

void m4_mul_v4(const m4_t m, v4_t v, v4_t* out);

void m4_identity(m4_t out);

void m4_translation(v3_t position, m4_t out);

void m4_rotation(float pitch, float yaw, float roll, m4_t out);

void look_at(v3_t position, v3_t direction, m4_t out);

void m4_model_matrix(v3_t position, v3_t eulers, v3_t scale, m4_t out);

void m4_normal_matrix(v3_t eulers, v3_t scale, m4_t out);

void m4_transposed(const m4_t in, m4_t out);

void m4_copy_m3(const m4_t in, m4_t out);

void m4_projection(float fov, float aspect_ratio, float near_plane, float far_plane, m4_t out);

char* m4_to_str(const m4_t m);

#endif
