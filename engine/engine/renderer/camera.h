#ifndef CAMERA_H
#define CAMERA_H

#include "maths/matrix4.h"
#include "maths/vector3.h"

typedef struct 
{
	float pitch;
	float yaw;
	float roll;

	v3_t direction;
	v3_t position;

} camera_t;

void calculate_view_matrix(const camera_t* camera, m4_t out);

#endif 