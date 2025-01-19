#ifndef CAMERA_H
#define CAMERA_H

#include "maths/matrix4.h"
#include "maths/vector3.h"

typedef struct 
{
	float pitch;
	float yaw;
	float roll;

	V3 direction;
	V3 position;

} Camera;

void calculate_view_matrix(const Camera* camera, M4 out);

#endif 