#ifndef UTILS_H
#define UTILS_H

#define _USE_MATH_DEFINES
#include <math.h>

#define PI (float)M_PI
#define PI_2 (float)M_PI_2

// TODO: Should this file be renamed to differentiate from utils folder?

inline float radians(float degrees)
{
	return (float)(degrees * PI) / 180.f;
}

inline float lerp(float a, float b, float t)
{
	return a + (b - a) * t;
}

inline void direction_to_eulers(const V3 direction, float* pitch, float* yaw)
{
	// Converts a direction to its euler angles, roll is 0 for a direction.


	/*
	camera->direction[0] = sinf(camera->yaw) * cosPitch;
    camera->direction[1] = sinf(camera->pitch);
    camera->direction[2] = cosf(camera->yaw) * cosPitch;


	asin(dir[0] / cosPitch) = acos(dir[2] / cosPitch)
	atan2f(

	pitch = asin(direction[1])
	cos(yaw) = acos(dir[2] / cosPitch)
	*/

	// Pitch
	//*pitch = atan2f(direction[1], sqrtf(direction[0] * direction[0] + direction[2] * direction[2]));

	// Yaw
	//*yaw = -atan2f(direction[0], -direction[2]);

	

}

#endif