#ifndef BOUNDING_SPHERE_H
#define BOUNDING_SPHERE_H

#include "vector3.h"

typedef struct bounding_sphere
{
    v3_t centre;
    float radius;

} bounding_sphere_t;

#endif