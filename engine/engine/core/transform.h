#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "maths/vector3.h"

typedef struct
{
    V3 position;
    V3 scale;
    V3 rotation;
} Transform;

#endif