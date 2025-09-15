#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "maths/vector3.h"

// TODO: Should this just go into components?

typedef struct
{
    V3 position;
    V3 scale;
    V3 rotation;
} Transform;

inline void Transform_init(Transform* t)
{
    t->position = (V3){ 0.f, 0.f, 0.f };
    t->scale = (V3){ 1.f,1.f,1.f };
    t->rotation = (V3){ 0.f, 0.f, 0.f };
}

#endif