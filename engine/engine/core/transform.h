#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "maths/vector3.h"

// TODO: Should this just go into components?

typedef struct
{
    v3_t position;
    v3_t scale;
    v3_t rotation;
} transform_t;

inline void transform_init(transform_t* t)
{
    t->position = (v3_t){ 0.f, 0.f, 0.f };
    t->scale = (v3_t){ 1.f,1.f,1.f };
    t->rotation = (v3_t){ 0.f, 0.f, 0.f };
}

#endif