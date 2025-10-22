#ifndef COLLISION_H
#define COLLISION_H

#include "maths/bounding_sphere.h"

#include <chds/vec.h>

typedef struct v3 v3_t;
typedef struct mesh_base mesh_base_t;
typedef struct physics physics_t;
typedef struct scene scene_t;

typedef enum
{
    COLLISION_SHAPE_ELLIPSOID,
    COLLISION_SHAPE_MESH
} collision_shape_type_t;

typedef struct
{
    const mesh_base_t* mb;
    chds_vec(v3_t) wsps;
} collision_mesh_t;

typedef struct
{
    // Tagged union for the type narrow phase shape.
    collision_shape_type_t type;
    union
    {
        collision_mesh_t mesh;
        v3_t ellipsoid;
        float radius;
    };

    // TODO: How do we ensure that these are set???? Just down to user???
    uint8_t dirty; // Recalculate world space positions
    uint8_t scale_dirty; // Recalculate bounding sphere radius

    // Broad phase shape.
    bounding_sphere_t bs;

} collision_shape_t;

// TODO: In the future this could contain some callback etc.
typedef struct collider
{
    collision_shape_t shape;

    // Ratio of relative velocity of separation to relative velocity of approach.
    // Determines collisions elasticity (how much energy loss)
    // 0 = completely inelastic, 1 = perfectly elastic (no energy loss)
    float restiution_coeff;

    // Ratio of frictional force to normal force pushing objects together.
    // 0 = no friction, 1 is as much friction as the normal force.
    float friction_coeff;

} collider_t;

void collider_init(collider_t* c);
void collider_destroy(collider_t* c);

void handle_collisions(physics_t* physics, scene_t* scene);

#endif