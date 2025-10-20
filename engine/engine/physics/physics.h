#ifndef PHYSICS_H
#define PHYSICS_H

#include "physics_frame.h"

#include "core/scene.h"

#include "maths/vector3.h"
#include "maths/bounding_sphere.h"
#include "maths/plane.h"

#include "common/status.h"

#include <cecs/ecs.h>

#include <chds/vec.h>

#include <stdint.h>

// TODO: Rename -> PhysicsSystem?
typedef struct
{
    // TODO: Should this contain a scene also?

    cecs_t* ecs;

    physics_frame_t frame;

    // Views
    cecs_view_id_t physics_view; // TODO: Rename physicsdata view?
    cecs_view_id_t moving_colliders_view;
    cecs_view_id_t static_colliders_view;
    cecs_view_id_t colliders_view;

} physics_t;

// TODO: Move to separate file?
typedef struct physics_data
{
    v3_t impulses; // Forces applied instantaneously.
    v3_t velocity;

    float mass;


    // TODO: TEMP
    uint8_t floating;
} physics_data_t;

void physics_data_init(physics_data_t* data);

status_t physics_init(physics_t* physics, cecs_t* ecs);
void physics_tick(physics_t* physics, scene_t* scene, float dt);

typedef enum
{
    COLLISION_SHAPE_ELLIPSOID,
    COLLISION_SHAPE_MESH
} collision_shape_type_t;

typedef struct
{
    mesh_base_t* mb;
    chds_vec(v3_t) wsps;
} collision_mesh_t;

typedef struct
{
    collision_shape_type_t type;

    union
    {
        collision_mesh_t mesh;
        v3_t ellipsoid;
    };

    // TODO: How do we ensure that these are set???? Just down to user???
    uint8_t dirty; // Recalculate world space positions
    uint8_t scale_dirty; // Recalculate bounding sphere radius

    // TODO: For broad phase. But surely we don't need this if the type isn't a Mesh?
    bounding_sphere_t bs;

} collision_shape_t;

typedef struct collider
{
    // TODO: In the future this could contain some callback etc.
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

#endif