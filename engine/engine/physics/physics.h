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

typedef struct physics
{
    // TODO: Should this contain a scene also?

    cecs_t* ecs;

    physics_frame_t frame;

    // Views
    cecs_view_id_t physics_view; // TODO: Rename physicsdata view?
    cecs_view_id_t moving_colliders_view;
    cecs_view_id_t static_colliders_view;
    cecs_view_id_t colliders_view;

    uint8_t max_collision_iters;

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

#endif