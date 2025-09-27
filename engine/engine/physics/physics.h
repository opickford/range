#ifndef PHYSICS_H
#define PHYSICS_H

#include "maths/vector3.h"

#include <cecs/ecs.h>

// TODO: Move to separate file?
typedef struct
{
    V3 force;
    V3 velocity;

    float mass;
} PhysicsData;

void PhysicsData_init(PhysicsData* data);

void Physics_tick(ECS* ecs, System* physics_system, System* collision_system, float dt);

#endif