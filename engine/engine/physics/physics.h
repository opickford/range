#ifndef PHYSICS_H
#define PHYSICS_H

#include "core/scene.h"

#include "maths/vector3.h"
#include "maths/bounding_sphere.h"

#include "utils/vector.h"

#include <cecs/ecs.h>

#include <stdint.h>

// TODO: Move to separate file?
typedef struct
{
    V3 force;
    V3 velocity;

    float mass;
} PhysicsData;

void PhysicsData_init(PhysicsData* data);

void Physics_tick(ECS* ecs, Scene* scene, System* physics_system, System* collision_system, float dt);


typedef struct
{
    // TODO: How do we ensure that these are set???? Just down to user???
    uint8_t dirty; // Recalculate world space positions
    uint8_t scale_dirty; // Recalculate bounding sphere radius

    Vector(V3) wsps;
    BoundingSphere bs;

} CollisionCache;

void CollisionCache_init(CollisionCache* cc);
void CollisionCache_destroy(CollisionCache* cc);

#endif