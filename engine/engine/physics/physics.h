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

void Physics_tick(ECS* ecs, Scene* scene, ViewID physics_view, ViewID collision_view, float dt);

typedef enum
{
    COLLISION_SHAPE_ELLIPSOID,
    COLLISION_SHAPE_MESH
} CollisionShapeType;

typedef struct
{
    MeshBase* mb;
    Vector(V3) wsps;
} CollisionMesh;

// TODO: I want to refactor this into something like a Collider?
typedef struct
{
    CollisionShapeType type;

    union
    {
        CollisionMesh mesh;
        V3 ellipsoid;
    };

    // TODO: How do we ensure that these are set???? Just down to user???
    uint8_t dirty; // Recalculate world space positions
    uint8_t scale_dirty; // Recalculate bounding sphere radius

    // TODO: For broad phase. But surely we don't need this if the type isn't a Mesh?
    BoundingSphere bs;

} CollisionShape;

typedef struct
{
    // TODO: In the future this could contain some callback etc.
    CollisionShape shape;
} Collider;

void Collider_init(Collider* c);
void Collider_destroy(Collider* c);

#endif