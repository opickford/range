#pragma once

#ifndef PHYSICS_FRAME_H
#define PHYSICS_FRAME_H

#include <cecs/archetype.h>

#include <chds/vector.h>

#include <string.h>

// TODO: Should this just store the necessary for each entity, i think so? Would make it easier to unpack right and less computations....
// TODO: An issue with pointers would be if a collision called a 
//       callback that created a new entity which could invalidate 
//       the pointers maybe? But not an issue for now.

typedef struct
{
    //MeshInstance* mi0;
    //Collider* c0;
    //PhysicsData* pd0;


    // Collider collides with a target.
    ArchetypeID collider_aid;
    int collider_offset;

    ArchetypeID target_aid;
    int target_offset;

} PotentialCollision;

typedef struct
{
    Vector(PotentialCollision) broad_phase_collisions;

} PhysicsFrame;

inline void PhysicsFrame_init(PhysicsFrame* pf)
{
    memset(pf, 0, sizeof(PhysicsFrame));
}

inline void PhysicsFrame_destroy(PhysicsFrame* pf)
{
    Vector_destroy(pf->broad_phase_collisions);
}

#endif