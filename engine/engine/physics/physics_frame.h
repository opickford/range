#pragma once

#include "utils/vector.h"


#include <cecs/archetype.h>

#include <string.h>

typedef struct
{
    // Collider collides with a target.
    ArchetypeID collider_aid;
    int collider_offset;

    ArchetypeID target_aid;
    int target_offset;

} PotentialCollision;

typedef struct
{
    // TODO: Here lies the issue with my vector implementation, no idea of how mcuh data taken.
    int num_potential_collisions;
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