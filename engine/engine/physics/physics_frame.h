#pragma once

#ifndef PHYSICS_FRAME_H
#define PHYSICS_FRAME_H

#include <maths/vector3.h>

#include <cecs/archetype.h>

#include <chds/vec.h>


#include <string.h>

typedef struct collider collider_t;
typedef struct physics_data physics_data_t;
typedef struct mesh_instance mesh_instance_t;
typedef struct transform transform_t;

// TODO: Should this just store the necessary for each entity, i think so? Would make it easier to unpack right and less computations....
// TODO: An issue with pointers would be if a collision called a 
//       callback that created a new entity which could invalidate 
//       the pointers maybe? But not an issue for now. I believe 
//       we would want to wait for all collisions to be resolved
//       before firing callbacks anyways.

// TODO: Should this be elsewhere?
typedef struct
{
    mesh_instance_t* mi0;
    collider_t* c0;
    physics_data_t* pd0;
    transform_t* t0;

    mesh_instance_t* mi1;
    collider_t* c1;
    physics_data_t* pd1;
    transform_t* t1;
    



    /*
    // collider_t collides with a target.
    cecs_archetype_id_t collider_aid;
    int collider_offset;

    cecs_archetype_id_t target_aid;
    int target_offset;
    */

} potential_collision_t;

typedef struct
{
    v3_t rel_vel;
    v3_t collision_normal;

    uint8_t hit;

    potential_collision_t pc; // TODO: TEMP: Just for the entity ptrs?

} collision_data_t;


//typedef struct
//{
//    int temp;
//} collision_t;

typedef struct
{
    chds_vec(potential_collision_t) broad_phase_collisions;
    // float dt;
    // TODO: Should this have dt? Simplifies passing args but not needed in some places.
    // TODO: Store collisions for discrete collision resolution.

    chds_vec(collision_data_t) collisions;

} physics_frame_t;

inline void physics_frame_init(physics_frame_t* pf)
{
    memset(pf, 0, sizeof(physics_frame_t));
}

inline void physics_frame_destroy(physics_frame_t* pf)
{
    chds_vec_destroy(pf->broad_phase_collisions);
}

#endif