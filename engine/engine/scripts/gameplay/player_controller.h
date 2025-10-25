#ifndef PLAYER_CONTROLLER_H
#define PLAYER_CONTROLLER_H

#include <core/engine.h>

#include <cecs/ecs.h>

uint8_t third_person = 1;

inline void player_controller(engine_t* engine, cecs_entity_id_t eid, float dt)
{

    if (engine->input_mode != INPUT_MODE_GAME) return;

    /*
    TODO: How can this be customisable?
    
    We really want this player controller to be provided by the engine.
    - Where/how can we provide it?

    Note, using this we see the choppy movement due to less physics updates,
    therefore, we need to lerp between positions for all entities.


    TODO:
    - Could accept an entity id.
    - Would this control the mouse as well?? 
    - This shouldn't be based off a camera direction right? more a player direction?
    
    */

    const static v3_t up = { 0, 1.f, 0 };
    const v3_t forward = v3_normalised((v3_t) { engine->renderer.camera.direction.x, 0.f, engine->renderer.camera.direction.z });
    const v3_t right = v3_normalised(cross(forward, up));

    physics_data_t* pd = cecs_get_component(engine->ecs, eid, COMPONENT_PHYSICS_DATA);
    transform_t* t = cecs_get_component(engine->ecs, eid, COMPONENT_TRANSFORM);

    // By multiplying by dt we're essentially converting the force into an impulse.
    const float speed = 20.f * dt * pd->mass;

    if (CSRGE_KEYDOWN(engine->window.keys['W']))
    {    
        v3_add_eq_v3(&pd->impulses, v3_mul_f(forward, speed));
    }
    if (CSRGE_KEYDOWN(engine->window.keys['S']))
    {
        v3_sub_eq_v3(&pd->impulses, v3_mul_f(forward, speed));
    }
    if (CSRGE_KEYDOWN(engine->window.keys['A']))
    {
        v3_sub_eq_v3(&pd->impulses, v3_mul_f(right, speed));
    }
    if (CSRGE_KEYDOWN(engine->window.keys['D']))
    {   
        v3_add_eq_v3(&pd->impulses, v3_mul_f(right, speed));
    }
    if (CSRGE_KEYDOWN(engine->window.keys[' ']))
    {

        //collider_t* c = cecs_get_component(engine->ecs, player_entity, COMPONENT_COLLIDER);

        
        // TODO: Only if colliding with something???
        const float jump_height = pd->mass;
        v3_add_eq_v3(&pd->impulses, v3_mul_f(up, jump_height));   
    }

    v3_t player_pos = v3_lerp(t->previous_position, t->position, engine->renderer.frame_data.physics_alpha);
    if (third_person)
    {
        const static float cam_dist = 4.f;
        const static float lateral_offset = 2.f;
        const static float vertical_offset = 2.f;


        v3_t pos = v3_sub_v3(player_pos, v3_mul_f(engine->renderer.camera.direction, cam_dist));
        v3_add_eq_v3(&pos, v3_mul_f(right, lateral_offset));
        v3_add_eq_v3(&pos, v3_mul_f(up, vertical_offset));
        engine->renderer.camera.position = pos;
    }
    else
    {
        engine->renderer.camera.position = player_pos;
    }
    
}

#endif