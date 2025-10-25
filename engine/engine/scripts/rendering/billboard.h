#ifndef BILLBOARD_H
#define BILLBOARD_H

#include <core/engine.h>

inline void update_billboard(engine_t* engine, cecs_entity_id_t eid, float dt)
{
    // TODO: Should be a billboard tag component which you add to an entity, then a system does this for you.
    // TODO: Could be an engine functon to create a billboard entity.
    // TODO: We really want to support transparency in textures for this.


    transform_t* transform = cecs_add_component(engine->ecs, eid, COMPONENT_TRANSFORM);

    v3_t dir = v3_sub_v3(transform->position, engine->renderer.camera.position);
    v3_normalise(&dir);

    direction_to_eulers(dir, &transform->rotation.x, &transform->rotation.y);
    transform->rotation.x = 0;
}

#endif