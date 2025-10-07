#ifndef COMPONENTS_H
#define COMPONENTS_H

#include "transform.h"
#include "mesh_instance.h"
#include "lights.h"

#include "physics/physics.h"

#include <cecs/ecs.h>

#define CORE_COMPONENTS_LIST          \
    X(MESH_INSTANCE, mesh_instance_t) \
    X(POINT_LIGHT, point_light_t)     \
    X(TRANSFORM, transform_t)         \
    X(PHYSICS_DATA, physics_data_t)   \
    X(COLLIDER, collider_t)

#define X(name, T) cecs_component_id_t COMPONENT_##name;
    CORE_COMPONENTS_LIST
#undef X

inline void core_components_init(cecs_t* ecs)
{
    // Register all core engine components.
#define X(name, T) \
    COMPONENT_##name = cecs_register_component(ecs, sizeof(T)); \
    printf("Registered Component: " #name " id: %d\n", COMPONENT_##name);

        CORE_COMPONENTS_LIST
#undef X
}


#endif