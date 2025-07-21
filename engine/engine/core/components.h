#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <cecs/ecs.h>

// TODO: Is XMacro pattern necessary here?


#define CORE_COMPONENTS_LIST \
    X(MeshInstance) \
    X(PointLight)

// TODO: Should this be like COMPONENT_ID_...?
#define X(T) ComponentID COMPONENT_##T;
    CORE_COMPONENTS_LIST
#undef X

inline void CoreComponents_init(ECS* ecs)
{       
    // Register all core engine components.
#define X(T) COMPONENT_##T = ECS_register_component(ecs, sizeof(T)); printf("Registered Component: " #T "id: %d\n", COMPONENT_##T);
    CORE_COMPONENTS_LIST            
#undef X
}


#endif