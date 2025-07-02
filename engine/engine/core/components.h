#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <cecs/ecs.h>

// TODO: Is XMacro pattern necessary here?


#define CORE_COMPONENTS_LIST \
    X(MeshInstance)

#define X(T) ComponentID COMPONENT_##T;
    CORE_COMPONENTS_LIST
#undef X

inline void CoreComponents_init(ECS* ecs)
{
        
    // Register all core engine components.
#define X(T) COMPONENT_##T = ECS_register_component(ecs, sizeof(T)); printf("Registered Component: " #T "id: %d\n", COMPONENT_##T);
    CORE_COMPONENTS_LIST            
#undef X

        //COMPONENT_MeshInstance = ECS_register_component(ecs, sizeof(COMPONENT_MeshInstance)); printf("Registered Component: " "MeshInstance" "id: %d\n", COMPONENT_MeshInstance);

}


#endif