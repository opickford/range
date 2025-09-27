#include "physics.h"

#include "core/components.h"

#include "maths/vector3.h"

static void apply_forces(ECS* ecs, System* physics_system)
{
    // TODO: Air resistance

    // Disable gravity for now.
    static V3 acceleration = { 0, 0, 0 };
    //static V3 acceleration = { 0, -9.8f, 0 };

    for (int si = 0; si < physics_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = physics_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        PhysicsData* physics_datas = Archetype_get_component_list(archetype, COMPONENT_PhysicsData);
        Transform* transforms = Archetype_get_component_list(archetype, COMPONENT_Transform);

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            PhysicsData* physics_data = &physics_datas[i];
            Transform* transform = &transforms[i];

            // F = MA
            v3_add_eq_v3(&physics_data->force, v3_mul_f(acceleration, physics_data->mass));

            // Apply force.
            v3_add_eq_v3(&physics_data->velocity, v3_mul_f(physics_data->force, 1.f / physics_data->mass));

            // TODO: Function for zeroing/filling?
            physics_data->force = (V3){ 0.f, 0.f, 0.f };
        }
    }
}

static void apply_velocities(ECS* ecs, System* physics_system, float dt)
{
    for (int si = 0; si < physics_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = physics_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        PhysicsData* physics_datas = Archetype_get_component_list(archetype, COMPONENT_PhysicsData);
        Transform* transforms = Archetype_get_component_list(archetype, COMPONENT_Transform);

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            PhysicsData* physics_data = &physics_datas[i];
            Transform* transform = &transforms[i];

            // Update position with new velocity.
            v3_add_eq_v3(&transform->position, v3_mul_f(physics_data->velocity, dt));
        }
    }
}

static void handle_collisions(ECS* ecs, System* collision_system, float dt)
{
    for (int si = 0; si < collision_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = collision_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        PhysicsData* physics_datas = Archetype_get_component_list(archetype, COMPONENT_PhysicsData);
        Transform* transforms = Archetype_get_component_list(archetype, COMPONENT_Transform);
        Transform* mis = Archetype_get_component_list(archetype, COMPONENT_MeshInstance);

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            PhysicsData* physics_data = &physics_datas[i];
            Transform* transform = &transforms[i];
            MeshInstance* mi = &mis[i];



        }
    }
}

void PhysicsData_init(PhysicsData* data)
{
    memset(data, 0, sizeof(PhysicsData));

    // Mass of 0 will cause divide by zero error.
    data->mass = 1.f;
}

void Physics_tick(ECS* ecs, System* physics_system, System* collision_system, float dt)
{
    apply_forces(ecs, physics_system);

    handle_collisions(ecs, collision_system, dt);

    apply_velocities(ecs, physics_system, dt);
}
