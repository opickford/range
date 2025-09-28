#include "physics.h"

#include "core/components.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"
#include "maths/vector_maths.h"


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

static void update_collision_caches(ECS* ecs, Scene* scene, System* collision_system)
{
    // TODO: Eventually we might want specific collision meshes to simplify it?
    // TODO: E.g. we don't need the higher detail meshes to collide with, but not needed for now.
    for (int si = 0; si < collision_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = collision_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        const Transform* transforms = Archetype_get_component_list(archetype, COMPONENT_Transform);
        const MeshInstance* mis = Archetype_get_component_list(archetype, COMPONENT_MeshInstance);
        CollisionCache* collision_caches = Archetype_get_component_list(archetype, COMPONENT_CollisionCache);

        // Update world space positions of entity, update centre of bounding sphere.
        for (int i = 0; i < archetype->entity_count; ++i)
        {
            CollisionCache* collision_cache = &collision_caches[i];

            // TODO: Where can we actually set this?????? after appling forces???
            //if (!collision_cache->dirty)
            //{
            //    continue;
            //}

            const Transform transform = transforms[i];
            const MeshInstance* mi = &mis[i];
            const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];

            Vector_reserve(collision_cache->wsps, mb->num_positions);

            // Calculate model matrix.
            // TODO: rotation or direction or eulers?
            M4 model_matrix;
            m4_model_matrix(transform.position, transform.rotation, transform.scale, model_matrix);

            // Update the centre of the bounding sphere.
            V4 world_space_centre;
            m4_mul_v4(
                model_matrix,
                v3_to_v4(mb->centre, 1.f),
                &world_space_centre);

            // Write out the bounding sphere.
            collision_cache->bs.centre = v4_xyz(world_space_centre);

            // Convert object space positions to world space.
            for (int j = 0; j < mb->num_positions; ++j)
            {
                const V4 osp = v3_to_v4(mb->object_space_positions.data[j], 1.f);
                V4 wsp;
                m4_mul_v4(model_matrix, osp, &wsp);

                collision_cache->wsps.data[j].x = wsp.x;
                collision_cache->wsps.data[j].y = wsp.y;
                collision_cache->wsps.data[j].z = wsp.z;
            }

            collision_cache->dirty = 0;
        }

        // Update bounding sphere radius
        for (int i = 0; i < archetype->entity_count; ++i)
        {
            CollisionCache* collision_cache = &collision_caches[i];

            if (!collision_cache->scale_dirty)
            {
                continue;
            }

            const Transform transform = transforms[i];
            const MeshInstance* mi = &mis[i];
            const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];

            BoundingSphere* bs = &collision_cache->bs;
            V3 centre = bs->centre;

            // Calculate the new radius of the mi's bounding sphere.
            float radius_squared = -1;

            for (int j = 0; j < mb->num_positions; ++j)
            {
                V3 v = collision_cache->wsps.data[j];
                V3 between = v3_sub_v3(v, centre);

                radius_squared = max(size_squared(between), radius_squared);
            }

            bs->radius = sqrtf(radius_squared);

            collision_cache->scale_dirty = 0;
        }
    }
}

static void handle_collisions(ECS* ecs, Scene* scene, System* collision_system, float dt)
{
    update_collision_caches(ecs, scene, collision_system);
    
    for (int si = 0; si < collision_system->num_archetypes; ++si)
    {
        const ArchetypeID archetype_id = collision_system->archetype_ids[si];
        Archetype* archetype = &ecs->archetypes[archetype_id];

        PhysicsData* physics_datas = Archetype_get_component_list(archetype, COMPONENT_PhysicsData);
        Transform* transforms = Archetype_get_component_list(archetype, COMPONENT_Transform);
        const MeshInstance* mis = Archetype_get_component_list(archetype, COMPONENT_MeshInstance);
        const CollisionCache* collision_caches = Archetype_get_component_list(archetype, COMPONENT_CollisionCache);

        for (int i = 0; i < archetype->entity_count; ++i)
        {
            PhysicsData* physics_data = &physics_datas[i];
            
            // If entity not moving, it cannot have collided with something.
            if (physics_data->velocity.x == 0 &&
                physics_data->velocity.y == 0 &&
                physics_data->velocity.z == 0)
            {
                continue;
            }

            // TODO: Calculate potential collisions broad phase??

            Transform* transform = &transforms[i];
            const MeshInstance* mi = &mis[i];
            const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
            const CollisionCache* collision_cache = &collision_caches[i];
            
            // TODO: This iteration is getting painfully messy. Look into iterators. or some macro

            // TODO: Broad phase check bounding sphere vs bounding sphere - change the colours of the mesh or something so i can visualise this.
            //       have it so one moves into another.
            for (int si1 = 0; si1 < collision_system->num_archetypes; ++si1)
            {
                const ArchetypeID archetype_id1 = collision_system->archetype_ids[si1];
                Archetype* archetype1 = &ecs->archetypes[archetype_id1];

                const MeshInstance* mis1 = Archetype_get_component_list(archetype1, COMPONENT_MeshInstance);
                const CollisionCache* collision_caches1 = Archetype_get_component_list(archetype1, COMPONENT_CollisionCache);

                for (int i = 0; i < archetype1->entity_count; ++i)
                {
                    MeshInstance* mi1 = &mis1[i];

                    // Don't collide with self.
                    if (mi1 == mi) continue;

                    const CollisionCache* collision_cache1 = &collision_caches1[i];
                    
                    const BoundingSphere bs0 = collision_cache->bs;
                    const BoundingSphere bs1 = collision_cache1->bs;

                    const float dist = size_squared(v3_sub_v3(bs1.centre, bs0.centre));
                    const float n = (bs0.radius + bs1.radius) * (bs0.radius + bs1.radius);

                    if (dist <= n)
                    {
                        // TODO: Collided!!! save to list or something?
                    }



                }
            }
            




            // TODO: We already calculate boudning sphere i n

            

        }
    }
}

void PhysicsData_init(PhysicsData* data)
{
    memset(data, 0, sizeof(PhysicsData));

    // Mass of 0 will cause divide by zero error.
    data->mass = 1.f;
}

void Physics_tick(ECS* ecs, Scene* scene, System* physics_system, System* collision_system, float dt)
{
    apply_forces(ecs, physics_system);

    handle_collisions(ecs, scene, collision_system, dt);

    apply_velocities(ecs, physics_system, dt);
}

void CollisionCache_init(CollisionCache* cc)
{
    memset(cc, 0, sizeof(CollisionCache));
    cc->dirty = 1;
    cc->scale_dirty = 1;
}

void CollisionCache_destroy(CollisionCache* cc)
{
    Vector_destroy(cc->wsps);
}
