#include "physics.h"

#include "core/components.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"
#include "maths/vector_maths.h"

#include "utils/logger.h"

// TODO: Reorganise functions here.

static void Physics_setup_views(Physics* physics)
{
    physics->physics_view = ECS_view(physics->ecs,
        COMPONENT_ID_TO_BITSET(COMPONENT_PhysicsData) | 
        COMPONENT_ID_TO_BITSET(COMPONENT_Transform),
        0);

    physics->colliders_view = ECS_view(physics->ecs,
        COMPONENT_ID_TO_BITSET(COMPONENT_Transform) | 
        COMPONENT_ID_TO_BITSET(COMPONENT_MeshInstance) | 
        COMPONENT_ID_TO_BITSET(COMPONENT_Collider),
        0);

    // TODO: Currently these views a MeshInstance is this correct?

    // Note, may not actually be moving, just has physicsdata so COULD be moving.
    physics->moving_colliders_view = ECS_view(physics->ecs,
        COMPONENT_ID_TO_BITSET(COMPONENT_PhysicsData) | 
        COMPONENT_ID_TO_BITSET(COMPONENT_Transform) |
        COMPONENT_ID_TO_BITSET(COMPONENT_MeshInstance) | 
        COMPONENT_ID_TO_BITSET(COMPONENT_Collider),
        0);

    // Static means no physics data.
    physics->static_colliders_view = ECS_view(physics->ecs,
        COMPONENT_ID_TO_BITSET(COMPONENT_Transform) | 
        COMPONENT_ID_TO_BITSET(COMPONENT_MeshInstance) | 
        COMPONENT_ID_TO_BITSET(COMPONENT_Collider),
        COMPONENT_ID_TO_BITSET(COMPONENT_PhysicsData));
}

Status Physics_init(Physics* physics, ECS* ecs)
{
    memset(physics, 0, sizeof(Physics));
    physics->ecs = ecs;

    Physics_setup_views(physics);

    PhysicsFrame_init(&physics->frame);

    return STATUS_OK;
}

static void apply_forces(Physics* physics)
{
    // TODO: Air resistance

    // Disable gravity for now.
    static V3 acceleration = { 0, 0, 0 };
    //static V3 acceleration = { 0, -9.8f, 0 };

    ViewIter it = ECS_view_iter(physics->ecs, physics->physics_view);

    while (ECS_view_iter_next(&it))
    {
        PhysicsData* physics_datas = ECS_get_column(it, COMPONENT_PhysicsData);
        Transform* transforms = ECS_get_column(it, COMPONENT_Transform);

        for (int i = 0; i < it.num_entities; ++i)
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

static void apply_velocities(Physics* physics, float dt)
{
    ViewIter it = ECS_view_iter(physics->ecs, physics->physics_view);
    while (ECS_view_iter_next(&it))
    {
        PhysicsData* physics_datas = ECS_get_column(it, COMPONENT_PhysicsData);
        Transform* transforms = ECS_get_column(it, COMPONENT_Transform);

        for (int i = 0; i < it.num_entities; ++i)
        {
            PhysicsData* physics_data = &physics_datas[i];
            Transform* transform = &transforms[i];

            // Update position with new velocity.
            v3_add_eq_v3(&transform->position, v3_mul_f(physics_data->velocity, dt));
        }
    }
}

static void update_collision_mesh_bounding_sphere(Collider* collider, const MeshInstance* mi, const Transform transform)
{
    const MeshBase* mb = collider->shape.mesh.mb;

    Vector_reserve(collider->shape.mesh.wsps, mb->num_positions);

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
    collider->shape.bs.centre = v4_xyz(world_space_centre);

    // Convert object space positions to world space.
    for (int j = 0; j < mb->num_positions; ++j)
    {
        const V4 osp = v3_to_v4(mb->object_space_positions[j], 1.f);
        V4 wsp;
        m4_mul_v4(model_matrix, osp, &wsp);

        collider->shape.mesh.wsps[j].x = wsp.x;
        collider->shape.mesh.wsps[j].y = wsp.y;
        collider->shape.mesh.wsps[j].z = wsp.z;
    }

    collider->shape.dirty = 0;

    // Update bounding sphere radius
    if (!collider->shape.scale_dirty)
    {
        return;
    }

    BoundingSphere* bs = &collider->shape.bs;
    V3 centre = bs->centre;

    // Calculate the new radius of the mi's bounding sphere.
    float radius_squared = -1;

    for (int j = 0; j < mb->num_positions; ++j)
    {
        V3 v = collider->shape.mesh.wsps[j];
        V3 between = v3_sub_v3(v, centre);

        radius_squared = max(size_squared(between), radius_squared);
    }

    bs->radius = sqrtf(radius_squared);

    collider->shape.scale_dirty = 0;
}

static void update_colliders(Physics* physics, Scene* scene)
{
    // TODO: Eventually we might want specific collision meshes to simplify it?
    // TODO: E.g. we don't need the higher detail meshes to collide with, but not needed for now.
    ViewIter it = ECS_view_iter(physics->ecs, physics->colliders_view);
    while (ECS_view_iter_next(&it))
    {
        const Transform* transforms = ECS_get_column(it, COMPONENT_Transform);
        const MeshInstance* mis = ECS_get_column(it, COMPONENT_MeshInstance);
        Collider* colliders = ECS_get_column(it, COMPONENT_Collider);

        // Update world space positions of entity, update centre of bounding sphere.
        for (int i = 0; i < it.num_entities; ++i)
        {
            Collider* collider = &colliders[i];

            // TODO: Where can we actually set this?????? after appling forces???
            //if (!collider->shape.dirty)
            //{
            //    continue;
            //}

            switch (collider->shape.type)
            {
            case COLLISION_SHAPE_MESH:
            {
                const Transform transform = transforms[i];
                const MeshInstance* mi = &mis[i];
                const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
                
                collider->shape.mesh.mb = mb;
                update_collision_mesh_bounding_sphere(collider, mi, transform);
                break;
            }
            case COLLISION_SHAPE_ELLIPSOID:
            {
                const Transform transform = transforms[i];
                collider->shape.bs.centre = transform.position;

                // TODO: not sure if this is really correct, could do for now.
                // TODO: I don't like this very much... e.g. if we used an ellipsoid collider on a 
                //       cube this wouldn't work, but then i suppose that's also wrong..... for now it's fine.
                collider->shape.bs.radius = max(collider->shape.ellipsoid.x, max(collider->shape.ellipsoid.y, collider->shape.ellipsoid.z));
                
                break;
            }
            default:
                log_error("Collision shape not defined!!\n");
                break;
            }
        }
    }
}

static void broad_phase(Physics* physics, Scene* scene, float dt)
{
    // Broad phase is purely bounding sphere tests
    // TODO: not ideal if the narrow phase is a sphere, but obv not an issue for now.

    Vector_clear(physics->frame.broad_phase_collisions);

    // TODO: Would be nicer to get the number of entities in a view (VIEW) or something instead?

    int num_entities = 0;
    {
        ViewIter it = ECS_view_iter(physics->ecs, physics->colliders_view);
        while (ECS_view_iter_next(&it))
        {
            num_entities += it.num_entities;
        }
    }

    // TODO: Idk what the calc is.
    Vector_reserve(physics->frame.broad_phase_collisions, num_entities * num_entities);

    /*
    TODO: A static mesh may not have physics data, therefore, the inner loop shouldn't use the
            collision view as it requires the physics data component.

    TODO: Refactor loop into separate ones. One for the moving vs moving and one for moving vs static

    this means physicsdata + collider vs physicsdata + collider and physicsdata + collider vs collider

    Note we can do this thing for the moving entities:
    for i in entities
        for j = i + 1 in entities

    */

    ViewIter it = ECS_view_iter(physics->ecs, physics->moving_colliders_view);
    while (ECS_view_iter_next(&it))
    {
        PhysicsData* physics_datas = ECS_get_column(it, COMPONENT_PhysicsData);
        Transform* transforms = ECS_get_column(it, COMPONENT_Transform);
        const MeshInstance* mis = ECS_get_column(it, COMPONENT_MeshInstance);
        const Collider* colliders = ECS_get_column(it, COMPONENT_Collider);

        for (int i = 0; i < it.num_entities; ++i)
        {
            PhysicsData* physics_data = &physics_datas[i];

            // TODO: Because of testing only past the current moving entity,
            //       we cannot have this check as it means a static mesh
            //       wouldn't 'collide' with a moving one...
            /*
            // If entity not moving, it cannot have collided with something.
            if (physics_data->velocity.x == 0 &&
                physics_data->velocity.y == 0 &&
                physics_data->velocity.z == 0)
            {
                continue;
            }*/

            Transform* transform = &transforms[i];
            const MeshInstance* mi = &mis[i];
            const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
            const Collider* collider = &colliders[i];

            // Iterate from the current entity, will this work????
            // TODO: This might not work if we don't resolve collisions
            //       properly between two entities!!!!!!! Otherwise only
            //       one will be updated!!!!!!
            ViewIter it_from_it0 = it; // TODO: Rename this....

            // Check each remaining moving entity
            do
            {
                // TODO: Detect.
                const MeshInstance* mis1 = ECS_get_column(it_from_it0, COMPONENT_MeshInstance);
                const Collider* colliders1 = ECS_get_column(it_from_it0, COMPONENT_Collider);
                PhysicsData* physics_datas1 = ECS_get_column(it_from_it0, COMPONENT_PhysicsData);

                // Start past current entity
                for (int j = i + 1; j < it_from_it0.num_entities; ++j)
                {
                    MeshInstance* mi1 = &mis1[j];

                    // TODO: Cannot collide with self for now, leave this as reminder
                    //       until the logic is validated.
                    // Don't collide with self.
                    //if (mi1 == mi) continue;

                    const Collider* collider1 = &colliders1[j];

                    // Account for entity1's velocity.
                    PhysicsData* physics_data1 = &physics_datas1[j];
                    const V3 rel_v = v3_sub_v3(physics_data->velocity, physics_data1->velocity);

                    BoundingSphere bs0 = collider->shape.bs;
                    const BoundingSphere bs1 = collider1->shape.bs;

                    // 'Sweep' sphere to account for velocities, otherwise we would miss
                    // collisions at low fps or high velocity. Instead of sweeping just
                    // scale and move bounding sphere.
                    v3_add_v3(bs0.centre, v3_mul_f(rel_v, 0.5f * dt));
                    bs0.radius += 0.5f * size(rel_v) * dt;

                    // Test for overlap.
                    const float dist = size_squared(v3_sub_v3(bs1.centre, bs0.centre));
                    const float n = (bs0.radius + bs1.radius) * (bs0.radius + bs1.radius);

                    if (dist <= n)
                    {
                        printf("potentially collidign WITH MOVING!\n");
                        // TODO:

                        // TODO: Write out potential collision, doesn't have to be perfect for now.
                        PotentialCollision pc = {
                            .collider_aid = it.aid,
                            .collider_offset = i,
                            .target_aid = it_from_it0.aid,
                            .target_offset = j
                        };

                        Vector_push_back(physics->frame.broad_phase_collisions, pc);
                    }
                }
            } while (ECS_view_iter_next(&it_from_it0));

            // Check with each static entity
            ViewIter sc_it = ECS_view_iter(physics->ecs, 
                physics->static_colliders_view);

            while (ECS_view_iter_next(&sc_it))
            {
                // TODO: Detect.
                const MeshInstance* mis1 = ECS_get_column(sc_it, COMPONENT_MeshInstance);
                const Collider* colliders1 = ECS_get_column(sc_it, COMPONENT_Collider);

                // Start past current entity
                for (int j = 0; j < sc_it.num_entities; ++j)
                {
                    MeshInstance* mi1 = &mis1[j];
                    const Collider* collider1 = &colliders1[j];

                    BoundingSphere bs0 = collider->shape.bs;
                    const BoundingSphere bs1 = collider1->shape.bs;

                    // 'Sweep' sphere to account for velocities, otherwise we would miss
                    // collisions at low fps or high velocity. Instead of sweeping just
                    // scale and move bounding sphere.
                    v3_add_v3(bs0.centre, v3_mul_f(physics_data->velocity, 0.5f * dt));
                    bs0.radius += 0.5f * size(physics_data->velocity) * dt;

                    // Test for overlap.
                    const float dist = size_squared(v3_sub_v3(bs1.centre, bs0.centre));
                    const float n = (bs0.radius + bs1.radius) * (bs0.radius + bs1.radius);

                    if (dist <= n)
                    {
                        printf("potentially collidign - WITH STATIC!\n");
                        // TODO:

                        // TODO: Write out potential collision, doesn't have to be perfect for now.
                        //PotentialCollision pc = {
                        //    .collider_aid = it_id,
                        //    .collider_offset = i,
                        //    .target_aid = it_id1,
                        //    .target_offset = j
                        //};
                        //physics_frame->broad_phase_collisions[physics_frame->num_potential_collisions++] = pc;
                    }
                }
            }
        }
    }
}

static void narrow_phase(Physics* physics, Scene* scene, float dt)
{    
    const int num_potential_collisions = Vector_size(physics->frame.broad_phase_collisions);
    for (int i = 0; i < num_potential_collisions; ++i)
    {
        PotentialCollision pc = physics->frame.broad_phase_collisions[i];

        /*
        Archetype* ca = &ecs->its[pc.collider_aid];

        Collider* collider = &(((Collider*)(ECS_get_column(ca, COMPONENT_Collider)))[pc.collider_offset]);
        Transform* collider_transform = &(((Transform*)(ECS_get_column(ca, COMPONENT_Transform)))[pc.collider_offset]);
        PhysicsData* collider_pd = &(((PhysicsData*)(ECS_get_column(ca, COMPONENT_PhysicsData)))[pc.collider_offset]);

        collider_pd->velocity = v3_uniform(0);

        Archetype* ta = &ecs->its[pc.target_aid];
        Collider* target_collider = &(((Collider*)(ECS_get_column(ta, COMPONENT_Collider)))[pc.target_offset]);
        Transform* target_transform = &(((Transform*)(ECS_get_column(ta, COMPONENT_Transform)))[pc.target_offset]);
        PhysicsData* target_pd = &(((PhysicsData*)(ECS_get_column(ta, COMPONENT_PhysicsData)))[pc.target_offset]);
        */
        //printf("%s\n", v3_to_str(target_pd->velocity));
        //v3_mul_eq_f(&target_pd->velocity, -1);
        //printf("%s\n\n", v3_to_str(target_pd->velocity));

    }
    /*
    for pair in collision pairs
        for face in mesh
            if thing collides with face
                add to resolve list
    
    */
}

static void resolve_collisions(Physics* physics, Scene* scene, float dt)
{

}

static void handle_collisions(Physics* physics, Scene* scene, float dt)
{
    /*
    TODO:
    Outline:

    - Broad phase gathers potential pairs of collisions
    - Narrow phase confirms collisions
    - Solver resolves contacts iteratively.

    */

    update_colliders(physics, scene);

    broad_phase(physics, scene, dt);

    narrow_phase(physics, scene, dt);



    /* TODO: What do I actually want to handle collisions for? Surely we're just thinking
             of players (potentially enemies) realistically for now. In the future it
             would be cool to simulate more, but that's outside the scope of what I need now.

             At some point we would want to handle collisions from some bounding area to
             faces or potentially bounding area to bounding area. This will depend on whatever
             the entity sets. - TODO: This could be a nice thing to setup now ahead of time.
    */

    // TODO: The broad phase should calculate pairs of collisions for the entire 
}

void PhysicsData_init(PhysicsData* data)
{
    memset(data, 0, sizeof(PhysicsData));

    // Mass of 0 will cause divide by zero error.
    data->mass = 1.f;
}

void Physics_tick(Physics* physics, Scene* scene, float dt)
{
    apply_forces(physics);

    handle_collisions(physics, scene, dt);

    apply_velocities(physics, dt);
}

void Collider_init(Collider* c)
{
    memset(c, 0, sizeof(Collider));
    c->shape.dirty = 1;
    c->shape.scale_dirty = 1;

    // TODO: Should we default to this??? Note we cannot set the ellipsoid
    //       as this would put random values into the other union members.
    c->shape.type = COLLISION_SHAPE_ELLIPSOID; 
}

void Collider_destroy(Collider* c)
{
    switch (c->shape.type)
    {
    case COLLISION_SHAPE_ELLIPSOID:
    {
        break;
    }
    case COLLISION_SHAPE_MESH:
    {
        Vector_destroy(c->shape.mesh.wsps);
        break;
    }
    default:
        break;
    }
}
