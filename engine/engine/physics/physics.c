#include "physics.h"

#include "physics_frame.h"

#include "core/components.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"
#include "maths/vector_maths.h"

#include "utils/logger.h"

static void apply_forces(ECS* ecs, ViewID physics_view)
{
    // TODO: Air resistance

    // Disable gravity for now.
    static V3 acceleration = { 0, 0, 0 };
    //static V3 acceleration = { 0, -9.8f, 0 };

    ViewIter it = ECS_view_iter(ecs, physics_view);

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

static void apply_velocities(ECS* ecs, ViewID physics_view, float dt)
{
    ViewIter it = ECS_view_iter(ecs, physics_view);
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

static void update_colliders(ECS* ecs, Scene* scene, ViewID collision_view)
{
    // TODO: Eventually we might want specific collision meshes to simplify it?
    // TODO: E.g. we don't need the higher detail meshes to collide with, but not needed for now.

    ViewIter it = ECS_view_iter(ecs, collision_view);
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

static void broad_phase(PhysicsFrame* physics_frame, ECS* ecs, Scene* scene, ViewID collision_view, float dt)
{
    // Broad phase is purely bounding sphere tests
    // TODO: not ideal if the narrow phase is a sphere, but obv not an issue for now.


    physics_frame->num_potential_collisions = 0;

    // TODO: Would be nicer to get the number of entities in a view (VIEW) or something instead?

    int num_entities = 0;
    {
        ViewIter it = ECS_view_iter(ecs, collision_view);
        while (ECS_view_iter_next(&it))
        {
            num_entities += it.num_entities;
        }
    }

    // TODO: Idk what the calc is.
    Vector_reserve(physics_frame->broad_phase_collisions, num_entities * num_entities);

    /*
    TODO: A static mesh may not have physics data, therefore, the inner loop shouldn't use the
            collision view as it requires the physics data component.

    TODO: Refactor loop into separate ones. One for the moving vs moving and one for moving vs static

    this means physicsdata + collider vs physicsdata + collider and physicsdata + collider vs collider

    For this we are obviously going to need two views, however, the no physicsdata would also give
    entities that have physicsdata. For this we are going to need to update the ECS to allow a view
    to say it must NOT have a component. So we will have to edit ComponentsSignature to have a not mask
    ComponentsBitset so we set that it must NOT have the PhysicsData, this should solve our issue.

    Note we can do this thing for the moving entities:
    for i in entities
        for j = i + 1 in entities



    */

    ViewIter it = ECS_view_iter(ecs, collision_view);
    while (ECS_view_iter_next(&it))
    {
        PhysicsData* physics_datas = ECS_get_column(it, COMPONENT_PhysicsData);
        Transform* transforms = ECS_get_column(it, COMPONENT_Transform);
        const MeshInstance* mis = ECS_get_column(it, COMPONENT_MeshInstance);
        const Collider* colliders = ECS_get_column(it, COMPONENT_Collider);

        for (int i = 0; i < it.num_entities; ++i)
        {
            PhysicsData* physics_data = &physics_datas[i];

            // If entity not moving, it cannot have collided with something.
            if (physics_data->velocity.x == 0 &&
                physics_data->velocity.y == 0 &&
                physics_data->velocity.z == 0)
            {
                continue;
            }

            Transform* transform = &transforms[i];
            const MeshInstance* mi = &mis[i];
            const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
            const Collider* collider = &colliders[i];

            ViewIter it1 = ECS_view_iter(ecs, collision_view);
            while (ECS_view_iter_next(&it1))
            {
                const MeshInstance* mis1 = ECS_get_column(it1, COMPONENT_MeshInstance);
                const Collider* colliders1 = ECS_get_column(it1, COMPONENT_Collider);
                PhysicsData* physics_datas1 = ECS_get_column(it1, COMPONENT_PhysicsData);

                // TODO: The entity should realistically define it's narrow phase, as it's narrow phase may simply be 
                //       sphere etc, but not needed for now as we only want player to face collision.
                for (int j = 0; j < it1.num_entities; ++j)
                {
                    MeshInstance* mi1 = &mis1[j];

                    // Don't collide with self.
                    if (mi1 == mi) continue;

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
                        // TODO:
                        /*
                        // TODO: Write out potential collision, doesn't have to be perfect for now.
                        PotentialCollision pc = {
                            .collider_aid = it_id,
                            .collider_offset = i,
                            .target_aid = it_id1,
                            .target_offset = j
                        };
                        physics_frame->broad_phase_collisions[physics_frame->num_potential_collisions++] = pc;
                        */

                    }
                }
            }
        }
    }
}

static void narrow_phase(PhysicsFrame* physics_frame, ECS* ecs, Scene* scene, ViewID collision_view, float dt)
{
    
    for (int i = 0; i < physics_frame->num_potential_collisions; ++i)
    {
        PotentialCollision pc = physics_frame->broad_phase_collisions[i];

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

static void resolve_collisions(ECS* ecs, Scene* scene, ViewID collision_view, float dt)
{

}

static void handle_collisions(PhysicsFrame* physics_frame, ECS* ecs, Scene* scene, ViewID collision_view, float dt)
{
    /*
    TODO:
    Outline:

    - Broad phase gathers potential pairs of collisions
    - Narrow phase confirms collisions
    - Solver resolves contacts iteratively.

    */

    update_colliders(ecs, scene, collision_view);

    broad_phase(physics_frame, ecs, scene, collision_view, dt);

    narrow_phase(physics_frame, ecs, scene, collision_view, dt);



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

// TODO: fix whatever rubbish this is.
PhysicsFrame physics_frame;
int initialised = 0;
void Physics_tick(ECS* ecs, Scene* scene, ViewID physics_view, ViewID collision_view, float dt)
{
    if (!initialised)
    {
        PhysicsFrame_init(&physics_frame);
    }

    apply_forces(ecs, physics_view);

    handle_collisions(&physics_frame, ecs, scene, collision_view, dt);

    apply_velocities(ecs, physics_view, dt);
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
