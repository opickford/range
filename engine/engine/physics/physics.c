#include "physics.h"

#include "core/components.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"
#include "maths/vector_maths.h"

#include "utils/logger.h"

// TODO: Reorganise functions here.

static void physics_setup_views(physics_t* physics)
{
    physics->physics_view = cecs_view(physics->ecs,
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_PHYSICS_DATA) | 
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_TRANSFORM),
        0);

    physics->colliders_view = cecs_view(physics->ecs,
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_TRANSFORM) | 
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_MESH_INSTANCE) | 
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_COLLIDER),
        0);

    // TODO: Currently these views a mesh_instance_t is this correct?

    // Note, may not actually be moving, just has physicsdata so COULD be moving.
    physics->moving_colliders_view = cecs_view(physics->ecs,
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_PHYSICS_DATA) | 
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_TRANSFORM) |
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_MESH_INSTANCE) | 
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_COLLIDER),
        0);

    // Static means no physics data.
    physics->static_colliders_view = cecs_view(physics->ecs,
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_TRANSFORM) | 
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_MESH_INSTANCE) | 
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_COLLIDER),
        CECS_COMPONENT_ID_TO_BITSET(COMPONENT_PHYSICS_DATA));
}

status_t physics_init(physics_t* physics, cecs_t* ecs)
{
    memset(physics, 0, sizeof(physics_t));
    physics->ecs = ecs;

    physics_setup_views(physics);

    physics_frame_init(&physics->frame);

    return STATUS_OK;
}

static void apply_forces(physics_t* physics)
{
    // TODO: Air resistance

    // Disable gravity for now.
    static v3_t acceleration = { 0, 0, 0 };
    //static v3_t acceleration = { 0, -9.8f, 0 };

    cecs_view_iter_t it = cecs_view_iter(physics->ecs, physics->physics_view);

    while (cecs_view_iter_next(&it))
    {
        physics_data_t* physics_datas = cecs_get_column(it, COMPONENT_PHYSICS_DATA);
        transform_t* transforms = cecs_get_column(it, COMPONENT_TRANSFORM);

        for (int i = 0; i < it.num_entities; ++i)
        {
            physics_data_t* physics_data = &physics_datas[i];
            transform_t* transform = &transforms[i];

            // F = MA
            v3_add_eq_v3(&physics_data->force, v3_mul_f(acceleration, physics_data->mass));

            // Apply force.
            v3_add_eq_v3(&physics_data->velocity, v3_mul_f(physics_data->force, 1.f / physics_data->mass));

            // TODO: Function for zeroing/filling?
            physics_data->force = (v3_t){ 0.f, 0.f, 0.f };
        }
    }
}

static void apply_velocities(physics_t* physics, float dt)
{
    cecs_view_iter_t it = cecs_view_iter(physics->ecs, physics->physics_view);
    while (cecs_view_iter_next(&it))
    {
        physics_data_t* physics_datas = cecs_get_column(it, COMPONENT_PHYSICS_DATA);
        transform_t* transforms = cecs_get_column(it, COMPONENT_TRANSFORM);

        for (int i = 0; i < it.num_entities; ++i)
        {
            physics_data_t* physics_data = &physics_datas[i];
            transform_t* transform = &transforms[i];

            // Update position with new velocity.
            v3_add_eq_v3(&transform->position, v3_mul_f(physics_data->velocity, dt));
        }
    }
}

static void update_collision_mesh_bounding_sphere(collider_t* collider, const mesh_instance_t* mi, const transform_t transform)
{
    const mesh_base_t* mb = collider->shape.mesh.mb;

    chds_vec_reserve(collider->shape.mesh.wsps, mb->num_positions);

    // Calculate model matrix.
    // TODO: rotation or direction or eulers?
    m4_t model_matrix;
    m4_model_matrix(transform.position, transform.rotation, transform.scale, model_matrix);

    // Update the centre of the bounding sphere.
    v4_t world_space_centre;
    m4_mul_v4(
        model_matrix,
        v3_to_v4(mb->centre, 1.f),
        &world_space_centre);

    // Write out the bounding sphere.
    collider->shape.bs.centre = v4_xyz(world_space_centre);

    // Convert object space positions to world space.
    for (int j = 0; j < mb->num_positions; ++j)
    {
        const v4_t osp = v3_to_v4(mb->object_space_positions[j], 1.f);
        v4_t wsp;
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

    bounding_sphere_t* bs = &collider->shape.bs;
    v3_t centre = bs->centre;

    // Calculate the new radius of the mi's bounding sphere.
    float radius_squared = -1;

    for (int j = 0; j < mb->num_positions; ++j)
    {
        v3_t v = collider->shape.mesh.wsps[j];
        v3_t between = v3_sub_v3(v, centre);

        radius_squared = max(size_squared(between), radius_squared);
    }

    bs->radius = sqrtf(radius_squared);

    collider->shape.scale_dirty = 0;
}

static void update_colliders(physics_t* physics, scene_t* scene)
{
    // TODO: Eventually we might want specific collision meshes to simplify it?
    // TODO: E.g. we don't need the higher detail meshes to collide with, but not needed for now.
    cecs_view_iter_t it = cecs_view_iter(physics->ecs, physics->colliders_view);
    while (cecs_view_iter_next(&it))
    {
        const transform_t* transforms = cecs_get_column(it, COMPONENT_TRANSFORM);
        const mesh_instance_t* mis = cecs_get_column(it, COMPONENT_MESH_INSTANCE);
        collider_t* colliders = cecs_get_column(it, COMPONENT_COLLIDER);

        // Update world space positions of entity, update centre of bounding sphere.
        for (int i = 0; i < it.num_entities; ++i)
        {
            collider_t* collider = &colliders[i];

            // TODO: Where can we actually set this?????? after appling forces???
            //if (!collider->shape.dirty)
            //{
            //    continue;
            //}

            switch (collider->shape.type)
            {
            case COLLISION_SHAPE_MESH:
            {
                const transform_t transform = transforms[i];
                const mesh_instance_t* mi = &mis[i];
                const mesh_base_t* mb = &scene->mesh_bases.bases[mi->mb_id];
                
                collider->shape.mesh.mb = mb;
                update_collision_mesh_bounding_sphere(collider, mi, transform);
                break;
            }
            case COLLISION_SHAPE_ELLIPSOID:
            {
                const transform_t transform = transforms[i];
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

static void broad_phase(physics_t* physics, scene_t* scene, float dt)
{
    // Broad phase is purely bounding sphere tests
    // TODO: not ideal if the narrow phase is a sphere, but obv not an issue for now.

    chds_vec_clear(physics->frame.broad_phase_collisions);

    // TODO: How can we get the number of entities in a nicer way?

    int num_entities = 0;
    {
        cecs_view_iter_t it = cecs_view_iter(physics->ecs, physics->colliders_view);
        while (cecs_view_iter_next(&it))
        {
            num_entities += it.num_entities;
        }
    }

    // TODO: Idk what the calc is.
    chds_vec_reserve(physics->frame.broad_phase_collisions, num_entities * num_entities);

    // TODO: COmment properly.
    /*
    TODO: A static mesh may not have physics data, therefore, the inner loop shouldn't use the
            collision view as it requires the physics data component.

    TODO: Refactor loop into separate ones. One for the moving vs moving and one for moving vs static

    this means physicsdata + collider vs physicsdata + collider and physicsdata + collider vs collider

    Note we can do this thing for the moving entities:
    for i in entities
        for j = i + 1 in entities

    */

    cecs_view_iter_t it = cecs_view_iter(physics->ecs, physics->moving_colliders_view);
    while (cecs_view_iter_next(&it))
    {
        physics_data_t* physics_datas = cecs_get_column(it, COMPONENT_PHYSICS_DATA);
        transform_t* transforms = cecs_get_column(it, COMPONENT_TRANSFORM);
        const mesh_instance_t* mis = cecs_get_column(it, COMPONENT_MESH_INSTANCE);
        const collider_t* colliders = cecs_get_column(it, COMPONENT_COLLIDER);

        for (int i = 0; i < it.num_entities; ++i)
        {
            physics_data_t* physics_data = &physics_datas[i];

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

            transform_t* transform = &transforms[i];
            const mesh_instance_t* mi = &mis[i];
            const mesh_base_t* mb = &scene->mesh_bases.bases[mi->mb_id];
            const collider_t* collider = &colliders[i];

            // Iterate from the current entity, will this work????
            // TODO: This might not work if we don't resolve collisions
            //       properly between two entities!!!!!!! Otherwise only
            //       one will be updated!!!!!!
            cecs_view_iter_t it_from_it0 = it; // TODO: Rename this....

            // Check each remaining moving entity
            do
            {
                // TODO: Detect.
                const mesh_instance_t* mis1 = cecs_get_column(it_from_it0, COMPONENT_MESH_INSTANCE);
                const collider_t* colliders1 = cecs_get_column(it_from_it0, COMPONENT_COLLIDER);
                physics_data_t* physics_datas1 = cecs_get_column(it_from_it0, COMPONENT_PHYSICS_DATA);

                // Start past current entity
                for (int j = i + 1; j < it_from_it0.num_entities; ++j)
                {
                    mesh_instance_t* mi1 = &mis1[j];

                    // TODO: Cannot collide with self for now, leave this as reminder
                    //       until the logic is validated.
                    // Don't collide with self.
                    //if (mi1 == mi) continue;

                    const collider_t* collider1 = &colliders1[j];

                    // Account for entity1's velocity.
                    physics_data_t* physics_data1 = &physics_datas1[j];
                    const v3_t rel_v = v3_sub_v3(physics_data->velocity, physics_data1->velocity);

                    bounding_sphere_t bs0 = collider->shape.bs;
                    const bounding_sphere_t bs1 = collider1->shape.bs;

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
                        potential_collision_t pc = {
                            .c0 = collider,
                            .mi0 = mi,
                            .pd0 = physics_data,
                            .c1 = collider1,
                            .mi1 = mi1,
                            .pd1 = physics_data1,
                        };

                        chds_vec_push_back(physics->frame.broad_phase_collisions, pc);
                    }
                }
            } while (cecs_view_iter_next(&it_from_it0));

            // Check with each static entity
            cecs_view_iter_t sc_it = cecs_view_iter(physics->ecs, 
                physics->static_colliders_view);

            while (cecs_view_iter_next(&sc_it))
            {
                // TODO: Detect.
                const mesh_instance_t* mis1 = cecs_get_column(sc_it, COMPONENT_MESH_INSTANCE);
                const collider_t* colliders1 = cecs_get_column(sc_it, COMPONENT_COLLIDER);

                // Start past current entity
                for (int j = 0; j < sc_it.num_entities; ++j)
                {
                    mesh_instance_t* mi1 = &mis1[j];
                    const collider_t* collider1 = &colliders1[j];

                    bounding_sphere_t bs0 = collider->shape.bs;
                    const bounding_sphere_t bs1 = collider1->shape.bs;

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
                        //potential_collision_t pc = {
                        //    .collider_aid = it_id,
                        //    .collider_offset = i,
                        //    .target_aid = it_id1,
                        //    .target_offset = j
                        //};
                        //physics_frame->broad_phase_collisions[physics_frame->num_potential_collisions++] = pc;

                        
                        potential_collision_t pc = {
                            .c0 = collider,
                            .mi0 = mi,
                            .pd0 = physics_data,
                            .c1 = collider1,
                            .mi1 = mi1,
                            .pd1 = 0, // static
                        };

                        chds_vec_push_back(physics->frame.broad_phase_collisions, pc);

                        
                    }
                }
            }
        }
    }
}

static void narrow_phase(physics_t* physics, scene_t* scene, float dt)
{    
    const int num_potential_collisions = chds_vec_size(physics->frame.broad_phase_collisions);
    for (int i = 0; i < num_potential_collisions; ++i)
    {
        potential_collision_t pc = physics->frame.broad_phase_collisions[i];
        
        // TODO: how to map collision inter

        switch (pc.c0->shape.type)
        {
        case COLLISION_SHAPE_MESH:
        {
            // TODO: Func?
            break;
        }
        case COLLISION_SHAPE_ELLIPSOID:
        {
            // TODO: Func?
            break;
        }
        default:
        {
            log_error("Collision interaction not defined!!");
            break;
        }
        }
    }
    /*
    for pair in collision pairs
        for face in mesh
            if thing collides with face
                add to resolve list
    
    */
}

static void resolve_collisions(physics_t* physics, scene_t* scene, float dt)
{

}

static void handle_collisions(physics_t* physics, scene_t* scene, float dt)
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

void physics_data_init(physics_data_t* data)
{
    memset(data, 0, sizeof(physics_data_t));

    // Mass of 0 will cause divide by zero error.
    data->mass = 1.f;
}

void physics_tick(physics_t* physics, scene_t* scene, float dt)
{
    apply_forces(physics);

    handle_collisions(physics, scene, dt);

    apply_velocities(physics, dt);
}

void collider_init(collider_t* c)
{
    memset(c, 0, sizeof(collider_t));
    c->shape.dirty = 1;
    c->shape.scale_dirty = 1;

    // TODO: Should we default to this??? Note we cannot set the ellipsoid
    //       as this would put random values into the other union members.
    c->shape.type = COLLISION_SHAPE_ELLIPSOID; 
}

void collider_destroy(collider_t* c)
{
    switch (c->shape.type)
    {
    case COLLISION_SHAPE_ELLIPSOID:
    {
        break;
    }
    case COLLISION_SHAPE_MESH:
    {
        chds_vec_destroy(c->shape.mesh.wsps);
        break;
    }
    default:
        break;
    }
}
