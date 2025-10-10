#include "physics.h"

#include "core/components.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"
#include "maths/vector_maths.h"
#include "maths/plane.h"

#include "utils/logger.h"
#include "utils/common.h"

// TODO: Reorganise functions here.

// Note, this is specifically for solving the intersection time of a swept sphere
// with a triangle, this will not return a negative root as cannot have negative time.
// NOTE: This could be made more general if needed elsewhere
static uint8_t lowest_root(float a, float b, float c, float max_r, float* r)
{
    // Use quadratic formula to solve equation, return the lowest root below
    // max_r.

    // Calculate discriminant for number of roots.
    float d = b * b - 4.0f * a * c;

    // Negative means no roots.
    if (d < 0.0f) return 0;

    
    // Note, x0 == x1 if discriminant == 0 but not a necessary optimisation
    // makes code very messy without much/any gain.
    float sqrt_d = sqrtf(d);

    // x = (-b +/- sqrt(discriminant)) / 2a
    float x0 = (-b - sqrt_d) / (2.f * a);
    float x1 = (-b + sqrt_d) / (2.f * a);

    if (x0 > x1)
    {
        SWAP(float, x0, x1);
    }

    // Return the lowest postive root below given max.
    if (x0 > 0 && x0 < max_r)
    {
        *r = x0;
        return 1;
    }

    if (x1 > 0 && x1 < max_r)
    {
        *r = x1;
        return 1;
    }

    // No valid solutions
    return 0;
}

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

static void apply_forces(physics_t* physics, float dt)
{
    static float elapsed = 0.0f;


    // TODO: Air resistance applied as continuous force? note even when touching surface
    //       as still have to push through air duh.

    // Disable gravity for now.
    // TODO: Define in physics world.
    static const v3_t gravity = { 0, -9.8f, 0 };
    
    cecs_view_iter_t it = cecs_view_iter(physics->ecs, physics->physics_view);

    while (cecs_view_iter_next(&it))
    {
        physics_data_t* physics_datas = cecs_get_column(it, COMPONENT_PHYSICS_DATA);
        transform_t* transforms = cecs_get_column(it, COMPONENT_TRANSFORM);

        for (int i = 0; i < it.num_entities; ++i)
        {
            physics_data_t* physics_data = &physics_datas[i];
            transform_t* transform = &transforms[i];

            // Sum continuous acceleration.
            v3_t total_acceleration = { 0 };
            v3_add_eq_v3(&total_acceleration, gravity);

            // TODO: FRICTION SHOULD ACTUALLY BE APPLIED AT COLLISION RESPONSE, THEN
            //       THE SURFACE CAN DECIDE HOW STICKY IT IS!!!!!!!!!!!!!!!!!!!!!!!!    
            //       define in collider? 

            // TODO: Apply friction
            // TODO: Cache .
            float speed = v3_size(physics_data->velocity);
            if (speed > 0.f)
            {
                // TODO: Get some data from surfaces colliding with?
                // TODO: Because i dont have friction, objects slide more when laggy,
                //       friction should fix this.
            }

            // TODO: essentially i need to figure out what will work best for them game. maybe just remove this for now
            //       honestly. i don't know whether we should use a physically force based drag formula that will take in mass,
            //       and other parameters (requires more tuning per mi but might feel better), or just a dampening factor on the 
            //       velocity which would require configuring the damping factor, but that's it.

            // Air resistance/drag, simple mass sensitive 
            if (speed > 0.f)
            {
                /*
                // TODO: TESTING
                // Doesn't take into account mass because it's just affecting the velocity.
                float drag_rate = 1.5f;
                float damping = expf(-drag_rate * dt);
                v3_mul_eq_f(&physics_data->velocity, damping);
                */

                float drag_k = 0.35; // TODO: Parameter would need to be tuned. e.g. for a 1kg, 0.01m sphere, 0.35 is way too high.
                //       probs just in physics_data_t? or we could just dampen velocity, but then everytihng would
                //       drop at the same speed.

                float drag_mag = drag_k * speed * speed / physics_data->mass;

                // drag_accel = v3_normalised(v) * -drag_mag
                v3_t drag_accel = v3_mul_f(v3_mul_f(physics_data->velocity, 1.f / speed), -drag_mag);

                v3_add_eq_v3(&total_acceleration, drag_accel);

            }

            // Integrate continuous acceleration over dt.
            v3_add_eq_v3(&physics_data->velocity, v3_mul_f(total_acceleration, dt));

            // Apply continuous forces, over time, note, dt!
            //v3_add_eq_v3(&physics_data->velocity, v3_mul_f(continuous_force, dt / physics_data->mass));

            // Apply impulses.
            v3_add_eq_v3(&physics_data->velocity, v3_mul_f(physics_data->impulses, 1.f / physics_data->mass));

            elapsed += dt;



            
            


            // Clear impulses/instantaneous forces.
            physics_data->impulses = (v3_t){ 0.f, 0.f, 0.f };
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

        radius_squared = max(v3_size_sqrd(between), radius_squared);
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
                transform_t* transforms1 = cecs_get_column(it_from_it0, COMPONENT_TRANSFORM);

                // Start past current entity
                for (int j = i + 1; j < it_from_it0.num_entities; ++j)
                {
                    mesh_instance_t* mi1 = &mis1[j];

                    // TODO: Cannot collide with self for now, leave this as reminder
                    //       until the logic is validated.
                    // Don't collide with self.
                    //if (mi1 == mi) continue;

                    const collider_t* collider1 = &colliders1[j];
                    transform_t* transform1 = &transforms1[j];

                    // Account for entity1's velocity.
                    physics_data_t* physics_data1 = &physics_datas1[j];
                    const v3_t rel_v = v3_sub_v3(physics_data->velocity, physics_data1->velocity);

                    bounding_sphere_t bs0 = collider->shape.bs;
                    const bounding_sphere_t bs1 = collider1->shape.bs;

                    // 'Sweep' sphere to account for velocities, otherwise we would miss
                    // collisions at low fps or high velocity. Instead of sweeping just
                    // scale and move bounding sphere.
                    v3_add_v3(bs0.centre, v3_mul_f(rel_v, 0.5f * dt));
                    bs0.radius += 0.5f * v3_size(rel_v) * dt;

                    // Test for overlap.
                    const float dist = v3_size_sqrd(v3_sub_v3(bs1.centre, bs0.centre));
                    const float n = (bs0.radius + bs1.radius) * (bs0.radius + bs1.radius);

                    if (dist <= n)
                    {
                        //printf("potentially collidign WITH MOVING!\n");
                        // TODO:
                        
                        // TODO: Write out potential collision, doesn't have to be perfect for now.
                        potential_collision_t pc = {
                            .c0 = collider,
                            .mi0 = mi,
                            .pd0 = physics_data,
                            .t0 = transform,
                            .c1 = collider1,
                            .mi1 = mi1,
                            .pd1 = physics_data1,
                            .t1 = transform1
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
                const transform_t* transforms1 = cecs_get_column(sc_it, COMPONENT_TRANSFORM);

                // Start past current entity
                for (int j = 0; j < sc_it.num_entities; ++j)
                {
                    mesh_instance_t* mi1 = &mis1[j];
                    const collider_t* collider1 = &colliders1[j];
                    transform_t* transform1 = &transforms1[j];

                    bounding_sphere_t bs0 = collider->shape.bs;
                    const bounding_sphere_t bs1 = collider1->shape.bs;

                    // 'Sweep' sphere to account for velocities, otherwise we would miss
                    // collisions at low fps or high velocity. Instead of sweeping just
                    // scale and move bounding sphere.
                    v3_add_v3(bs0.centre, v3_mul_f(physics_data->velocity, 0.5f * dt));
                    bs0.radius += 0.5f * v3_size(physics_data->velocity) * dt;

                    // Test for overlap.
                    const float dist = v3_size_sqrd(v3_sub_v3(bs1.centre, bs0.centre));
                    const float n = (bs0.radius + bs1.radius) * (bs0.radius + bs1.radius);

                    if (dist <= n)
                    {
                        //printf("potentially collidign - WITH STATIC!\n");
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
                            .t0 = transform,
                            .c1 = collider1,
                            .mi1 = mi1,
                            .pd1 = 0, // static
                            .t1 = transform1
                        };

                        chds_vec_push_back(physics->frame.broad_phase_collisions, pc);

                        
                    }
                }
            }
        }
    }
}

// P1


static void test_point_with_ellipsoid(v3_t p, v3_t start, v3_t vel, float a, float *t, uint8_t* found_collision, v3_t* collision_point)
{
    float b = 2.f * dot(vel, v3_sub_v3(start, p));
    float c = v3_size_sqrd(v3_sub_v3(start, p)) - 1.f;

    float new_t = 0.f;
    if (lowest_root(a, b, c, *t, &new_t))
    {
        *t = new_t;
        *found_collision = 1;
        *collision_point = p;
    }
}

static void test_edge_with_ellipsoid(v3_t p0, v3_t p1, v3_t start, v3_t vel, float vel_size_sqrd, float* t, uint8_t* found_collision, v3_t* collision_point)
{
    // TODO: Refactor.
    v3_t edge = v3_sub_v3(p1, p0);
    v3_t baseToVertex = v3_sub_v3(p0, start);

    float edgeSquaredLength = v3_size_sqrd(edge);
    float edgeDotVelocity = dot(edge, vel);
    float edgeDotBaseToVertex = dot(edge, baseToVertex);

    // Calculate parameters for equation
    float a = edgeSquaredLength * -vel_size_sqrd +
        edgeDotVelocity * edgeDotVelocity;
    float b = edgeSquaredLength * (2 * dot(vel, baseToVertex)) -
        2.0 * edgeDotVelocity * edgeDotBaseToVertex;
    float c = edgeSquaredLength * (1 - v3_size_sqrd(baseToVertex)) +
        edgeDotBaseToVertex * edgeDotBaseToVertex;
    // Does the swept sphere collide against infinite edge?

    float new_t;
    if (lowest_root(a, b, c, *t, &new_t)) 
    {
        // Check if intersection is within line segment:
        float f = (edgeDotVelocity * new_t - edgeDotBaseToVertex) / edgeSquaredLength;
        if (f >= 0.0 && f <= 1.0) {
            // intersection took place within segment.
            *t = new_t;
            *found_collision = 1;
            *collision_point = v3_add_v3(p0, v3_mul_f(edge, f));
        }
    }
}

#include "core/globals.h"
#include <assert.h>
#include <float.h>
static void narrow_ellipsoid_vs_mi(physics_t* physics, scene_t* scene, potential_collision_t pc, float dt)
{

    // TODO: Note this is continuous detection, do we realy need this??? it feels kinda free here with this
    //       method? 

    // https://www.peroxide.dk/papers/collision/collision.pdf

    // TODO: We need to comment/enforce how an mi cannot be the thing that is colliding into something?
    //       we may want to allow this, but internally we will always be testing against the mi.
    //      
    //       Not really sure how to go about that but can solve it when needed? ^^^
    // TODO: Assert in broad phase to ensure this doesn't happen?

    // TODO: TEMP: TESTING ONLY VS STATIC
    // Seems to be working for moving as well, but have no repsonse for that one, also not taking their motion
    // into account with the velocity. TODO: Should do rel velocity.
    assert(pc.pd1 == 0);

    // c0 is ellipsoid collider !
    collider_t* ellipsoid_collider = pc.c0;
    physics_data_t* ellipsoid_pd = pc.pd0;
    transform_t* ellipsoid_transform = pc.t0;
    v3_t e_pos = ellipsoid_transform->position;

    v3_t ellipsoid = ellipsoid_collider->shape.ellipsoid;
    v3_t inv_ellipsoid = (v3_t){
        1.f / ellipsoid.x, 
        1.f / ellipsoid.y,
        1.f / ellipsoid.z
    };

    // Scale velocity by dt to calculate how much the entity will move this frame.
    // TODO: When we update the velocity should we take into account the remaining or not...
    v3_t vel = ellipsoid_pd->velocity;
    v3_t dir = v3_normalised(vel);
    v3_t e_vel = v3_mul_f(vel, dt);
    float vel_size = v3_size(e_vel);
    v3_t e_dir = v3_normalised(e_vel);

    v3_mul_eq_v3(&e_vel, inv_ellipsoid);

    //v3_t vel = pc.pd0->velocity;
    v3_t e_start_pos = e_pos;
    v3_mul_eq_v3(&e_start_pos, inv_ellipsoid);
    
    collider_t* mi_collider = pc.c1;
    //physics_data_t* mi_pd = pc.pd1;
    //mesh_instance_t* mi = pc.mi1; // TODO: Does this do anything for us here??

    mesh_base_t* mb = mi_collider->shape.mesh.mb;
    v3_t* wsps = mi_collider->shape.mesh.wsps;

    // TODO: We will have to deal with both velocities after getting it working one way first at least.

    // TODO: Rneame stuff.
    float earliest_t = 1.f;
    float nearest_dist = 0.f;
    v3_t nearest_collision_point = { 0 };
    uint8_t found_collision = 0;
    v3_t collision_face_normal = { 0 };

    // TODO: Slap this in a function that returns the collision?
    for (int i = 0; i < mb->num_faces * 3; i += 3)
    {
        // TODO: These are just wrong???
        v3_t p0 = wsps[mb->position_indices[i]];
        v3_t p1 = wsps[mb->position_indices[i + 1]];
        v3_t p2 = wsps[mb->position_indices[i + 2]];

        // Convert into ellipsoid space
        v3_mul_eq_v3(&p0, inv_ellipsoid);
        v3_mul_eq_v3(&p1, inv_ellipsoid);
        v3_mul_eq_v3(&p2, inv_ellipsoid);

        plane_t plane = plane_from_points(p0, p1, p2);

        // Backface culling, velocity and normal should face towards each other!
        // dot(A,B) = |A||B|cos(theta), note, we only care about sign.
        // If angle is > 90, they are facing each other so use: cos(>90) < 0
        if (dot(plane.normal, e_dir) >= 0)
        {
            continue;
        }

        float d = signed_distance(&plane, e_start_pos);
        float dot_normal_velocity = dot(plane.normal, e_vel);

        // Determine what time the sphere collides with the plane (if it does).
        uint8_t embedded_in_plane = 0;
        float t0 = 0.0f;
        float t1 = 1.0f;

        // If sphere moving parallel to the plane, e.g. 90deg between normal and velocity
        //if (dot_normal_velocity == 0.0f)
        if (fabsf(dot_normal_velocity) < FLT_EPSILON) 
        {
            if (fabs(d) >= 1.0f) 
            {
                // Sphere is not embedded in the plane so no collision.
                continue;
            }
            else 
            {
                // Sphere is embedded in plane, so it intersects in the whole range [0..1]
                embedded_in_plane = 1;
            }
        }
        else
        {
            // Calculate collision interval.
            t0 = (-1.0f - d) / dot_normal_velocity;
            t1 = (1.0f - d) / dot_normal_velocity;

            // TODO: This is giving a huge range for t0,t1? because dot_normal_velocity is so low?

            // Sort ascending.
            if (t0 > t1) SWAP(float, t0, t1);

            // Check that at least one result is within range:
            if (t0 > 1.0f || t1 < 0.0f) continue;

            // Clamp times of collision.
            if (t0 < 0.0f) t0 = 0.0f;
            if (t1 < 0.0f) t1 = 0.0f;
            if (t0 > 1.0f) t0 = 1.0f;
            if (t1 > 1.0f) t1 = 1.0f;
        }

        // At this point we know the sphere intersects the plane sometime between 
        // t0 and t1. Calculate the actual time.
        uint8_t found_new_collision = 0;
        v3_t collision_point = { 0 };
        float t = 1.0f;


        // Check if the collision is inside the triangle. This must occur at t0 as this is when the
        // sphere rests on the front side of the triangle plane. Can only happen if not embedded.
        // This will always happen before a collision with a point or an edge. - TODO: Why?

        if (!embedded_in_plane)
        {
            // As we are in ellipsoid space, we are working with a unit sphere (r = 1), 
            // so calculate the point on the sphere that would collide with the triangle,
            // then add velocity * t0 to get intersection.
            v3_t plane_intersection_p = v3_add_v3(v3_sub_v3(e_start_pos, plane.normal), v3_mul_f(e_vel, t0));
            
            if (point_in_triangle(plane_intersection_p, p0, p1, p2))
            {
                found_new_collision = 1;
                t = t0;
                collision_point = plane_intersection_p;
            }
        }
        
        // Haven't found collision so sweep sphere against points and edges of triangle.
        if (!found_new_collision)
        {
            float vel_size_sqrd = v3_size_sqrd(e_vel);

            // For each vertex or edge a quadratic equation must be solved, parameterise
            // equation to: a * t^2 + b * t + c = 0

            // Check against points
            float a = vel_size_sqrd;
            test_point_with_ellipsoid(p0, e_start_pos, e_vel, a, &t, &found_new_collision, &collision_point);
            test_point_with_ellipsoid(p1, e_start_pos, e_vel, a, &t, &found_new_collision, &collision_point);
            test_point_with_ellipsoid(p2, e_start_pos, e_vel, a, &t, &found_new_collision, &collision_point);
            
            // p1p0
            test_edge_with_ellipsoid(p0, p1, e_start_pos, e_vel, vel_size_sqrd, &t, &found_new_collision, &collision_point);
            test_edge_with_ellipsoid(p0, p2, e_start_pos, e_vel, vel_size_sqrd, &t, &found_new_collision, &collision_point);
            test_edge_with_ellipsoid(p2, p1, e_start_pos, e_vel, vel_size_sqrd, &t, &found_new_collision, &collision_point);
        }

        // TODO: Figure out nearest face.
        if (found_new_collision)
        {
            float dist = vel_size * t;

            if (!found_collision || dist < nearest_dist)
            {
                earliest_t = t;
                nearest_dist = dist;
                nearest_collision_point = collision_point;
                collision_face_normal = plane.normal;
                found_collision = 1;
            }
        }
    }

    // TODO: Should collisions be resolved elsewhere? idk, for now just stay like this 
    //       but this may change as we introduce more interactions.

    if (found_collision)
    {
        // TODO: Collision response, should set position and update velocity?
        //       how do we know when to slide vs bounce or whatever? that must depend on wall property or mass or something??
        
        float unit_scale = 1.f / 100.0f;
        float very_close_dist = 0.005f * unit_scale;

        v3_t V = v3_mul_f(v3_normalised(vel), nearest_dist - very_close_dist);
        ellipsoid_transform->position = v3_add_v3(ellipsoid_transform->position, V);


        v3_t normal = v3_mul_v3(collision_face_normal, ellipsoid);
        v3_normalise(&normal);
        printf("%s\n", v3_to_str(collision_face_normal));

        v3_t applied_vel = v3_mul_f(vel, earliest_t);

        // ellipsoid_pd->velocity = v3_sub_v3(vel, v3_mul_f(normal, dot(vel, normal)));
        ellipsoid_pd->velocity = v3_mul_f(normal, v3_size(applied_vel));// v3_sub_v3(vel, v3_mul_f(normal, dot(vel, normal)));
        
        //printf("%s\n", v3_to_str(ellipsoid_pd->velocity));
        //printf("%f - %s\n ", g_elapsed, v3_to_str(v3_mul_v3(collision_point, ellipsoid)));
        //ellipsoid_pd->velocity = (v3_t){ 0 };
    }
}

static void narrow_ellipsoid_vs_ellipsoid(physics_t* physics, scene_t* scene, potential_collision_t pc, float dt)
{
    // NOTE: Currently this means they must already collide as the broad phase is 
    //       sphere vs sphere.

    //log_error("narrow_ellipsoid_vs_ellipsoid not implemented!!\n");
    //assert(0);
     
}

static void narrow_phase(physics_t* physics, scene_t* scene, float dt)
{    
    const int num_potential_collisions = chds_vec_size(physics->frame.broad_phase_collisions);
    for (int i = 0; i < num_potential_collisions; ++i)
    {
        potential_collision_t pc = physics->frame.broad_phase_collisions[i];

        // Sort shapes ascending so we only have to solve A vs B, never B vs A.
        if (pc.c0->shape.type > pc.c1->shape.type)
        {
            SWAP(collider_t*, pc.c0, pc.c1);
            SWAP(mesh_instance_t*, pc.mi0, pc.mi1);
            SWAP(physics_data_t*, pc.pd0, pc.pd1);
            SWAP(physics_data_t*, pc.t0, pc.t1);
        }

        // TODO: Could sort by shape value (int) then only have to have to solve A vs B never B vs A.
        switch (pc.c0->shape.type)
        {
        case COLLISION_SHAPE_ELLIPSOID:
        {
            switch (pc.c1->shape.type)
            {
            case COLLISION_SHAPE_ELLIPSOID: narrow_ellipsoid_vs_ellipsoid(physics, scene, pc, dt);  break;
            case COLLISION_SHAPE_MESH: narrow_ellipsoid_vs_mi(physics, scene, pc, dt); break;
            default:
            {
                break;
            }
            }

            break;
        }
        case COLLISION_SHAPE_MESH:
        {
            // TODO: a mesh should never be colliding into something right???? So this is kinda invalid???
            log_error("Mesh instance collinding is not supported! Note, we should not have got here.\n"); break;

            switch (pc.c1->shape.type)
            {
            case COLLISION_SHAPE_MESH: log_error("Mesh instance vs mesh instance is not supported, note we shouldn't have got here!!!.\n"); break;
            default: break;
            }

            break;
        }
        default:
        {
            log_error("Collider not supported!");
            break;
        }
        }
    }
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

    // TODO: Do we actually want to be running continuous detection? Should I not 
    //       just use discrete 1/60 e.g?  Roblox and unity use this apparently.
    /*
    We should just do physics at 60fps essentially. If something is really fast
    we will miss it but who cares.
    
    */

    update_colliders(physics, scene);

    broad_phase(physics, scene, dt);

    narrow_phase(physics, scene, dt);
}

void physics_data_init(physics_data_t* data)
{
    memset(data, 0, sizeof(physics_data_t));

    // Mass of 0 will cause divide by zero error.

    // TODO: What unit is this?
    //data->mass = 1.f; 
    data->mass = 100.f; 
}

void physics_tick(physics_t* physics, scene_t* scene, float dt)
{
    apply_forces(physics, dt);

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
