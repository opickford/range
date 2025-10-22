#include "collision.h"

#include "physics_frame.h"

#include "core/mesh_base.h"
#include "core/transform.h"
#include "core/components.h"

#include "maths/vector_maths.h"
#include "maths/matrix4.h"
#include "maths/bounding_sphere.h"

#include "utils/logger.h"

#include <chds/vec.h>

#include <cecs/ecs.h>

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
        for (uint32_t i = 0; i < it.num_entities; ++i)
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

static void broad_phase(physics_t* physics, scene_t* scene)
{
    // Broad phase is purely bounding sphere tests
    // TODO: not ideal if the narrow phase is a sphere, but obv not an issue for now.
    //       but should work this out at some point ^

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

        for (uint32_t i = 0; i < it.num_entities; ++i)
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
                for (uint32_t j = i + 1; j < it_from_it0.num_entities; ++j)
                {
                    const mesh_instance_t* mi1 = &mis1[j];

                    // TODO: Cannot collide with self for now, leave this as reminder
                    //       until the logic is validated.
                    // Don't collide with self.
                    //if (mi1 == mi) continue;

                    const collider_t* collider1 = &colliders1[j];
                    transform_t* transform1 = &transforms1[j];

                    // Account for entity1's velocity.
                    physics_data_t* physics_data1 = &physics_datas1[j];

                    bounding_sphere_t bs0 = collider->shape.bs;
                    const bounding_sphere_t bs1 = collider1->shape.bs;

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
                for (uint32_t j = 0; j < sc_it.num_entities; ++j)
                {
                    mesh_instance_t* mi1 = &mis1[j];
                    const collider_t* collider1 = &colliders1[j];
                    transform_t* transform1 = &transforms1[j];

                    bounding_sphere_t bs0 = collider->shape.bs;
                    const bounding_sphere_t bs1 = collider1->shape.bs;

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

static void resolve_single_collision(v3_t rel_vel, v3_t collision_normal, float penetration_depth, physics_data_t* a_pd, physics_data_t* b_pd, transform_t* a_t, transform_t* b_t,
    float e, float mu)
{
    // Resolve a collision between objects A and B.
    const float slop = 0; // TODO: Do we want slop? This is actually stopping us from separating fully right????
    const float vel_along_n = dot(rel_vel, collision_normal);

    // Avoid jitter by moving slightly away from collision point.
    //float penetration_depth = max(0.f, vel_along_n * (1.f - collision_time) * dt);
    penetration_depth = max(penetration_depth - slop, 0.f);

    // If inv mass = 0, the object will not gain velocity
    const float a_inv_mass = a_pd->mass > 0.f ? (1.f / a_pd->mass) : 0.f;
    const float b_inv_mass = b_pd->mass > 0.f ? (1.f / b_pd->mass) : 0.f;

    const float total_inv_mass = a_inv_mass + b_inv_mass;
    if (total_inv_mass <= 0.f) return;

    // Separate objects.
    v3_t correction = v3_mul_f(collision_normal, penetration_depth / total_inv_mass);

    v3_add_eq_v3(&a_t->position, v3_mul_f(correction, a_inv_mass));
    v3_sub_eq_v3(&b_t->position, v3_mul_f(correction, b_inv_mass));

    // Normal impulse, apply restitution.

    // Ignore if separating as the current velocity wouldn't cause them 
    // to re-collide, therefore, no normal force and friction.
    // TODO: I don't think we will actually ever have this situatio???
    if (vel_along_n > 0.f) return;

    float j = -(1.f + e) * vel_along_n;
    j /= total_inv_mass;

    v3_t normal_impulse = v3_mul_f(collision_normal, j);

    // Tangential impulse, apply Coulomb friction.
    v3_t tangent = v3_sub_v3(rel_vel, v3_mul_f(collision_normal, dot(collision_normal, rel_vel)));
    const float size = v3_size(tangent);
    if (size > 0.f) v3_mul_eq_f(&tangent, 1.f / size);

    // Calculate size of friction force.
    float jt = -dot(rel_vel, tangent) / total_inv_mass;

    // Clamp friction to realistic maximum directly proportional to the normal force.
    // Normal impulse pushes objects apart, tangiential impulse slides along
    // surface but never exceeds coeff of friction * normal impulse. (directly proportional)

    const float max_friction = mu * fabsf(j);
    if (jt < -max_friction) jt = -max_friction;
    else if (jt > max_friction) jt = max_friction;

    v3_t friction_impulse = v3_mul_f(tangent, jt);

    // Apply impulses.
    v3_t total_impulse = v3_add_v3(normal_impulse, friction_impulse);

    v3_add_eq_v3(&a_pd->velocity, v3_mul_f(total_impulse, a_inv_mass));
    v3_sub_eq_v3(&b_pd->velocity, v3_mul_f(total_impulse, b_inv_mass));

    /* TODO: Clamp low velocities??? Should this test not be sqrd?
    if (v3_size_sqrd(a_pd->velocity) < 0.0001f)
    {
        a_pd->velocity = v3_uniform(0.f);
    }
    if (v3_size_sqrd(b_pd->velocity) < 0.0001f)
    {
        b_pd->velocity = v3_uniform(0.f);
    }*/
}

static void unit_sphere_tri_edge_collision(v3_t p0, v3_t p1, v3_t centre, uint8_t* collided, v3_t* collision_point, float* penetration_depth)
{
    // Essentially performing ray triangle collision but between t0 and t1,
    // the start and ends of the line.

    /*

    How formula is derived:

    To find the closest point on the ray to the centre, we want to minimise
    the distance between the points.

    Define some variables

        C = sphere centre.

        L(t) = P + Dt, where P is the start of the ray, D is the direction and t
                       is the time.

        Q = L(t), a point on the ray.

        x = || Q - C ||, x is the distance between a point on the ray and C.

    Remove the size from the equation

        size = sqrt(dot(A, B))

        x^2 = dot(Q - C, Q - C)

    Expand Q

        x^2 = dot( (P + Dt - C), (P + Dt - C) )

        A = P - C, vector between start point and centre.

        x^2 = dot( (A + Dt), (A + Dt) )

    Use distributive property of dot product

        x^2 = dot(A, A) + 2(A, Dt) + dot(Dt, Dt)

    Factor out scalars, we now have a quadratic.

        x^2 = dot(A, A) + 2t * dot(A, D) + (t^2) * dot(D, D)

    Differentiate for t

        d/Dt = 2 * dot(A, D) + 2t * dot(D, D)

    Solve for minimum t by setting d/Dt = 0, this works because
    in the quadratic equation at^2 + bt + c, the a term 'dot(A, D)'
    is always positive, meaning we have a U shaped parabola, so
    from our start point, we are decreasing, so we get the min.

        0 = 2 * dot(A, D) + 2t * dot(D, D)

        2t * dot(D, D) = -2 * dot(A, D)

        2t = (-2 * dot(A, D)) / dot(D, D)
        t = -dot(A, D) / dot(D, D)

        t = -dot(P - C, D) / dot(D, D)
        t = dot(C - P, D) / dot(D, D)

    If 0 <= t <= 1, the point is on the line. Note, if the
    closest point was outside of this range, the sphere could
    be colliding with the vertex or face instead.

    If the distance between the point and centre is less than the
    radius (1.f as we're working with a unit sphere),
    then we have a collision.

    */

    v3_t p1p0 = v3_sub_v3(p1, p0);
    float t = dot(p1p0, v3_sub_v3(centre, p0)) / dot(p1p0, p1p0);

    if (t >= 0 && t <= 1)
    {
        v3_t tmp_collision_point = v3_add_v3(p0, v3_mul_f(p1p0, t));

        if (v3_size_sqrd(v3_sub_v3(tmp_collision_point, centre)) <= 1.f)
        {
            // Only return a collision if the distance is greater than the
            // inputted penetration depth.
            float tmp_penetration_depth = 1.f - v3_size_sqrd(v3_sub_v3(tmp_collision_point, centre));
            if (!*collided || tmp_penetration_depth > *penetration_depth)
            {
                *collided = 1;
                *penetration_depth = tmp_penetration_depth;
                *collision_point = tmp_collision_point;
            }
        }
    }
}

static void unit_sphere_tri_vertex_collision(v3_t p, v3_t centre, uint8_t* collided, v3_t* collision_point, float* penetration_depth)
{
    // If distance between vertex and sphere centre is less than radius, we are colliding.
    float d = v3_size_sqrd(v3_sub_v3(p, centre));

    // Choose the point of furthest collision, this should ensure we 
    // push the ellipsoid fully out of the triangle!
    if (d < 1.f)
    {
        // Only return a collision if the distance is greater than the
        // inputted penetration depth.
        float tmp_penetration_depth = 1.f - d;
        if (!*collided || (tmp_penetration_depth > *penetration_depth))
        {
            *collided = 1;
            *collision_point = p;
            *penetration_depth = 1.f - d;
        }
    }
}

static void narrow_ellipsoid_vs_mi(physics_t* physics, scene_t* scene, potential_collision_t pc)
{
    // Inspired from: https://www.peroxide.dk/papers/collision/collision.pdf


    // TODO: We need to comment/enforce how an mi cannot be the thing that is colliding into something?
    //       we may want to allow this, but internally we will always be testing against the mi.
    //      
    //       Not really sure how to go about that but can solve it when needed? ^^^
    // TODO: Assert in broad phase to ensure this doesn't happen?

    physics_data_t* mi_pd = pc.pd1;
    transform_t* mi_transform = pc.t1;

    // c0 is ellipsoid collider !
    collider_t* ellipsoid_collider = pc.c0;
    physics_data_t* ellipsoid_pd = pc.pd0;
    transform_t* ellipsoid_transform = pc.t0;
    v3_t e_pos = ellipsoid_transform->position;

    v3_t ellipsoid = ellipsoid_collider->shape.ellipsoid;
    v3_t inv_ellipsoid = v3_inv(ellipsoid);

    // TODO: Rename vars and tidy this all up.

    // Calculate relative velocity between objects.
    v3_t vel = v3_sub_v3(ellipsoid_pd->velocity, mi_pd->velocity);

    v3_t e_start_pos = e_pos;
    v3_mul_eq_v3(&e_start_pos, inv_ellipsoid);

    collider_t* mi_collider = pc.c1;

    const mesh_base_t* mb = mi_collider->shape.mesh.mb;
    v3_t* wsps = mi_collider->shape.mesh.wsps;

    // TODO: Rneame stuff.

    // Per mi face, test for collision and pick the furthest point of collision 
    // on the face from the radius, this should ensure we push the ellipsoid 
    // fully out of the triangle!
    for (int i = 0; i < mb->num_faces * 3; i += 3)
    {
        v3_t p0 = wsps[mb->position_indices[i]];
        v3_t p1 = wsps[mb->position_indices[i + 1]];
        v3_t p2 = wsps[mb->position_indices[i + 2]];

        // Convert into ellipsoid space
        v3_mul_eq_v3(&p0, inv_ellipsoid);
        v3_mul_eq_v3(&p1, inv_ellipsoid);
        v3_mul_eq_v3(&p2, inv_ellipsoid);

        plane_t tri_plane = plane_from_points(p0, p1, p2);

        // Backface culling, velocity and normal should face towards each other!
        // dot(A,B) = |A||B|cos(theta), note, we only care about sign.
        // If angle is < 90, they not facing each other so use: cos(>90) < 0
        if (dot(vel, tri_plane.normal) >= 0) continue;

        float D = signed_distance(&tri_plane, e_start_pos);
        float dist = fabsf(D);

        uint8_t collided = 0;
        v3_t collision_point = { 0 };

        float penetration_depth = 0;

        // Test unit sphere with triangle plane.
        if (dist <= 1.f)
        {
            v3_t plane_collision_point = v3_sub_v3(e_start_pos, v3_mul_f(tri_plane.normal, D));

            if (point_in_triangle(plane_collision_point, p0, p1, p2))
            {
                collided = 1;
                collision_point = plane_collision_point;
                penetration_depth = 1.f - dist;
            }
        }

        // Test triangle vertices.
        unit_sphere_tri_vertex_collision(p0, e_start_pos, &collided, &collision_point, &penetration_depth);
        unit_sphere_tri_vertex_collision(p1, e_start_pos, &collided, &collision_point, &penetration_depth);
        unit_sphere_tri_vertex_collision(p2, e_start_pos, &collided, &collision_point, &penetration_depth);

        // Test against tri edges.
        unit_sphere_tri_edge_collision(p0, p1, e_start_pos, &collided, &collision_point, &penetration_depth);
        unit_sphere_tri_edge_collision(p0, p2, e_start_pos, &collided, &collision_point, &penetration_depth);
        unit_sphere_tri_edge_collision(p2, p1, e_start_pos, &collided, &collision_point, &penetration_depth);

        if (collided)
        {
            v3_t actual_plane_collision_point = v3_mul_v3(collision_point, ellipsoid);

            v3_t collision_normal = v3_sub_v3(ellipsoid_transform->position, actual_plane_collision_point);
            v3_normalise(&collision_normal);

            v3_t n = v3_add_v3(collision_point, collision_normal);
            v3_mul_eq_v3(&n, ellipsoid);

            float dist = v3_size(v3_sub_v3(n, e_pos));

            collision_data_t cd = { 0 };
            cd.collision_normal = collision_normal;
            cd.hit = 1;
            cd.rel_vel = vel;
            cd.penetration_depth = penetration_depth;
            cd.pc = pc;
            chds_vec_push_back(physics->frame.collisions, cd);
        }
    }
}

static void narrow_sphere_vs_sphere(physics_t* physics, scene_t* scene, potential_collision_t pc)
{
    // NOTE: Currently this means they must already collide as the broad phase is 
    //       sphere vs sphere.

    // We are just treating ellipsoids as spheres here.

    physics_data_t* a_pd = pc.pd0;
    physics_data_t* b_pd = pc.pd1;

    v3_t a_pos = pc.t0->position;
    v3_t b_pos = pc.t1->position;

    v3_t a_ellipsoid = pc.c0->shape.ellipsoid;
    v3_t b_ellipsoid = pc.c1->shape.ellipsoid;

    // Approximate ellipsoid with sphere.
    float a_radius = max(a_ellipsoid.x, max(a_ellipsoid.y, a_ellipsoid.z));
    float b_radius = max(b_ellipsoid.x, max(b_ellipsoid.y, b_ellipsoid.z));

    v3_t rel_v = v3_sub_v3(a_pd->velocity, b_pd->velocity);
    v3_t rel_p = v3_sub_v3(a_pos, b_pos);

    v3_t n = v3_normalised(rel_p);

    v3_t a_deepest = v3_add_v3(a_pos, v3_mul_f(n, -a_radius));
    v3_t b_deepest = v3_add_v3(b_pos, v3_mul_f(n, b_radius));

    float penetration_depth = v3_size(v3_sub_v3(a_deepest, b_deepest));

    collision_data_t cd = { 0 };
    cd.rel_vel = rel_v;
    cd.penetration_depth = penetration_depth;
    cd.collision_normal = n;
    cd.hit = 1;
    cd.pc = pc;

    chds_vec_push_back(physics->frame.collisions, cd);
}

static void narrow_phase(physics_t* physics, scene_t* scene)
{
    chds_vec_clear(physics->frame.collisions);

    const int num_potential_collisions = (int)chds_vec_size(physics->frame.broad_phase_collisions);
    for (int i = 0; i < num_potential_collisions; ++i)
    {
        potential_collision_t pc = physics->frame.broad_phase_collisions[i];

        // Sort shapes ascending so we only have to solve A vs B, never B vs A.
        if (pc.c0->shape.type > pc.c1->shape.type)
        {
            SWAP(collider_t*, pc.c0, pc.c1);
            SWAP(mesh_instance_t*, pc.mi0, pc.mi1);
            SWAP(physics_data_t*, pc.pd0, pc.pd1);
            SWAP(transform_t*, pc.t0, pc.t1);
        }

        collision_data_t cd = { 0 };

        // TODO: Could sort by shape value (int) then only have to have to solve A vs B never B vs A.
        switch (pc.c0->shape.type)
        {
        case COLLISION_SHAPE_ELLIPSOID:
        {
            switch (pc.c1->shape.type)
            {
            case COLLISION_SHAPE_ELLIPSOID: narrow_sphere_vs_sphere(physics, scene, pc);  break;
            case COLLISION_SHAPE_MESH: narrow_ellipsoid_vs_mi(physics, scene, pc); break;
            default:
            {
                log_error("Collision between COLLISION_SHAPE_ELLIPSOID and %d not supported.\n", (int)pc.c1->shape.type);
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

static void resolve_collisions(physics_t* physics, scene_t* scene)
{
    // Currently just iterating through the collisions, seems to handle everything well 
    // enough for now, can refactor this to solve collisions together in the future.
    const int num_collisions = (int)chds_vec_size(physics->frame.collisions);

    for (int i = 0; i < num_collisions; ++i)
    {
        collision_data_t cd = physics->frame.collisions[i];

        // Combine coefficients by taking averages.
        float e = max(0.f, (cd.pc.c0->restiution_coeff + cd.pc.c1->restiution_coeff) / 2.f);
        float mu = max(0.f, (cd.pc.c0->friction_coeff + cd.pc.c1->friction_coeff) / 2.f);

        resolve_single_collision(cd.rel_vel, cd.collision_normal, cd.penetration_depth, cd.pc.pd0, cd.pc.pd1, cd.pc.t0, cd.pc.t1, e, mu);
    }
}

uint8_t handle_collisions(physics_t* physics, scene_t* scene)
{
    /*
    Outline:

    - Broad phase gathers potential pairs of collisions
    - Narrow phase confirms collisions
    - Solver resolves contacts iteratively.

    */

    update_colliders(physics, scene);

    // Detect collisions
    {
        broad_phase(physics, scene);
        narrow_phase(physics, scene);
    }

    resolve_collisions(physics, scene);

    // Return true if we processed more collisions, signals to keep
    // iterating.
    return (uint8_t)((int)chds_vec_size(physics->frame.collisions) > 0.f);
}

void collider_init(collider_t* c)
{
    memset(c, 0, sizeof(collider_t));
    c->shape.dirty = 1;
    c->shape.scale_dirty = 1;

    // TOOD: Experiment for good defaults.
    c->friction_coeff = 0.25f;
    c->restiution_coeff = 0.5f;

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
