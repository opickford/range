#include "physics.h"

#include "collision.h"

#include "core/components.h"

#include "maths/vector3.h"
#include "maths/matrix4.h"
#include "maths/vector_maths.h"
#include "maths/plane.h"

#include "utils/logger.h"
#include "utils/common.h"

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

    physics->max_collision_iters = 3;

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

        for (uint32_t i = 0; i < it.num_entities; ++i)
        {
            physics_data_t* physics_data = &physics_datas[i];

            if (physics_data->mass == 0.f) continue;

            // TODO: Is this correct? 
            // Ignore objects that shouldn't be affected by forces.
            // TODO: Instead of this maybe we could do have a floating flag.
            

            transform_t* transform = &transforms[i];

            // Sum continuous acceleration.
            v3_t total_acceleration = { 0 };

            if (!physics_data->floating)
            {
                v3_add_eq_v3(&total_acceleration, gravity);
            }

            // TODO: FRICTION SHOULD ACTUALLY BE APPLIED AT COLLISION RESPONSE, THEN
            //       THE SURFACE CAN DECIDE HOW STICKY IT IS!!!!!!!!!!!!!!!!!!!!!!!!    
            //       define in collider? 

            // TODO: Apply friction
            // TODO: Cache .
            float speed = v3_size(physics_data->velocity);

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

                float drag_k = 0.35f; // TODO: Parameter would need to be tuned. e.g. for a 1kg, 0.01m sphere, 0.35 is way too high.
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

        for (uint32_t i = 0; i < it.num_entities; ++i)
        {
            physics_data_t* physics_data = &physics_datas[i];
            transform_t* transform = &transforms[i];

            // Update position with new velocity.
            v3_add_eq_v3(&transform->position, v3_mul_f(physics_data->velocity, dt));
        }
    }
}

void physics_data_init(physics_data_t* data)
{
    memset(data, 0, sizeof(physics_data_t));

    // Mass of 0 means the object cannot be pushed.
    data->mass = 1.f; 
}

void physics_tick(physics_t* physics, scene_t* scene, float dt)
{
    // TODO: Comment how all this works.

    apply_forces(physics, dt);

    apply_velocities(physics, dt);

    for (uint8_t i = 0; i < physics->max_collision_iters; ++i)
    {
        if (!handle_collisions(physics, scene))
        {
            break;
        }

    }

}

