#include <engine/core/engine.h>

#include <engine/core/globals.h>
#include <engine/core/canvas.h>
#include <engine/core/transform.h>

#include <engine/utils/common.h>

#include <engine/maths/vector3.h>
#include <engine/common/status.h>

float* directions;

mesh_base_id_t sphere_base;
mesh_base_id_t cube_base;
mesh_base_id_t map_base;
mesh_base_id_t monkey_base;
mesh_base_id_t bowl_base;
mesh_base_id_t terrain_base;

cecs_entity_id_t map_entity;
cecs_entity_id_t monkey_entity;
cecs_entity_id_t player_entity;

uint8_t third_person = 1;

static void player_controller(engine_t* engine, float dt)
{

    if (engine->input_mode != INPUT_MODE_GAME) return;

    /*
    TODO: How can this be customisable?
    
    We really want this player controller to be provided by the engine.
    - Where/how can we provide it?

    Note, using this we see the choppy movement due to less physics updates,
    therefore, we need to lerp between positions for all entities.


    TODO:
    - Could accept an entity id.
    - Would this control the mouse as well?? 
    - This shouldn't be based off a camera direction right? more a player direction?
    
    */

    const static v3_t up = { 0, 1.f, 0 };
    const v3_t right = v3_normalised(cross(engine->renderer.camera.direction, up));

    physics_data_t* pd = cecs_get_component(engine->ecs, player_entity, COMPONENT_PHYSICS_DATA);
    transform_t* t = cecs_get_component(engine->ecs, player_entity, COMPONENT_TRANSFORM);

    // By multiplying by dt we're essentially converting the force into an impulse.
    const float speed = 20.f * dt * pd->mass;

    if (CSRGE_KEYDOWN(engine->window.keys['W']))
    {    
        const v3_t forward = v3_mul_f(engine->renderer.camera.direction, speed);
        v3_add_eq_v3(&pd->impulses, forward);
    }
    if (CSRGE_KEYDOWN(engine->window.keys['S']))
    {
        const v3_t forward = v3_mul_f(engine->renderer.camera.direction, speed);
        v3_sub_eq_v3(&pd->impulses, forward);
    }
    if (CSRGE_KEYDOWN(engine->window.keys['A']))
    {
        v3_sub_eq_v3(&pd->impulses, v3_mul_f(right, speed));
    }
    if (CSRGE_KEYDOWN(engine->window.keys['D']))
    {   
        v3_add_eq_v3(&pd->impulses, v3_mul_f(right, speed));
    }
    if (CSRGE_KEYDOWN(engine->window.keys[' ']))
    {

        //collider_t* c = cecs_get_component(engine->ecs, player_entity, COMPONENT_COLLIDER);

        
        // TODO: Only if colliding with something???
        const float jump_height = pd->mass;
        v3_add_eq_v3(&pd->impulses, v3_mul_f(up, jump_height));
        
        
    }

    if (third_person)
    {
        const static float cam_dist = 4.f;
        const static float lateral_offset = 2.f;
        const static float vertical_offset = 2.f;

        v3_t pos = v3_sub_v3(t->position, v3_mul_f(engine->renderer.camera.direction, cam_dist));
        v3_add_eq_v3(&pos, v3_mul_f(right, lateral_offset));
        v3_add_eq_v3(&pos, v3_mul_f(up, vertical_offset));
        engine->renderer.camera.position = pos;
    }
    else
    {
        engine->renderer.camera.position = t->position;
    }
    
}

void create_map(engine_t* engine)
{
    resources_load_texture(&engine->resources, "C:/Users/olive/source/repos/csrge/res/textures/landscape.bmp");
    resources_load_texture(&engine->resources, "C:/Users/olive/source/repos/csrge/res/textures/fortnite_peter.bmp");
    
    // TODO: Map like a csgo 1v1 map, just a floor and scoreboard and maybe some obstacles.
    
    // TODO: Also for models, the scene could be global? I think MBs should be 
    //       global?
    scene_t* scene = &engine->scene;
    status_t status = scene_init(scene);
    
    // Load some mesh bases.
    sphere_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[sphere_base], "C:/Users/olive/source/repos/csrge/res/models/sphere.obj");

    cube_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[cube_base], "C:/Users/olive/source/repos/csrge/res/models/cube.obj");

    map_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[map_base], "C:/Users/olive/source/repos/csrge/res/models/physics_test_map.obj");

    monkey_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[monkey_base], "C:/Users/olive/source/repos/csrge/res/models/suzanne.obj");

    bowl_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[bowl_base], "C:/Users/olive/source/repos/csrge/res/models/bowl.obj");

    terrain_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[terrain_base], "C:/Users/olive/source/repos/csrge/res/models/terrain.obj");

    // Create an entity
    {
        cecs_entity_id_t cube_entity = cecs_create_entity(engine->ecs);
        map_entity = cube_entity;

        // Add a mesh_instance_t component.
        mesh_instance_t* mi = cecs_add_component(engine->ecs, cube_entity, COMPONENT_MESH_INSTANCE);
        mesh_instance_init(mi, &scene->mesh_bases.bases[terrain_base]);

        mi->texture_id = 0;

        transform_t* transform = cecs_add_component(engine->ecs, cube_entity, COMPONENT_TRANSFORM);
        transform_init(transform);
        transform->scale = v3_uniform(10);


        // TODO: currently testing with static.
        //physics_data_t* physics_data = cecs_add_component(engine->ecs, cube_entity, COMPONENT_PHYSICS_DATA);
        //physics_data_init(physics_data);
        //physics_data->force = (v3_t){ 0,0,1 };

        collider_t* collider = cecs_add_component(engine->ecs, cube_entity, COMPONENT_COLLIDER);
        collider_init(collider);

        collider->shape.type = COLLISION_SHAPE_MESH;

        //physics_data_t* pd = cecs_add_component(engine->ecs, cube_entity, COMPONENT_PHYSICS_DATA);
        //physics_data_init(pd);
        //pd->mass = 0.f; // TODO: TEMP: Isn't moved by other things?
        //pd->floating = 1;
    }

    // Create player
    {
        player_entity = cecs_create_entity(engine->ecs);

        // Add a mesh_instance_t component.
        mesh_instance_t* mi = cecs_add_component(engine->ecs, player_entity, COMPONENT_MESH_INSTANCE);
        mesh_instance_init(mi, &scene->mesh_bases.bases[sphere_base]);

        mi->texture_id = 1;

        transform_t* transform = cecs_add_component(engine->ecs, player_entity, COMPONENT_TRANSFORM);
        transform_init(transform);
        transform->scale = (v3_t){ 1,2,1 };
        transform->position = (v3_t){ 0,3,0 };

        collider_t* collider = cecs_add_component(engine->ecs, player_entity, COMPONENT_COLLIDER);
        collider_init(collider);

        collider->shape.type = COLLISION_SHAPE_ELLIPSOID;
        collider->shape.ellipsoid = (v3_t){ 1,2,1 };

        physics_data_t* pd = cecs_add_component(engine->ecs, player_entity, COMPONENT_PHYSICS_DATA);
        physics_data_init(pd);
        pd->mass = 50.f;
    }
    
    // MONKEY
    /*
    {
        cecs_entity_id_t cube_entity = cecs_create_entity(engine->ecs);
        monkey_entity = cube_entity;

        // Add a mesh_instance_t component.
        mesh_instance_t* mi = cecs_add_component(engine->ecs, cube_entity, COMPONENT_MESH_INSTANCE);
        mesh_instance_init(mi, &scene->mesh_bases.bases[monkey_base]);

        mi->texture_id = 0;

        transform_t* transform = cecs_add_component(engine->ecs, cube_entity, COMPONENT_TRANSFORM);
        transform_init(transform);
        transform->scale = v3_uniform(1);

        // TODO: currently testing with static.
        //physics_data_t* physics_data = cecs_add_component(engine->ecs, cube_entity, COMPONENT_PHYSICS_DATA);
        //physics_data_init(physics_data);
        //physics_data->force = (v3_t){ 0,0,1 };

        collider_t* collider = cecs_add_component(engine->ecs, cube_entity, COMPONENT_COLLIDER);
        collider_init(collider);

        collider->shape.type = COLLISION_SHAPE_MESH;
        //collider->shape.ellipsoid = v3_uniform(1.f);

        physics_data_t* pd = cecs_add_component(engine->ecs, cube_entity, COMPONENT_PHYSICS_DATA);
        physics_data_init(pd);
    }*/
  

    /*
    // TODO: TEMP: Currently setting the spawned cube to have an ellipsoid collider, but this is just for the 
    //             broad phase which is using the sphere.
    // Normally to calculate bounding sphere of square we would half the scale, however,
    // the input .obj goes from -1 to 1, so the length is 2.
    v3_t half_sqrd = v3_mul_v3(transform->scale, transform->scale);
    float radius = sqrtf(half_sqrd.x + half_sqrd.y + half_sqrd.z);
    collider->shape.ellipsoid = v3_uniform(radius);
    printf("%f\n", radius);
    */
    //scene->ambient_light = v3_uniform(1.f);
    scene->ambient_light = v3_uniform(0.1f);
    
    scene->bg_colour = 0x11111111;
    
    // TODO: Should the camera be part of the scene??
    engine->renderer.camera.position = (v3_t) { 0, 0, 10.f };
}

void engine_on_init(engine_t* engine)
{
    // TODO: Should really be init by the engine!!!
    g_draw_normals = 0;
    g_debug_shadows = 0;
    g_debug_velocities = 0;

    create_map(engine);
}

void engine_on_update(engine_t* engine, float dt)
{
    {
        physics_data_t* pd = cecs_get_component(engine->ecs, map_entity, COMPONENT_PHYSICS_DATA);
        //pd->velocity = (v3_t){ 0.f, 0.f, -1.f };
    }

    {
        physics_data_t* pd = cecs_get_component(engine->ecs, player_entity, COMPONENT_PHYSICS_DATA);
        transform_t* t = cecs_get_component(engine->ecs, player_entity, COMPONENT_TRANSFORM);

        v3_t dir = v3_normalised(v3_sub_v3(engine->renderer.camera.position, t->position));
        
        //v3_add_eq_v3(&pd->impulses, v3_mul_f(dir, 0.01));

    }

    player_controller(engine, dt);
}

void engine_on_keyup(engine_t* engine, WPARAM wParam)
{

    switch (wParam)
    {
    case VK_F1:
    {
        mesh_base_id_t mb_ids[2] = { cube_base, sphere_base };

        mesh_base_id_t mb_id = mb_ids[(int)(random_float() + 0.5f)];

        scene_t* scene = &engine->scene;

        v3_t colour =
        {
            random_float(),
            random_float(),
            random_float()
        };

        
        const camera_t* camera = &engine->renderer.camera;

        const v3_t pos = v3_add_v3(camera->position, v3_mul_f(camera->direction, 10.f * (random_float() + 1)));

        cecs_entity_id_t cube_entity = cecs_create_entity(engine->ecs);
        mesh_instance_t* mi = cecs_add_component(engine->ecs, cube_entity, COMPONENT_MESH_INSTANCE);
        mesh_instance_init(mi, &scene->mesh_bases.bases[sphere_base]);
        mesh_instance_set_albedo(mi, &scene->mesh_bases.bases[sphere_base], colour);

        transform_t* transform = cecs_add_component(engine->ecs, cube_entity, COMPONENT_TRANSFORM);
        transform_init(transform);
        transform->position = pos;

        break;
    }
    case VK_F2:
    {
        //g_draw_normals = !g_draw_normals;
        scene_t* scene = &engine->scene;
        
        break;
    }
    case VK_F3:
    {
        v3_t colour =
        {
            random_float(),
            random_float(),
            random_float()
        };
        
        // TODO: Should camera be an entity component or entity or leave it? Not sure.
        //       Fine for now.
        camera_t* camera = &engine->renderer.camera;
        
        const v3_t pos = v3_add_v3(camera->position, v3_mul_f(camera->direction, 10.f * (random_float() + 1)));

        cecs_entity_id_t e = cecs_create_entity(engine->ecs);
        point_light_t* pl = cecs_add_component(engine->ecs, e, COMPONENT_POINT_LIGHT);
        pl->position = pos;
        pl->colour = colour;
        pl->strength = 1.f;

        break;
    }
    case VK_F4:
    {
        mesh_instance_t* mi = cecs_get_component(engine->ecs, 0, COMPONENT_MESH_INSTANCE);
        mesh_instance_destroy(mi);
        cecs_remove_component(engine->ecs, 0, COMPONENT_MESH_INSTANCE);

        cecs_destroy_entity(engine->ecs, 0);

        break;
    }
    case VK_F5:
    {
        engine->renderer.camera.position = (v3_t){ -5, 0, 2 };
        engine->renderer.camera.direction = (v3_t){ 1,0,0 };

        // TODO: Function for this.
        engine->renderer.camera.pitch = asinf(engine->renderer.camera.direction.y);
        engine->renderer.camera.yaw = atan2f(engine->renderer.camera.direction.x, engine->renderer.camera.direction.z);

        break;
    }
    case VK_F6:
    {
        g_debug_velocities = !g_debug_velocities;
        break;
    }
    case VK_F7:
    {
        cecs_entity_id_t cube_entity = cecs_create_entity(engine->ecs);
        map_entity = cube_entity;

        // Add a mesh_instance_t component.
        mesh_instance_t* mi = cecs_add_component(engine->ecs, cube_entity, COMPONENT_MESH_INSTANCE);
        mesh_instance_init(mi, &engine->scene.mesh_bases.bases[sphere_base]);

        transform_t* transform = cecs_add_component(engine->ecs, cube_entity, COMPONENT_TRANSFORM);
        transform_init(transform);
        transform->scale = v3_uniform(3);
        transform->position = (v3_t){ 0, 10, 0 };

        collider_t* collider = cecs_add_component(engine->ecs, cube_entity, COMPONENT_COLLIDER);
        collider_init(collider);

        collider->shape.type = COLLISION_SHAPE_ELLIPSOID;
        collider->shape.ellipsoid = transform->scale;


        physics_data_t* pd = cecs_add_component(engine->ecs, cube_entity, COMPONENT_PHYSICS_DATA);
        physics_data_init(pd);
        pd->mass = 100.f;

        break;
    }
    case VK_F8:
    {
        for (int i = 0; i < 10; ++i)
        {
            cecs_entity_id_t cube_entity = cecs_create_entity(engine->ecs);
            map_entity = cube_entity;

            // Add a mesh_instance_t component.
            mesh_instance_t* mi = cecs_add_component(engine->ecs, cube_entity, COMPONENT_MESH_INSTANCE);
            mesh_instance_init(mi, &engine->scene.mesh_bases.bases[sphere_base]);

            v3_t colour =
            {
                random_float(),
                random_float(),
                random_float()
            };

            float size = random_float();

            mesh_instance_set_albedo(mi, &engine->scene.mesh_bases.bases[sphere_base], colour);

            transform_t* transform = cecs_add_component(engine->ecs, cube_entity, COMPONENT_TRANSFORM);
            transform_init(transform);
            transform->scale = v3_uniform(size);
            transform->position = (v3_t){ 0, transform->scale.x + i * (transform->scale.x * 2), 0};

            collider_t* collider = cecs_add_component(engine->ecs, cube_entity, COMPONENT_COLLIDER);
            collider_init(collider);

            collider->shape.type = COLLISION_SHAPE_ELLIPSOID;
            collider->shape.ellipsoid = transform->scale;

            physics_data_t* pd = cecs_add_component(engine->ecs, cube_entity, COMPONENT_PHYSICS_DATA);
            physics_data_init(pd);
            pd->mass = 1.f;
        }
        
        break;
    }
    case 'H':
    {
        third_person = !third_person;
        break;
    }
    
    }
}

void engine_on_lmbdown(engine_t* engine)
{
    scene_t* scene = &engine->scene;

    // Create an entity
    cecs_entity_id_t cube_entity = cecs_create_entity(engine->ecs);

    //mesh_base_t* mb = &scene->mesh_bases.bases[cube_base];
    mesh_base_t* mb = &scene->mesh_bases.bases[sphere_base];

    // Add a mesh_instance_t component.
    cecs_add_component(engine->ecs, cube_entity, COMPONENT_MESH_INSTANCE);
    mesh_instance_t* mi = cecs_get_component(engine->ecs, cube_entity,
        COMPONENT_MESH_INSTANCE);
    mesh_instance_init(mi, mb);


    v3_t colour =
    {
        random_float(),
        random_float(),
        random_float()
    };
    mesh_instance_set_albedo(mi, mb, colour);

    cecs_add_component(engine->ecs, cube_entity, COMPONENT_TRANSFORM);
    transform_t* transform = cecs_get_component(engine->ecs, cube_entity, COMPONENT_TRANSFORM);
    transform_init(transform);
    transform->position = v3_add_v3(engine->renderer.camera.position, v3_mul_f(engine->renderer.camera.direction, 3));

    //transform->scale = v3_uniform(0.1);
    //transform->scale = (v3_t){ 0.5f,2,1 };

    physics_data_t* physics_data = cecs_add_component(engine->ecs, cube_entity, COMPONENT_PHYSICS_DATA);
    physics_data_init(physics_data);

    physics_data->mass = 1.f;

    // TODO: Dt?
    float speed = physics_data->mass * 20.f;
    v3_add_eq_v3(&physics_data->impulses, v3_mul_f(engine->renderer.camera.direction, speed));

    // TODO: Must remember that the pointers go invalid quick, should specifiy this in cecs!!!!

    collider_t* collider = cecs_add_component(engine->ecs, cube_entity, COMPONENT_COLLIDER);
    collider_init(collider);
    collider->shape.ellipsoid = transform->scale;
}

int main()
{
	engine_t engine;
	if (STATUS_OK == engine_init(&engine, 800, 600))
	{
		engine_run(&engine);
	}
	
	engine_destroy(&engine);

	return 0;
}
