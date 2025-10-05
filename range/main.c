#include <engine/core/engine.h>

#include <engine/core/globals.h>
#include <engine/core/canvas.h>
#include <engine/core/transform.h>

#include <engine/utils/common.h>

#include <engine/maths/vector3.h>
#include <engine/common/status.h>

float* directions;

MeshBaseID sphere_base;
MeshBaseID cube_base;

void create_map(Engine* engine)
{
    resources_load_texture(&engine->resources, "C:/Users/olive/source/repos/range/res/textures/rickreal.bmp");
    
    // TODO: Map like a csgo 1v1 map, just a floor and scoreboard and maybe some obstacles.
    
    // TODO: Also for models, the scene could be global? I think MBs should be 
    //       global?
    Scene* scene = &engine->scene;
    Status status = scene_init(scene);
    
    // Load some mesh bases.
    sphere_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[sphere_base], "C:/Users/olive/source/repos/range/res/models/sphere.obj");

    cube_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[cube_base], "C:/Users/olive/source/repos/range/res/models/cube.obj");

    // Create an entity
    EntityID cube_entity = ECS_create_entity(engine->ecs);

    // Add a MeshInstance component.
    MeshInstance* mi = ECS_add_component(engine->ecs, cube_entity, COMPONENT_MeshInstance);
    MeshInstance_init(mi, &scene->mesh_bases.bases[cube_base]);

    mi->texture_id = 0;

    Transform* transform = ECS_add_component(engine->ecs, cube_entity, COMPONENT_Transform);
    Transform_init(transform);

    PhysicsData* physics_data = ECS_add_component(engine->ecs, cube_entity, COMPONENT_PhysicsData);
    PhysicsData_init(physics_data);
    //physics_data->force = (V3){ 0,0,1 };

    Collider* collider = ECS_add_component(engine->ecs, cube_entity, COMPONENT_Collider);
    Collider_init(collider);
    
    collider->shape.type = COLLISION_SHAPE_MESH;

    /*
    // TODO: TEMP: Currently setting the spawned cube to have an ellipsoid collider, but this is just for the 
    //             broad phase which is using the sphere.
    // Normally to calculate bounding sphere of square we would half the scale, however,
    // the input .obj goes from -1 to 1, so the length is 2.
    V3 half_sqrd = v3_mul_v3(transform->scale, transform->scale);
    float radius = sqrtf(half_sqrd.x + half_sqrd.y + half_sqrd.z);
    collider->shape.ellipsoid = v3_uniform(radius);
    printf("%f\n", radius);
    */
    scene->ambient_light = v3_uniform(1.f);
    
    scene->bg_colour = 0x11111111;
    
    // TODO: Should the camera be part of the scene??
    engine->renderer.camera.position = (V3) { 0, 0, 10.f };
}

void engine_on_init(Engine* engine)
{
    g_draw_normals = 0;
    g_debug_shadows = 0;

    create_map(engine);
}

void engine_on_update(Engine* engine, float dt)
{
    return;
}

void engine_on_keyup(Engine* engine, WPARAM wParam)
{

    switch (wParam)
    {
    case VK_F1:
    {
        MeshBaseID mb_ids[2] = { cube_base, sphere_base };

        MeshBaseID mb_id = mb_ids[(int)(random_float() + 0.5f)];


        Scene* scene = &engine->scene;

        V3 colour =
        {
            random_float(),
            random_float(),
            random_float()
        };

        
        const Camera* camera = &engine->renderer.camera;

        const V3 pos = v3_add_v3(camera->position, v3_mul_f(camera->direction, 10.f * (random_float() + 1)));

        EntityID cube_entity = ECS_create_entity(engine->ecs);
        MeshInstance* mi = ECS_add_component(engine->ecs, cube_entity, COMPONENT_MeshInstance);
        MeshInstance_init(mi, &scene->mesh_bases.bases[sphere_base]);
        MeshInstance_set_albedo(mi, &scene->mesh_bases.bases[sphere_base], colour);

        Transform* transform = ECS_add_component(engine->ecs, cube_entity, COMPONENT_Transform);
        Transform_init(transform);
        transform->position = pos;

        break;
    }
    case VK_F2:
    {
        //g_draw_normals = !g_draw_normals;
        Scene* scene = &engine->scene;
        
        break;
    }
    case VK_F3:
    {
        V3 colour =
        {
            random_float(),
            random_float(),
            random_float()
        };
        
        // TODO: Should camera be an entity component or entity or leave it? Not sure.
        //       Fine for now.
        Camera* camera = &engine->renderer.camera;
        
        const V3 pos = v3_add_v3(camera->position, v3_mul_f(camera->direction, 10.f * (random_float() + 1)));

        EntityID e = ECS_create_entity(engine->ecs);
        PointLight* pl = ECS_add_component(engine->ecs, e, COMPONENT_PointLight);
        pl->position = pos;
        pl->colour = colour;
        pl->strength = 1.f;

        break;
    }
    case VK_F4:
    {
        MeshInstance* mi = ECS_get_component(engine->ecs, 0, COMPONENT_MeshInstance);
        MeshInstance_destroy(mi);
        ECS_remove_component(engine->ecs, 0, COMPONENT_MeshInstance);

        ECS_destroy_entity(engine->ecs, 0);

        break;
    }
    }
}

void engine_on_lmbdown(Engine* engine)
{
    Scene* scene = &engine->scene;

    // Create an entity
    EntityID cube_entity = ECS_create_entity(engine->ecs);

    MeshBase* mb = &scene->mesh_bases.bases[cube_base];

    // Add a MeshInstance component.
    ECS_add_component(engine->ecs, cube_entity, COMPONENT_MeshInstance);
    MeshInstance* mi = ECS_get_component(engine->ecs, cube_entity,
        COMPONENT_MeshInstance);
    MeshInstance_init(mi, mb);


    V3 colour =
    {
        random_float(),
        random_float(),
        random_float()
    };
    MeshInstance_set_albedo(mi, mb, colour);

    ECS_add_component(engine->ecs, cube_entity, COMPONENT_Transform);
    Transform* transform = ECS_get_component(engine->ecs, cube_entity, COMPONENT_Transform);
    Transform_init(transform);
    transform->position = v3_add_v3(engine->renderer.camera.position, v3_mul_f(engine->renderer.camera.direction, 3));
    transform->rotation = engine->renderer.camera.direction;
    transform->scale = v3_uniform(0.1f);

    PhysicsData* physics_data = ECS_add_component(engine->ecs, cube_entity, COMPONENT_PhysicsData);
    PhysicsData_init(physics_data);

    // TODO: Dt?
    float speed = 20;
    physics_data->force = v3_mul_f(engine->renderer.camera.direction, speed);

    // TODO: Must remember that the pointers go invalid quick, should specifiy this in cecs!!!!

    Collider* collider = ECS_add_component(engine->ecs, cube_entity, COMPONENT_Collider);
    Collider_init(collider);
    collider->shape.ellipsoid = transform->scale;
}

int main()
{
	Engine engine;
	if (STATUS_OK == engine_init(&engine, 800, 600))
	{
		engine_run(&engine);
	}
	
	engine_destroy(&engine);

	return 0;
}
