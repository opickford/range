#include <engine/engine.h>

#include <engine/globals.h>
#include <engine/canvas.h>

#include <engine/utils/common.h>

float* directions;

MeshBaseID sphere_base;

// TODO: Switch to C compiler.
// TODO: All of my code is C++ style C, why... Not good.

void create_map(Engine* engine)
{
    // TODO: Map like a fortnite 1v1 map, just a floor and scoreboard and maybe some obstacles.

    // TODO: Also for models, the scene could be global?
    Scene* scene = &engine->scenes[0];
    Status status = scene_init(scene);

    // TODO: Should really have a helper for this sort of thing?
    // TODO: Do we even want to support multiple scenes.. probably not.
    engine->current_scene_id = 0;
    ++engine->scenes_count;
    
    sphere_base = mesh_bases_add(&scene->mesh_bases);
    mesh_base_from_obj(&scene->mesh_bases.bases[sphere_base], "C:/Users/olive/source/repos/range/res/models/cube.obj");

    /*
    const int n = 10;

    float offset = 1;
    float x = -n * 0.5f * offset;
    float z = 0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            MeshInstanceID sphere_instance = mesh_instances_add(&scene->mesh_instances);
            scene_mesh_instance_set_base(scene, sphere_instance, sphere_base);
            scene->mesh_instances.instances[sphere_instance].position.x = x;
            scene->mesh_instances.instances[sphere_instance].position.z = z;
            z -= offset;
        }

        z = 0;
        x += offset;
    }
    */



    
    //resources_load_texture(&engine->resources, "C:/Users/olive/source/repos/range/res/textures/rickreal.bmp");

    /*
    mb_from_obj(&scene->models, &engine->renderer.buffers, "C:/Users/olive/source/repos/range/res/models/cube.obj");
    mi_create(&scene->models, &engine->renderer.buffers, 0, 1);
    mi_set_transform(&scene->models, 0, (V3) { 0, 0, -5 }, (V3) { 0, 0, 0 }, (V3) { 5, 5, 5 });
    */
    //scene->models.mis_texture_ids[0] = 0;

    //scene->ambient_light = (V3){ 0.1,0.1,0.1 };
    scene->ambient_light = (V3){ 1,1,1 };

    // TODO: Should the camera be part of the scene??
    engine->renderer.camera.position = (V3) { 0, 0, 10.f };
    
    // TODO: Also this should just be done by a flag so at the start of the render,
    //       the buffers are resized. Or even we check each time.
    //render_buffers_resize(&engine->renderer.buffers);
}

void engine_on_init(Engine* engine)
{
    g_draw_normals = 0;
    g_debug_shadows = 0;

    create_map(engine);

    /*
    // Create a scene
    Scene* scene = &engine->scenes[0];
    Status status = scene_init(scene);

    scene->ambient_light.x = 0.1f;
    scene->ambient_light.y = 0.1f;
    scene->ambient_light.z = 0.1f;

    if (STATUS_OK != status)
    {
        log_error("Failed to scene_init because of %s", status_to_str(status));
        return;
    }

    engine->current_scene_id = 0;
    ++engine->scenes_count;
    
    // Setup scene for shadow testing.
    // TODO: Could be nice to have a wrapper so I dont need to include the buffers param?
    load_model_base_from_obj(&scene->models, &engine->renderer.buffers, "C:/Users/olive/source/repos/range/range/res/models/cube.obj");
    load_model_base_from_obj(&scene->models, &engine->renderer.buffers, "C:/Users/olive/source/repos/range/range/res/models/suzanne.obj");
    
    V3 eulers = { 0, 0, 0 };

    create_model_instances(&scene->models, &engine->renderer.buffers, 0, 1);
    V3 plane_pos = { 0, 0, -4 };
    V3 plane_scale = { 5.f, 0.1f, 10.f };
    mi_set_transform(&scene->models, 0, plane_pos, eulers, plane_scale);

    if (0)
    {

        create_model_instances(&scene->models, &engine->renderer.buffers, 0, 2);
        V3 pos0 = { -1, 1, 3 };
        V3 pos1 = { 1, 1, 3 };
        
        V3 scale = { 0.5, 1, 0.5 };
        mi_set_transform(&scene->models, 1, pos0, eulers, scale);
        mi_set_transform(&scene->models, 2, pos1, eulers, scale);
    }
    else
    {
        
        create_model_instances(&scene->models, &engine->renderer.buffers, 1, 1);
        V3 pos = { 0, 1, 0 };

        V3 scale = { 1, 1, 1 };
        mi_set_transform(&scene->models, 1, pos, eulers, scale);  
        
    }

    V3 pl_pos0 = { 0, 2, 14 };
    V3 pl_col0 = { 1, 1, 1 };
    point_lights_create(&scene->point_lights, &engine->renderer.buffers, pl_pos0, pl_col0, 50.f);

    engine->renderer.camera.position.z = 20;

    // TODO: Maybe this is something that should be called after making changes to the models and lights?
    // TODO: But either way, we need to find a better way of doing this automatically because otherwise
    //       it will definitely cause some mistakes.
    render_buffers_resize(&engine->renderer.buffers);*/
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
        Scene* scene = &engine->scenes[engine->current_scene_id];

        V3 colour =
        {
            random_float(),
            random_float(),
            random_float()
        };

        const Camera* camera = &engine->renderer.camera;

        const V3 pos = v3_add_v3(camera->position, v3_mul_f(camera->direction, 10.f * (random_float() + 1)));

        MeshInstanceID sphere_instance = mesh_instances_add(&scene->mesh_instances);
        scene_mesh_instance_set_base(scene, sphere_instance, sphere_base);
        scene_mesh_instance_set_albedo(scene, sphere_instance, colour);
        mesh_instances_get(&scene->mesh_instances, sphere_instance)->position = pos;
        


        /*
        point_lights_create(&scene->point_lights, &engine->renderer.buffers, engine->renderer.camera.position, colour, 1);

        if (directions)
        {
             float* temp = realloc(directions, (size_t)scene->point_lights.count * 3 * sizeof(float));
             if (temp)
             {
                 directions = temp;
             }
        }
        else
        {
            directions = malloc((size_t)scene->point_lights.count * 3 * sizeof(float));
        }

        if (!directions)
        {
            printf("!directions.\n");
            return;
        }

        int i = (scene->point_lights.count - 1) * 3;
        directions[i] = engine->renderer.camera.direction.x;
        directions[i + 1] = engine->renderer.camera.direction.y;
        directions[i + 2] = engine->renderer.camera.direction.z;




        render_buffers_resize(&engine->renderer.buffers);
        */

        break;
    }
    case VK_F2:
    {
        g_draw_normals = !g_draw_normals;
        break;
    }
    case VK_F3:
    {
        Scene* scene = &engine->scenes[engine->current_scene_id];

        if (scene->mesh_instances.count > 0)
            mesh_instances_remove(&scene->mesh_instances, 
                scene->mesh_instances.id_to_index[scene->mesh_instances.count - 1]);

 
        /*
        Scene* scene = &engine->scenes[engine->current_scene_id];
        engine->renderer.camera.position.x = scene->point_lights.world_space_positions[0];
        engine->renderer.camera.position.y = scene->point_lights.world_space_positions[1];
        engine->renderer.camera.position.z = scene->point_lights.world_space_positions[2];
        */
        break;
    }
    case VK_F4:
    {
        g_debug_shadows = !g_debug_shadows;
        break;
    }
    case VK_F5:
    {
        /*
        Scene* scene = &engine->scenes[engine->current_scene_id];
        scene->point_lights.attributes[0] = 0.f;
        scene->point_lights.attributes[1] = 0.f;
        scene->point_lights.attributes[2] = 0.f;
        */
        break;
    }
    }
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
