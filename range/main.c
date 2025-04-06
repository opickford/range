#include <engine/engine.h>

#include <engine/globals.h>
#include <engine/canvas.h>

#include <engine/utils/common.h>

float* directions;

// TODO: Switch to C compiler.
// TODO: All of my code is C++ style C, why... Not good.

void create_map(Engine* engine)
{
    // TODO: Map like a fortnite 1v1 map, just a floor and scoreboard and maybe some obstacles.

    // Setup scene.
    Scene* scene = &engine->scenes[0];
    Status status = scene_init(scene);
    RenderBuffers* rbs = &engine->renderer.buffers;

    engine->current_scene_id = 0;
    ++engine->scenes_count;

    scene->ambient_light = (V3){ 0,0,0 };

    // Actually define the scene.
    mb_from_obj(&scene->models, rbs, "C:/Users/olive/source/repos/range/res/models/cube.obj");
    mi_create(&scene->models, rbs, 0, 1);

    // Create centre cube
    mi_set_transform(&scene->models, 0, (V3) { 0, 0, -5 }, (V3) { 0, 0, 0 }, (V3) { 1, 1, 1 });

    
    /*
    TODO: How can this mi load it's lightmap file?

    We need to tie some uv coordinates and the lightmap to it. So maybe I need to export a scene from the lightmap system?

    each mi now needs a new uv coordinate channel, so just duplicate the current uv data, but this time it won't be instanced.
    it will be 3 uvs per face. then we need another texture channel.

    so TODO:

    - make a new uv channel
    - make a new texture id etc

    how do we load the lightmap now?

    - all we're doing is loading a lightmap texture and uv coordinates. 
    - we can add these as args for mi_create, or do a mi_set_lightmap?
        - keep separate for now, so do mi_set_lightmap(uvs, lightmap_png_name)

    how to automate this:

    - the lightmap tool could export a txt file with the uv pairs on new lines as well as the lightmap png,
      with the mi id, then we can load it like mi_lightmap_id, mi_lightmap_uvs_id.txt!!!!!!!!!!!!!!!


    Good plan I think :D

    TODO: Actually, if the meshes are scaled, this will no longer work. 

    
    
    
    */

    // Test point light
    point_lights_create(&scene->point_lights, rbs, (V3) { 0, 0, -3 }, (V3) { 1, 1, 1 }, 100.f);

    // TODO: Also this should just be done by a flag so at the start of the render,
    //       the buffers are resized. Or even we check each time.
    render_buffers_resize(&engine->renderer.buffers);
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
        Scene* scene = &engine->scenes[engine->current_scene_id];

        V3 colour =
        {
            random_float(),
            random_float(),
            random_float()
        };

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
        engine->renderer.camera.position.x = scene->point_lights.world_space_positions[0];
        engine->renderer.camera.position.y = scene->point_lights.world_space_positions[1];
        engine->renderer.camera.position.z = scene->point_lights.world_space_positions[2];
        break;
    }
    case VK_F4:
    {
        g_debug_shadows = !g_debug_shadows;
        break;
    }
    case VK_F5:
    {
        Scene* scene = &engine->scenes[engine->current_scene_id];
        scene->point_lights.attributes[0] = 0.f;
        scene->point_lights.attributes[1] = 0.f;
        scene->point_lights.attributes[2] = 0.f;
        break;
    }
    }
}

int main()
{
	Engine engine;
	if (STATUS_OK == engine_init(&engine, 1600, 900))
	{
		engine_run(&engine);
	}
	
	engine_destroy(&engine);

	return 0;
}
