#include <stdio.h>

#include <engine/engine.h>

void create_map(Engine* engine)
{
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

    // Test point light
    point_lights_create(&scene->point_lights, rbs, (V3) { 0, 0, -3 }, (V3) { 1, 1, 1 }, 1.f);



    // TODO: Also this should just be done by a flag so at the start of the render,
    //       the buffers are resized. Or even we check each time.
    render_buffers_resize(&engine->renderer.buffers);
}

void engine_on_init(Engine* engine)
{
    create_map(engine);
}

void engine_on_update(Engine* engine, float dt) { }

void engine_on_keyup(Engine* engine, WPARAM wParam) { }


#include "light_mapping.h"

int main()
{
	Engine engine;
    Status status = engine_init(&engine, 1600, 900);
    if (STATUS_OK != status)
    {
        log_error("Failed to initialise engine\n");
    }

    do_light_mapping(&engine);

	engine_destroy(&engine);

	return 0;
}