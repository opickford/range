#include "engine.h"

#include "globals.h"
#include "components.h"

#include "utils/logger.h"
#include "utils/timer.h"

#include "common/colour.h"

#include <Windows.h>

#include <stdio.h>
#include <stdlib.h>

// Internal helpers.
static void Engine_setup_ecs(Engine* engine)
{
    /*
    
    TODO: Honestly my systems at the moment are more just collections of 
    archetypes that match their bitset, the 'system' is not really right?
    it's kinda more of a view of the entities with certain components?
    
    */

    ECS* ecs = &engine->ecs;
    ECS_init(ecs);

    // Setup core components, MeshInstance etc.
    CoreComponents_init(ecs);

    // TODO: How does the render system also process lights.

    // Setup systems.
    {
        engine->render_system_id = ECS_register_system(ecs);
        System* system = &engine->ecs.systems[engine->render_system_id];
        system->components_bitset = COMPONENT_ID_TO_BITSET(COMPONENT_MeshInstance) | 
                                    COMPONENT_ID_TO_BITSET(COMPONENT_Transform);
    }

    {
        engine->lighting_system_id = ECS_register_system(ecs);
        System* system = &engine->ecs.systems[engine->lighting_system_id];
        system->components_bitset = COMPONENT_ID_TO_BITSET(COMPONENT_PointLight);

    }

    {
        engine->physics_system_id = ECS_register_system(ecs);
        System* system = &engine->ecs.systems[engine->physics_system_id];
        system->components_bitset = COMPONENT_ID_TO_BITSET(COMPONENT_PhysicsData) |
                                    COMPONENT_ID_TO_BITSET(COMPONENT_Transform);
    }

    {
        engine->collision_system_id = ECS_register_system(ecs);
        System* system = &engine->ecs.systems[engine->collision_system_id];
        system->components_bitset = COMPONENT_ID_TO_BITSET(COMPONENT_PhysicsData) |
                                    COMPONENT_ID_TO_BITSET(COMPONENT_Transform) |
                                    COMPONENT_ID_TO_BITSET(COMPONENT_MeshInstance);
    }
}



Status engine_init(Engine* engine, int window_width, int window_height)
{
    log_info("Initialising the engine.");
    memset(engine, 0, sizeof(Engine));

    // Set some default settings.
    engine->upscaling_factor = 1;
    engine->handle_input = 0;

    // Initialise the renderer.
    Status status = renderer_init(&engine->renderer, (int)(window_width / engine->upscaling_factor), (int)(window_height / engine->upscaling_factor));
    if (STATUS_OK != status)
    {
        log_error("Failed to renderer_init because of %s", Status_to_str(status));
        return status;
    }

    // Initialise the window.
    status = window_init(&engine->window, &engine->renderer.target.canvas, (void*)engine, window_width, window_height);
    if (STATUS_OK != status)
    {
        log_error("Failed to window_init because of %s", Status_to_str(status));
        return status;
    }

    // Set window callbacks.
    engine->window.on_resize = &engine_on_resize;
    engine->window.on_keyup = &engine_process_keyup;
    engine->window.on_lmbdown = &engine_process_lmbdown;

    // Initialise the UI.
    status = ui_init(&engine->ui, &engine->renderer.target.canvas);
    if (STATUS_OK != status)
    {
        log_error("Failed to ui_init because of %s", Status_to_str(status));
        return status;
    }

    // Initialise the resources.
    resources_init(&engine->resources);

    // Initialise a random seed. TODO: Might not want this for testing.
    srand((unsigned int)time(NULL));

    // Initialise the ECS.
    ECS_init(&engine->ecs);

    // Initialise components and systems.
    Engine_setup_ecs(engine);

    log_info("Fired engine_on_init event.");
    engine_on_init(engine);

    log_info("Engine successfully initialised.");

    return STATUS_OK;
}

void engine_run(Engine* engine)
{
    // TEMP
    LARGE_INTEGER frequency = { 0 };
    QueryPerformanceFrequency(&frequency);

    LARGE_INTEGER startTime = { 0 };
    LARGE_INTEGER endTime = { 0 };

    // TODO: Struct of perf data?
    int fps = 0;
    float dt = 0;
    int draw_ui_ms = 0;

    // Start the timers.
    QueryPerformanceCounter(&startTime);

    float dt_counter = 0;

    char fps_str[64] = "";
    char dir_str[64] = "";
    char pos_str[64] = "";

    
    char process_messages_str[64] = "";
    char handle_input_str[64] = "";
    char rt_clear_str[64] = "";
    char render_str[64] = "";
    char ui_draw_str[64] = "";
    char display_str[64] = "";
    char update_str[64] = "";
    char vertices_str[64] = "";
    

    // TODO: This can be a ui add text function
    int y = 10;
    const int TEXT_SCALE = 1;
    int h = engine->ui.font.char_height * TEXT_SCALE;
    
    engine->ui.text[engine->ui.text_count++] = text_create(fps_str, 10, engine->ui.text_count * h + 10, COLOUR_LIME, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(dir_str, 10, engine->ui.text_count * h + 10, COLOUR_RED, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(pos_str, 10, engine->ui.text_count * h + 10, COLOUR_RED, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(process_messages_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(handle_input_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(rt_clear_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(render_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(ui_draw_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(display_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(update_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);
    engine->ui.text[engine->ui.text_count++] = text_create(vertices_str, 10, engine->ui.text_count * h + 10, COLOUR_WHITE, TEXT_SCALE);

    engine->running = 1;
    while (engine->running)
    {
        Timer t = timer_start();

        // Process the application window messages.
        if (!window_process_messages())
        {
            // Break on WM_QUIT.
            break;
        }

        snprintf(process_messages_str, sizeof(process_messages_str), "ProcMsgs: %d", timer_get_elapsed(&t));
       
        M4 view_matrix;
        calculate_view_matrix(&engine->renderer.camera, view_matrix);

        // Apply physics
        Physics_tick(&engine->ecs, &engine->ecs.systems[engine->physics_system_id], &engine->ecs.systems[engine->collision_system_id], dt);

        // Clear the canvas.
        timer_restart(&t);
        render_target_clear(&engine->renderer.target, engine->scene.bg_colour);
        snprintf(rt_clear_str, sizeof(rt_clear_str), "RTClear: %d", timer_get_elapsed(&t));

        // Render scene.
        timer_restart(&t);
        
        // TODO: Render system..............

        // TODO: How will this work as a system?

            
        render(&engine->ecs, 
            &engine->ecs.systems[engine->render_system_id], 
            &engine->ecs.systems[engine->lighting_system_id], 
            &engine->renderer, 
            &engine->scene, 
            &engine->resources,
            view_matrix);
        
        snprintf(render_str, sizeof(render_str), "Render: %d", timer_get_elapsed(&t));
         
        // Handle any keyboard/mouse input.
        timer_restart(&t);
        engine_handle_input(engine, dt);
        snprintf(handle_input_str, sizeof(handle_input_str), "HandleInput: %d", timer_get_elapsed(&t));

        // Fire the engine update event.
        timer_restart(&t);
        engine_on_update(engine, dt);
        snprintf(update_str, sizeof(update_str), "UpdateEvent: %d", timer_get_elapsed(&t));

        // Draw ui elements.
        timer_restart(&t);
        ui_draw(&engine->ui, engine->upscaling_factor);
        snprintf(ui_draw_str, sizeof(render_str), "DrawUI: %d", draw_ui_ms);
        draw_ui_ms = timer_get_elapsed(&t); // Must be done a frame late.

        // Update the display.
        timer_restart(&t);
        window_display(&engine->window);
        snprintf(display_str, sizeof(display_str), "Display: %d", timer_get_elapsed(&t));

        // Calculate performance.
        QueryPerformanceCounter(&endTime);
        dt = (float)(endTime.QuadPart - startTime.QuadPart) / frequency.QuadPart;

        dt_counter += dt;

        startTime = endTime;

        fps = (int)(1.0f / dt);
        
        if (dt_counter > 2)
        {
            snprintf(fps_str, sizeof(fps_str), "FPS: %d", fps);
            dt_counter = 0;
        }

        snprintf(dir_str, sizeof(dir_str), "DIR: %.2f %.2f %.2f", engine->renderer.camera.direction.x, engine->renderer.camera.direction.y, engine->renderer.camera.direction.z);
        snprintf(pos_str, sizeof(pos_str), "POS: %.2f %.2f %.2f", engine->renderer.camera.position.x, engine->renderer.camera.position.y, engine->renderer.camera.position.z);
        
        int total_faces = 0;
        int mis_count = 0;
        const ECS* ecs = &engine->ecs;
        const System* render_system = &ecs->systems[engine->render_system_id];
        for (int si = 0; si < render_system->num_archetypes; ++si)
        {
            const ArchetypeID archetype_id = render_system->archetype_ids[si];
            Archetype* archetype = &ecs->archetypes[archetype_id];

            MeshInstance* mis = Archetype_get_component_list(archetype, COMPONENT_MeshInstance);

            for (int i = 0; i < archetype->entity_count; ++i)
            {
                MeshInstance* mi = &mis[i];

                const Scene* scene = &engine->scene;
                const MeshBase* mb = &scene->mesh_bases.bases[mi->mb_id];
                total_faces += mb->num_faces;
                ++mis_count;
            }
        }

        snprintf(vertices_str, sizeof(vertices_str), "VERTICES: %d", total_faces * 3);
    }
}

void engine_destroy(Engine* engine)
{
    // TODO: this stuff.
    ui_destroy(&engine->ui);
    window_destroy(&engine->window);
}

void engine_handle_input(Engine* engine, float dt)
{
    Camera* camera = &engine->renderer.camera;

    if (!engine->handle_input)
    {
        return;
    }

    // TODO: Could move this to an input handler or something. 
    //       Not sure if necessary.
    POINT mouse_position;
    GetCursorPos(&mouse_position);

    int rel_x = engine->window.mouse_dx;
    int rel_y = engine->window.mouse_dy;

    // Reset the delta mouse position as we're now responding to it.
    engine->window.mouse_dx = 0;
    engine->window.mouse_dy = 0;

    // Calculate the camera movement in radians.
    float sens = 0.05f; // TODO: Set this somewhere. 
                        //       Maybe this is where the input manager comes in. If this is to be an actual game engine, user needs control.
    float dy = radians(rel_x * sens);
    float dp = radians(rel_y * sens);

    camera->yaw -= dy;

    // Reset yaw after full 360 spin, accounting for amount past.
    const float TWO_PI = PI * 2;
    if (camera->yaw < -TWO_PI)
    {
        camera->yaw += TWO_PI;
    }
    else if (camera->yaw > TWO_PI)
    {
        camera->yaw -= TWO_PI;
    }

    camera->pitch -= dp;

    // Clamp pitch so we don't look past 90 degrees behind the camera.
    float maxPitch = PI_2 - 0.001f;
    camera->pitch = min(max(camera->pitch, -maxPitch), maxPitch);

    float cosPitch = cosf(camera->pitch);

    // Calculate the camera's direction.
    camera->direction.x = sinf(camera->yaw) * cosPitch;
    camera->direction.y = sinf(camera->pitch);
    camera->direction.z = cosf(camera->yaw) * cosPitch;
    normalise(&camera->direction);

    // TODO: How do I make the engine actually m/s?

    // Direct position changes must be multipled by dt.
    const float SPEED = 10.f;
    float meters_per_second = SPEED * dt;

    // Process keyboard input.
    BYTE keys[256];
    if (!GetKeyboardState(keys))
    {
        // TODO: Handle error?
        log_error("Failed to get keyboard state.");
        return;
    }

    const int KeyDown = 0x80;

    if (keys['W'] & KeyDown)
    {
        v3_add_eq_v3(&camera->position, v3_mul_f(camera->direction, meters_per_second));
    }
    if (keys['S'] & KeyDown)
    {
        v3_sub_eq_v3(&camera->position, v3_mul_f(camera->direction, meters_per_second));
    }
    if (keys['A'] & KeyDown)
    {
        V3 up = { 0, 1, 0 };
        V3 right = normalised(cross(camera->direction, up));

        v3_sub_eq_v3(&camera->position, v3_mul_f(right, meters_per_second));
    }
    if (keys['D'] & KeyDown)
    {
        V3 up = { 0, 1, 0 };
        V3 right = normalised(cross(camera->direction, up));

        v3_add_eq_v3(&camera->position, v3_mul_f(right, meters_per_second));
    }
    if (keys[VK_LSHIFT] & KeyDown)
    {
        camera->position.y -= meters_per_second;
    }
    if (keys[VK_SPACE] & KeyDown)
    {
        camera->position.y += meters_per_second;
    }
}

// Window events
static void engine_on_resize(void* ctx)
{
    Engine* engine = (Engine*)ctx;

    Status status = renderer_resize(&engine->renderer, 
        (int)(engine->window.width / engine->upscaling_factor), 
        (int)(engine->window.height / engine->upscaling_factor));

    // TODO: Feels wrong setting the window bitmap dimensions here instead
    //       of in the window. But need the upscaling information. Maybe
    //       the window could handle the upscaling factor instead?

    engine->window.bitmap.bmiHeader.biWidth = engine->renderer.target.canvas.width;
    engine->window.bitmap.bmiHeader.biHeight = -engine->renderer.target.canvas.height;

    if (STATUS_OK != status)
    {
        // TODO: What to do here?
    }
}

static void engine_process_keyup(void* ctx, WPARAM wParam) 
{
    Engine* engine = (Engine*)ctx;

    // Handle any engine specific keybinds, if
    // not engine specific, pass to user callback.
    switch (wParam)
    {
    case VK_TAB:
    {
        ShowCursor(engine->handle_input);
        engine->handle_input = !engine->handle_input;

        if (engine->handle_input)
        {
            RECT rect = { 0 };
            GetClientRect(engine->window.hwnd, &rect);

            POINT center = { 0 };
            center.x = (rect.left + rect.right) / 2;
            center.y = (rect.top + rect.bottom) / 2;

            // Reset the cursor to the center of the screen.
            ClientToScreen(engine->window.hwnd, &center);
            SetCursorPos(center.x, center.y);

            // Reset the delta mouse movement.
            engine->window.mouse_dx = 0;
            engine->window.mouse_dy = 0;

            // Restrict the cursor to the center of the screen.
            RECT cursor_area = { 0 };
            cursor_area.left = center.x;
            cursor_area.top = center.y;
            cursor_area.right = center.x;
            cursor_area.bottom = center.y;

            ClipCursor(&cursor_area);
        }
        else
        {
            // Release the cursor so it can move freely..
            ClipCursor(NULL);
        }

        break;
    }
    case VK_ESCAPE:
    {
        engine->running = 0;
        PostQuitMessage(0);

        break;
    }
    default:
    {
        // Fire the user defined callback.
        // TODO: Could be nice to not use WPARAM for this.
        engine_on_keyup(engine, wParam);
    }
    }
}

static void engine_process_lmbdown(void* ctx)
{
    Engine* engine = (Engine*)ctx;

    engine_on_lmbdown(ctx);
}