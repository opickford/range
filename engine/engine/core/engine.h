#ifndef ENGINE_H
#define ENGINE_H

#include "scene.h"

#include "window.h"
#include "resources.h"
#include "components.h"

#include "ui/ui.h"

// TODO: These need to be refactored.
#include "renderer/renderer.h"
#include "renderer/render.h"

#include "common/status.h"

#include <cecs/ecs.h>

typedef enum input_mode
{
    INPUT_MODE_UI,
    INPUT_MODE_NOCLIP,
    INPUT_MODE_GAME,
    INPUT_MODE_INVALID
} input_mode_t;

typedef struct
{
    cecs_t* ecs;

    // TODO: Defining these in the engine struct feels awful, should be a static?
    cecs_view_id_t render_view_id;
    cecs_view_id_t lighting_view_id;
    cecs_view_id_t physics_view_id;

    cecs_view_id_t moving_collider_view_id;
    cecs_view_id_t static_collider_view_id;

	// engine_t components.
	window_t window;
	ui_t ui;
	renderer_t renderer;
	resources_t resources; // Works fine for now, potentially something to refactor.
    physics_t physics;

	// Scene data. - I don't think the engine needs to manage multiple.
    // TODO: Only manage one scene.
	scene_t scene;

    // TODO: Stuff like this is private?

	// engine_t settings
	int running;

	// TODO: Move these somewhere?
    input_mode_t input_mode;
	float upscaling_factor;

	// TODO: Allow the user to set callbacks just like the window class.

} engine_t;

// Main API
status_t engine_init(engine_t* engine, int window_width, int window_height);

void engine_run(engine_t* engine);

void engine_destroy(engine_t* engine);

// Public engine events that the game should define.
void engine_on_init(engine_t* engine);

void engine_on_update(engine_t* engine, float dt);

void engine_on_keyup(engine_t* engine, WPARAM wParam);

void engine_on_lmbdown(engine_t* engine);

// Internal functions


// TODO: Some sort of input handler? Fine here for now.
void engine_handle_input(engine_t* engine, float dt);

// Private window events.
static void engine_on_resize(void* ctx);

static void engine_process_keyup(void* ctx, WPARAM wParam);
static void engine_process_lmbdown(void* ctx); // TODO: I don't really like this naming 'lmbdown'

#endif