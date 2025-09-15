#ifndef ENGINE_H
#define ENGINE_H

#define MAX_SCENES 3

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

typedef struct
{
    ECS ecs;

    // TODO: Defining these in the engine struct feels awful, should be a static?
    SystemID render_system_id;
    SystemID lighting_system_id;

	// Engine components.
	Window window;
	UI ui;
	Renderer renderer;
	Resources resources; // Works fine for now, potentially something to refactor.

	// Scene data. - I don't think the engine needs to manage multiple.
    // TODO: Only manage one scene.
	Scene scene;

    // TODO: Stuff like this is private?

	// Engine settings
	int running;

	// TODO: Move these somewhere?
	int handle_input;
	float upscaling_factor;

	// TODO: Allow the user to set callbacks just like the window class.

} Engine;

// Main API
Status engine_init(Engine* engine, int window_width, int window_height);

void engine_run(Engine* engine);

void engine_destroy(Engine* engine);

// Public engine events that the game should define.
void engine_on_init(Engine* engine);

void engine_on_update(Engine* engine, float dt);

void engine_on_keyup(Engine* engine, WPARAM wParam);

// Internal functions


// TODO: Some sort of input handler? Fine here for now.
void engine_handle_input(Engine* engine, float dt);

// Private window events.
static void engine_on_resize(void* ctx);

static void engine_process_keyup(void* ctx, WPARAM wParam);

#endif