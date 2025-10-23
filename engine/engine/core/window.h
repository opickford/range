#ifndef WINDOW_H
#define WINDOW_H

#include "canvas.h"

#include "common/status.h"

#include <Windows.h>

// TODO: Do we like the RANGE_ prefix? should we use something else?
// TODO: RANGE_ is an existing game engine name!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// TODO: I like the name SCOPE for the game engine (scope creep.)

#define RANGE_WINDOW_CLASS "range_window_class"
#define RANGE_WINDOW_TITLE "range"

// TODO: Should this be abstracted to some input file/functionality?
#define RANGE_WND_KEYDOWN 0x80
#define RANGE_KEYDOWN(k) ((k) & RANGE_WND_KEYDOWN)


// TODO: Is all of this really necessary.

// TODO: Top of file comments.

// window_t message handling.
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

int window_process_messages();

typedef struct
{
	canvas_t* canvas;
	
	HWND hwnd;
	HDC hdc;
	BITMAPINFO bitmap;

	int width;
	int height;

	// Event callbacks
	void* ctx; // Set to engine_t* so we can use it in the callbacks.
	void (*on_resize)(void*);
	void (*on_keyup)(void*, WPARAM);
	void (*on_lmbdown)(void*);
	
	// Relative mouse movement from raw input.
	int mouse_dx, mouse_dy;

    // Keyboard state
    BYTE keys[256];
	
} window_t;

status_t window_init(window_t* window, canvas_t* canvas, void* ctx, int width, int height);

void window_display(window_t* window);

void window_destroy(window_t* window);


#endif