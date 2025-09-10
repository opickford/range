#ifndef WINDOW_H
#define WINDOW_H

#include "canvas.h"

#include "common/status.h"

#include <Windows.h>

#define RANGE_WINDOW_CLASS "range_window_class"
#define RANGE_WINDOW_TITLE "range"

// TODO: Is all of this really necessary.

// TODO: Top of file comments.

// Window message handling.
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

int window_process_messages();

typedef struct
{
	Canvas* canvas;
	
	HWND hwnd;
	HDC hdc;
	BITMAPINFO bitmap;

	int width;
	int height;

	// Event callbacks
	void* ctx; // Set to Engine* so we can use it in the callbacks.
	void (*on_resize)(void*);
	void (*on_keyup)(void*, WPARAM);
	
	// Relative mouse movement from raw input.
	int mouse_dx, mouse_dy;
	
} Window;

Status window_init(Window* window, Canvas* canvas, void* ctx, int width, int height);

void window_display(Window* window);

void window_destroy(Window* window);


#endif