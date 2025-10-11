#ifndef UI_H
#define UI_H


#include "font.h"
#include "text.h"

#include "core/canvas.h"

#include "common/status.h"

// TODO: Why have i done this.
#define MAX_TEXT 100

typedef struct
{
	// Global ui_t
	canvas_t* canvas;
	font_t font;

	// Specific widgets
	int text_count;
	text_t text[MAX_TEXT];

} ui_t;

status_t ui_init(ui_t* ui, canvas_t* canvas);

void ui_draw(ui_t* ui, float upscaling_factor);

void ui_destroy(ui_t* ui);

// Functions for adding ui_t widgets?

#endif