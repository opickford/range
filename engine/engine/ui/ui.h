#ifndef UI_H
#define UI_H

#include "canvas.h"

#include "font.h"
#include "text.h"

#include "common/status.h"

#define MAX_TEXT 10

typedef struct
{
	// Global UI
	Canvas* canvas;
	Font font;

	// Specific widgets
	int text_count;
	Text text[MAX_TEXT];

} UI;

Status ui_init(UI* ui, Canvas* canvas);

void ui_draw(UI* ui, float upscaling_factor);

void ui_destroy(UI* ui);

// Functions for adding UI widgets?

#endif