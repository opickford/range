#ifndef FONT_H
#define FONT_H

#include "core/canvas.h"
#include "common/status.h"

typedef struct
{
	int char_width;
	int char_height;
	int chars_per_row;

	const char* defined_chars;

    canvas_t atlas;

} font_t;

status_t font_init(font_t* font);

int font_get_char_offset(font_t* font, char c);

#endif