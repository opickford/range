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

    Canvas atlas;

} Font;

Status Font_init(Font* font);

int Font_get_char_offset(Font* font, char c);

#endif