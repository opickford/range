#ifndef FONT_H
#define FONT_H

#include "common/status.h"

typedef struct
{
	int char_width;
	int char_height;
	int chars_per_row;

	int bitmap_width;

	const char* defined_chars;
	unsigned int* pixels;

} Font;

Status font_init(Font* font);

int font_get_char_index(Font* font, char c);

#endif