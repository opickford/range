#ifndef TEXT_H
#define TEXT_H

#include "font.h"
#include "core/canvas.h"

typedef struct
{
	char* text;
	int x;
	int y;
	int colour;
	int scale;

} text_t;

text_t text_create(char* text, int x, int y, int colour, int scale);

void text_draw(canvas_t* canvas, text_t* text, font_t* font, float upscaling_factor);

#endif