#include "text.h"

#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

Text text_create(char* str, int x, int y, int colour, int scale)
{
	// TODO: text_init instead and modify pointer?
	Text text = {
		.text = str,
		.x = x,
		.y = y,
		.colour = colour,
		.scale = scale
	};

	return text;
}

void text_draw(Canvas* canvas, Text* text, Font* font, float upscaling_factor)
{
	// TODO: All ui should be scaled the same so we can at least keep 
	// consistency with the ui.
	// TODO: I want the text to stay as similar as possible when upscaling.
	
	int scale = (int)fmaxf(text->scale / upscaling_factor, 1.f);

	int scaled_x = (int)(text->x / upscaling_factor);
	int scaled_y = (int)(text->y / upscaling_factor);

	int start_pixel = scaled_y * canvas->width + scaled_x;

	// For scaling the font, we draw the pixel as a rect with
	// width and height = scale.
    int* canvas_data = canvas->pixels.data;

	int count = 0;
	int row = 0;
	for (int ci = 0; ci < strlen(text->text); ++ci)
	{
		char c = text->text[ci];
		// For a space, simply move to the next character.
		if (c == ' ')
		{
			count++;
			continue;
		}

		// For a newline, move to the start of the next row.
		if (c == '\n')
		{
			count = 0;
			row++;
			continue;
		}

		// Calculate the start pixel of the char we're drawing, 
		// taking into account for the number of characters drawn and 
		// the space between them.
		int row_offset = (row * (font->char_height + 1) * scale) * canvas->width;
		int target_row = start_pixel + (font->char_width + 1) * scale * count + row_offset;

		int src_offset = Font_get_char_offset(font, c);

		if (src_offset == -1)
		{
			++count;
			continue;
		}
			
		// For each pixel in the char bitmap, if it has colour, scale and 
		// draw it.
		for (int y = 0; y < font->char_height; ++y)
		{
			int target_pixel = target_row;

			for (int x = 0; x < font->char_width; ++x)
			{
				int i = target_pixel;

				// Only draw the filled in pixels from the source.
				if (font->atlas.pixels.data[src_offset++])
				{
					for (int y2 = 0; y2 < scale; ++y2)
					{
						for (int x2 = 0; x2 < scale; ++x2)
						{
							canvas_data[i + x2] = text->colour;
						}

						i += canvas->width;
					}
				}

				// Move the target pixel along by the pixel width.
				target_pixel += scale;
			}

			// Increment the target row.
			target_row += canvas->width * scale;
		}
		count++;
	}
}