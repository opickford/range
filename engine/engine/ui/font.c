#include "font.h"

#include "utils/logger.h"

#include <Windows.h>

Status font_init(Font* font)
{
	// Initialise the font struct.
	memset(font, 0, sizeof(Font));

	// Initialise the character data.
	font->defined_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()[]{}<>/*:#%!?.,'\"@&$";
	font->chars_per_row = 13;
	font->char_width = 5;
	font->char_height = 9;

    // TODO: Stop hardcoding this.
    Status status = Canvas_init_from_bitmap(&font->atlas, "C:/Users/olive/source/repos/range/res/fonts/minogram_6x10_font.bmp");
    if (STATUS_OK != status)
    {
        log_error("Failed to load font bitmap, status: %s", Status_to_str(status));
        return status;
    }

	return STATUS_OK;
}

int font_get_char_index(Font* font, char c)
{
	// TODO: Switch to a map of precalculated offsets possibly?
	int defined = 0;
	int charIndex;

	for (charIndex = 0; charIndex < strlen(font->defined_chars); ++charIndex)
	{
		if (font->defined_chars[charIndex] == c)
		{
			defined = 1;
			break;
		}
	}

	if (!defined)
	{
		log_error("Char not defined: %c", c);
		return -1;
	}

	// Calculate the position of the char on the bitmap.
	int cy = (charIndex / font->chars_per_row);
	int cx = charIndex - font->chars_per_row * cy;

	// Add one to the charHeight as there is 1 pixel between rows.
	int rowOffset = cy * (font->char_height + 1) * font->atlas.width;

	// Add one to the charWidth as there is 1 pixel between characters.
	int colOffset = cx * (font->char_width + 1);

	return rowOffset + colOffset;
}
