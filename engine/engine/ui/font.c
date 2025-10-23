#include "font.h"

#include "utils/logger.h"

#include <Windows.h>

static int get_char_index(const font_t* font, char c)
{
    // TODO: Switch to a map of precalculated offsets possibly?
    int defined = 0;
    int char_index;

    for (char_index = 0; char_index < strlen(font->defined_chars); ++char_index)
    {
        if (font->defined_chars[char_index] == c)
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

    return char_index;
}

static int get_initial_atlas_char_offset(const font_t* font, 
    const int initial_atlas_width, char c)
{
    // Returns the offset to the char in the given font atlas, therefore,
    // requires the initial atlas width.

    int char_index = get_char_index(font, c);
    if (-1 == char_index) return -1;

    // Calculate the position of the char on the bitmap.
    int cy = (char_index / font->chars_per_row);
    int cx = char_index - font->chars_per_row * cy;

    // Add one to the charHeight as there is 1 pixel between rows.
    int rowOffset = cy * (font->char_height + 1) * initial_atlas_width;

    // Add one to the charWidth as there is 1 pixel between characters.
    int colOffset = cx * (font->char_width + 1);

    return rowOffset + colOffset;
}

status_t font_init(font_t* font)
{
    // TODO: There are strict rules about the input axis, should write these
    //       out properly!!! pixel gap between rows/cols etc
    
    // Initialise the font struct.
	memset(font, 0, sizeof(font_t));

	// Initialise the character data.
    // TODO: Allow this to be set?
	font->defined_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()[]{}<>/*:#%!?.,'\"@&$";
	font->chars_per_row = 13;
	font->char_width = 5;
	font->char_height = 9;

    // Read the given font atlas.
    canvas_t src;
    // TODO: Stop hardcoding this path.
    status_t status = canvas_init_from_bitmap(&src, "C:/Users/olive/source/repos/csrge/res/fonts/minogram_6x10_font.bmp");
    if (STATUS_OK != status)
    {
        log_error("Failed to load font atlas bitmap, status: %s", status_to_str(status));
        return status;
    }

    // The atlas will be transformed to a flat array with no padding.
    status = canvas_init(&font->atlas, 
        (size_t)(font->char_width * font->char_height) * strlen(font->defined_chars),
        1);

    // TODO: Should this sort of boilerplate error handling code simply be assertions?
    //       Or potentially some define like RETURN_IF_FAILED(msg, status, ...)?
    if (STATUS_OK != status)
    {
        log_error("Failed to init font atlas canvas, status: %s", status_to_str(status));
        return status;
    }

    // Rotate the atlas so that a whole char is continuous in memory, also
    // removes any unnecessary rows/cols between src atlas chars.
    int out = 0;
    for (int i = 0; i < strlen(font->defined_chars); ++i)
    {
        int src_row = get_initial_atlas_char_offset(font, src.width,
            font->defined_chars[i]);

        for (int j = 0; j < font->char_height; ++j)
        {
            for (int k = 0; k < font->char_width; ++k)
            {
                font->atlas.pixels[out++] = src.pixels[src_row + k];
            }

            // Increment row
            src_row += src.width;
        }
    }

	return STATUS_OK;
}

int font_get_char_offset(font_t* font, char c)
{
    // TODO: Switch to a map of precalculated offsets possibly?
    int char_index = get_char_index(font, c);
    if (-1 == char_index) return -1;
    
    // Chars are stored simply in a flat array.
    return font->char_width * font->char_height * char_index;
}
