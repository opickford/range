#include "canvas.h"

#include "common/status.h"

#include "utils/logger.h"
#include "utils/memory_utils.h"

#include <Windows.h>

#include <stdlib.h>
#include <string.h>

Status canvas_init(Canvas* canvas, int width, int height)
{
	memset(canvas, 0, sizeof(Canvas));

	canvas->width = width;
	canvas->height = height;
	canvas->pixels = malloc((size_t)width * height * sizeof(unsigned int) * 4);

	if (!canvas->pixels)
	{
		log_error("Failed to allocate memory for canvas pixels.");
		return STATUS_ALLOC_FAILURE;
	}

	return STATUS_OK;
}

Status canvas_init_from_bitmap(Canvas* canvas, const char* file)
{
    // TODO: Do we need to use Windows.h here? If we're only loading
    //       bitmaps, the image loading code should be quite simple.

    // Initialise the texture.
    memset(canvas, 0, sizeof(Canvas));
        
    // Try load the bitmap.
    HBITMAP h_bitmap = (HBITMAP)LoadImageA(
        NULL,
        file,
        IMAGE_BITMAP,
        0, 0,
        LR_LOADFROMFILE
    );

    if (!h_bitmap)
    {
        log_error("Failed to load texture bitmap.");
        return STATUS_FAILURE;
    }

    // Get bitmap properties.
    BITMAP bitmap = { 0 };
    GetObject(h_bitmap, sizeof(BITMAP), &bitmap);

    canvas->width = bitmap.bmWidth;
    canvas->height = bitmap.bmHeight;

    // Create a compatible DC.
    HDC hdc = GetDC(NULL);
    HDC mem_hdc = CreateCompatibleDC(hdc);

    HGDIOBJ prev = SelectObject(mem_hdc, h_bitmap);

    // Set bitmap info.
    BITMAPINFO bmi;
    memset(&bmi, 0, sizeof(BITMAPINFO));
    bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bmi.bmiHeader.biWidth = bitmap.bmWidth;
    bmi.bmiHeader.biHeight = -bitmap.bmHeight;
    bmi.bmiHeader.biPlanes = 1;
    bmi.bmiHeader.biBitCount = bitmap.bmBitsPixel;
    bmi.bmiHeader.biCompression = BI_RGB; // Uncompressed RGB.

    // Allocate memory for pixels
    Status status = resize_int_buffer(&canvas->pixels, bitmap.bmWidthBytes * bitmap.bmHeight);
    if (STATUS_OK != status)
    {
        return status;
    }

    // Get the pixels buffer.
    GetDIBits(mem_hdc, h_bitmap, 0, bitmap.bmHeight, canvas->pixels, &bmi, DIB_RGB_COLORS);

    // Cleanup.
    if (!DeleteObject(h_bitmap))
    {
        log_warn("Potential memory leak. Failed to DeleteObject when canvas_load_from_bitmap.");
    }

    SelectObject(mem_hdc, prev);
    if (!DeleteDC(mem_hdc))
    {
        log_warn("Potential memory leak. Failed to DeleteDC when canvas_load_from_bitmap.");
    }

    if (!ReleaseDC(NULL, hdc))
    {
        log_warn("Potential memory leak. Failed to ReleaseDC when canvas_load_from_bitmap.");
    }

    return STATUS_OK;
}

Status canvas_resize(Canvas* canvas, int width, int height)
{
	// Check the size has changed.
	if (canvas->width == width && canvas->height == height)
	{
		return STATUS_OK;
	}
	
	// Allocate memory for the new array.
	// TODO: Use my memory allocating helpers for this.
	unsigned int* new_pixels = realloc(canvas->pixels, (size_t)width * height * sizeof(unsigned int));

	// Check the allocation worked.
	if (!new_pixels)
	{
		log_error("Failed to reallocate memory for canvas pixels on resize.");
		return STATUS_ALLOC_FAILURE;
	}

	// Update the canvas.
	canvas->pixels = new_pixels;
	canvas->width = width;
	canvas->height = height;

	return STATUS_OK;
}

void canvas_fill(Canvas* canvas, const unsigned int colour)
{
	// TODO: Look for some sort of blit or fill function 
	const int length = canvas->width * canvas->height;	
	unsigned int* ptr = canvas->pixels;
	
	unsigned int i = length;
	
	while (i)
	{
		*ptr = colour;
		--i;
		++ptr;
	}
}

void canvas_draw(const Canvas* source, Canvas* target, int x_offset, int y_offset)
{
    for (int y = 0; y < source->height; ++y)
    {
        for (int x = 0; x < source->height; ++x)
        {
            target->pixels[(y + y_offset) * target->width + x + x_offset] = source->pixels[y * source->width + x];
        }
    }
}

void canvas_destroy(Canvas* canvas)
{
	free(canvas->pixels);
	canvas->pixels = 0;

	free(canvas);
	canvas = 0; // TODO: Do we want to do this?
}