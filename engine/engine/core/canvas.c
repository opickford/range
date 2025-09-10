#include "canvas.h"

#include "common/status.h"

#include "utils/logger.h"

#include <Windows.h>

#include <stdlib.h>
#include <string.h>

Status canvas_init(Canvas* canvas, int width, int height)
{
	memset(canvas, 0, sizeof(Canvas));

	canvas->width = width;
	canvas->height = height;

    const size_t length = (size_t)width * height * 4;
    Vector_reserve(canvas->pixels, length);
	
	if (!canvas->pixels.data)
	{
		log_error("Failed to allocate memory for canvas pixels.");
		return STATUS_ALLOC_FAILURE;
	}
    
    // Clear the allocated memory otherwise we might get artificats.
    memset(canvas->pixels.data, 0, length * sizeof(*canvas->pixels.data));

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
    Vector_reserve(canvas->pixels, bitmap.bmWidthBytes * bitmap.bmHeight);
    if (!canvas->pixels.data)
    {
        return STATUS_ALLOC_FAILURE;
    }

    // Get the pixels buffer.
    GetDIBits(mem_hdc, h_bitmap, 0, bitmap.bmHeight, canvas->pixels.data, &bmi, DIB_RGB_COLORS);

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

Status canvas_write_to_bitmap(const Canvas* canvas, const char* file)
{
    // TODO: TEMP: Copied from: https://stackoverflow.com/a/55504419

    // TODO: BMP store data in BGR format.

    static unsigned char infoHeader[] = {
    0,0,0,0, /// header size
    0,0,0,0, /// image width
    0,0,0,0, /// image height - TODO: This should be negative
    0,0, /// number of color planes
    0,0, /// bits per pixel
    0,0,0,0, /// compression
    0,0,0,0, /// image size
    0,0,0,0, /// horizontal resolution
    0,0,0,0, /// vertical resolution
    0,0,0,0, /// colors in color table
    0,0,0,0, /// important color count
    };

    const int width = canvas->width;
    const int height = canvas->height;

    const int bytesPerPixel = 4; /// red, green, blue
    const int fileHeaderSize = 14;
    const int infoHeaderSize = 40;

    infoHeader[0] = (unsigned char)(infoHeaderSize);
    infoHeader[4] = (unsigned char)(width);
    infoHeader[5] = (unsigned char)(width >> 8);
    infoHeader[6] = (unsigned char)(width >> 16);
    infoHeader[7] = (unsigned char)(width >> 24);
    infoHeader[8] = (unsigned char)(height * -1);
    infoHeader[9] = (unsigned char)(height * -1 >> 8);
    infoHeader[10] = (unsigned char)(height * -1 >> 16);
    infoHeader[11] = (unsigned char)(height * -1 >> 24);
    infoHeader[12] = (unsigned char)(1);
    infoHeader[14] = (unsigned char)(bytesPerPixel * 8);

    const int pitch = canvas->width * 4; // 4 bytes per pixel?
    

    unsigned char padding[3] = { 0, 0, 0 };
    int paddingSize = (4 - (pitch) % 4) % 4;

    int fileSize = fileHeaderSize + infoHeaderSize + (pitch + paddingSize) * height;

    static unsigned char fileHeader[] = {
        0,0, /// signature
        0,0,0,0, /// image file size in bytes
        0,0,0,0, /// reserved
        0,0,0,0, /// start of pixel array
    };

    fileHeader[0] = (unsigned char)('B');
    fileHeader[1] = (unsigned char)('M');
    fileHeader[2] = (unsigned char)(fileSize);
    fileHeader[3] = (unsigned char)(fileSize >> 8);
    fileHeader[4] = (unsigned char)(fileSize >> 16);
    fileHeader[5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = (unsigned char)(fileHeaderSize + infoHeaderSize);

    unsigned char* image = (unsigned char*)canvas->pixels.data;

    FILE* imageFile = fopen(file, "wb");

    fwrite(fileHeader, 1, fileHeaderSize, imageFile);
    fwrite(infoHeader, 1, infoHeaderSize, imageFile);

    int i;
    for (i = 0; i < height; i++) {
        fwrite(image + (i * pitch), bytesPerPixel, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
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
    size_t length = width * height;
    Vector_reserve(canvas->pixels, length);

	// Check the allocation worked.
	if (!canvas->pixels.data)
	{
		log_error("Failed to reallocate memory for canvas pixels on resize.");
		return STATUS_ALLOC_FAILURE;
	}

	// Update the canvas.
	canvas->width = width;
	canvas->height = height;

	return STATUS_OK;
}

void canvas_fill(Canvas* canvas, const unsigned int colour)
{
	// TODO: Look for some sort of blit or fill function 
	const int length = canvas->width * canvas->height;	
	unsigned int* ptr = canvas->pixels.data;
	
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
    int* source_data = source->pixels.data;
    int* target_data = target->pixels.data;

    for (int y = 0; y < source->height; ++y)
    {
        for (int x = 0; x < source->height; ++x)
        {
            target_data[(y + y_offset) * target->width + x + x_offset] = source_data[y * source->width + x];
        }
    }
}

void canvas_destroy(Canvas* canvas)
{
    Vector_destroy(canvas->pixels);

	free(canvas);
	canvas = 0; // TODO: Do we want to do this?
}