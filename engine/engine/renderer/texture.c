#include "texture.h"

#include "common/colour.h"
#include "utils/logger.h"

#include <Windows.h>

Status texture_load_from_bmp(Texture* texture, const char* file)
{
    // TODO: Do we need to use Windows.h here? If we're only loading
    //       bitmaps, the image loading code could be quite simple.

    // Initialise the texture.
    memset(texture, 0, sizeof(Texture));

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

    texture->width = bitmap.bmWidth;
    texture->height = bitmap.bmHeight;

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

    // Read the image as an array of ints
    int* temp_pixels = 0;
    Status status = resize_int_buffer(&temp_pixels, bitmap.bmWidthBytes * bitmap.bmHeight);
    if (STATUS_OK != status || !temp_pixels)
    {
        return status;
    }

    // Get the pixels buffer.
    GetDIBits(mem_hdc, h_bitmap, 0, bitmap.bmHeight, temp_pixels, &bmi, DIB_RGB_COLORS);

    // Cleanup.
    if (!DeleteObject(h_bitmap))
    {
        log_warn("Potential memory leak. Failed to DeleteObject when texture_load_from_bmp.");
    }

    SelectObject(mem_hdc, prev);
    if (!DeleteDC(mem_hdc))
    {
        log_warn("Potential memory leak. Failed to DeleteDC when texture_load_from_bmp.");
    }

    if (!ReleaseDC(NULL, hdc))
    {
        log_warn("Potential memory leak. Failed to ReleaseDC when texture_load_from_bmp.");
    }

    
    //resize_float_buffer(&texture->data, texture->width * texture->height * 3);
    resize_int_buffer(&texture->data, texture->width * texture->height * 3);

    // TODO: Gotta test the texture a few ways, do we want r,g,b uint8? 3 floats? what.
    for (int i = 0; i < texture->width * texture->height; ++i)
    {

        int index = i * 3;
        int colour = temp_pixels[i];


        //int r, g, b;

        //unpack_int_rgb_to_ints(colour, &texture->data[index], &texture->data[index + 1], &texture->data[index + 2]);
        //float r, g, b;

        //unpack_int_rgb_to_floats(colour, &r, &g, &b);
        //printf("%f %f %f\n", r, g, b);

        unpack_int_rgb_to_floats(colour, &texture->data[index], &texture->data[index + 1], &texture->data[index + 2]);
    }

    free(temp_pixels);

    return STATUS_OK;
}

void texture_destroy(Texture* texture)
{
    free(texture->data);
    texture->data = 0;

    free(texture);
    texture = 0; // TODO: Do we want to do this?
}