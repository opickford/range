#ifndef RENDER_TARGET_H
#define RENDER_TARGET_H

#include "core/canvas.h"

#include "utils/logger.h"

#include <stdlib.h>
#include <string.h>

// TODO: Consider global rendertarget. We're only ever going to have 1.

// Contains the different buffers needed for rendering.
typedef struct
{
	Canvas canvas;
	float* depth_buffer;

} RenderTarget;

// TODO: .c file?
inline Status render_target_init(RenderTarget* rt, const int width, const int height)
{
    memset(rt, 0, sizeof(RenderTarget));

    Status status = Canvas_init(&rt->canvas, width, height);
    if (STATUS_OK != status)
    {
        return status;
    }

    // TODO: Refactor to be a DepthBuffer.
    rt->depth_buffer = malloc((size_t)width * height * sizeof(float));

    if (!rt->depth_buffer)
    {
        log_error("Failed to allocate memory for depth buffer.");
        return STATUS_ALLOC_FAILURE;
    }

    return STATUS_OK;
}

inline Status render_target_resize(RenderTarget* rt, int width, int height)
{
    Status status = canvas_resize(&rt->canvas, width, height);

    if (STATUS_OK != status)
    {
        return status;
    }

    // Allocate memory for the new array.
    // TODO: Use my memory allocating helpers for this.
    float* new_db = realloc(rt->depth_buffer, (size_t)width * height * sizeof(float));

    // Check the allocation worked.
    if (!new_db)
    {
        log_error("Failed to reallocate memory for rt depth buffer on resize.");
        return STATUS_ALLOC_FAILURE;
    }

    // Update the depth buffer.
    rt->depth_buffer = new_db;

    return STATUS_OK;
}

inline void render_target_clear(RenderTarget* rt, uint32_t bg_colour)
{
    const int length = rt->canvas.width * rt->canvas.height;

    uint32_t* canvas_ptr = rt->canvas.pixels;
    float* depth_buffer_ptr = rt->depth_buffer;

    // Splitting these loops is much faster, guessing because of cache.
    // After extensive testing of how to fill these, can't figure out a better
    // way than this.
    for (int i = 0; i < length; ++i)
    {
        canvas_ptr[i] = bg_colour;
    }

    for (int i = 0; i < length; ++i)
    {
        depth_buffer_ptr[i] = 1.f; // Projection should map z to a max of 1.
    }
}


inline void render_target_destroy(RenderTarget* rt)
{
    canvas_destroy(&rt->canvas);

    free(rt->depth_buffer);
    rt->depth_buffer = 0;

    // TODO: Do i need to do rt = 0; here?? Not sure.
}

#endif