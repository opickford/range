#include "lights.h"

#include "strides.h"

#include "renderer/render_buffers.h"

#include "utils/memory_utils.h"
#include "utils/logger.h"

#include <string.h>
#include <stdlib.h>

Status PointLight_init(PointLight* light)
{
    memset(light, 0, sizeof(PointLight));
    light->strength = 1.f;

    return STATUS_OK;
}

void PointLight_destroy(PointLight* light)
{
}


Status ShadowCastingPointLight_init(ShadowCastingPointLight* light)
{
    memset(light, 0, sizeof(ShadowCastingPointLight));

    // TODO: Create shadow maps.
    //		 All temporary for now.
    // TODO: No idea what size would be best.
    const int RES = 256;

    // TODO: Need cubemap for point light.
    // TODO: Potentially a depth buffer the same aspect ratio as the window could
    //		 give better results. Not sure just a thought.

    // Create a new depth buffer for the new light.
    // TODO: How do we create the 6 maps?
    DepthBuffer* temp = realloc(light->depth_maps, (size_t)light->count * sizeof(DepthBuffer));
    if (!temp)
    {
        log_error("Failed to realloc for point_lights->depth_maps.");
        // TODO: Return status.
        return;
    }
    point_lights->depth_maps = temp;

    depth_buffer_init(&point_lights->depth_maps[point_lights->count - 1], RES, RES);
    depth_buffer_fill(&point_lights->depth_maps[point_lights->count - 1], 1.f);
}


