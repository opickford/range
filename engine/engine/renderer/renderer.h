#ifndef RENDERER_H
#define RENDERER_H

#include "render_target.h"
#include "render_settings.h"
#include "frame_data.h"
#include "camera.h"

#include "common/status.h"

typedef struct
{
	render_target_t target;
	render_settings_t settings;
	camera_t camera;


	frame_data_t frame_data;
	
} renderer_t;

status_t renderer_init(renderer_t* renderer, int width, int height);
status_t renderer_rev3_size(renderer_t* renderer, int width, int height);

#endif