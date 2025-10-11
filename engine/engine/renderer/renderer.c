#include "renderer.h"

status_t renderer_init(renderer_t* renderer, int width, int height)
{
	// Initialise the render target.
	status_t status = render_target_init(&renderer->target, width, height);
	if (STATUS_OK != status)
	{
		return status;
	}

	// Initialise the render settings.
	memset(&renderer->settings, 0, sizeof(render_settings_t));

	renderer->settings.fov = 90.f;
	renderer->settings.near_plane = 1.f;
	renderer->settings.far_plane = 100.f;

	update_projection_m4(&renderer->settings, width / (float)height);

	// Create a camera.
	memset(&renderer->camera, 0, sizeof(camera_t));
	renderer->camera.direction.x = 0;
	renderer->camera.direction.y = 0;
	renderer->camera.direction.z = -1.f;

	renderer->camera.position.x = 0;
	renderer->camera.position.y = 0;
	renderer->camera.position.z = 0;

	renderer->camera.yaw = PI; // TODO: Gotta be able to set one or the other or have them update together...

	// Create the view frustum.
	view_frustum_init(&renderer->settings.view_frustum, renderer->settings.near_plane, renderer->settings.far_plane, renderer->settings.fov,
		renderer->target.canvas.width / (float)(renderer->target.canvas.height));

	return STATUS_OK;
}

status_t renderer_rev3_size(renderer_t* renderer, int width, int height)
{
	// Resize the render target.
	status_t status = render_target_rev3_size(&renderer->target, width, height);
	if (STATUS_OK != status)
	{
		return status;
	}

	// Update the projection matrix.
	update_projection_m4(&renderer->settings, width / (float)height);

	// Recreate the view frustum.
	view_frustum_init(&renderer->settings.view_frustum, renderer->settings.near_plane, renderer->settings.far_plane, renderer->settings.fov,
		renderer->target.canvas.width / (float)(renderer->target.canvas.height));


	return STATUS_OK;
}
