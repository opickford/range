#ifndef RENDER_SETTINGS_H
#define RENDER_SETTINGS_H

#include "frustum_culling.h"

#include "core/canvas.h"

#include "maths/matrix4.h"
#include "maths/utils.h"

// TODO: render_settings_t is a bit off. We also want to store the view frustum
//		 as this won't change unless fov changes. This stuff is slightly 
//		 different to the fov/near/farplane.
//			
typedef struct
{
	float fov;
	float near_plane; 

	// TODO: ^^^ Will have to experiment. We want this as large as possible without clipping too close. Apparently
	//	   Large outdoor scene: nearplane 0.1-1, far 500-1000/
	//	   Small indoor scene: nearplane 0.1-0.01, 50-100. 
	//float nearPlane = near_plane_dist * -1.f;

	float far_plane;

	// TODO: Should these go to the renderer_t?
	m4_t projection_matrix;
	view_frustum_t view_frustum; // TODO: Definitely should go in the renderer.

} render_settings_t;

inline void update_projection_m4(render_settings_t* rs, float aspect_ratio)
{
	// TODO: Get rid of this helper function.
	m4_projection(rs->fov, aspect_ratio, rs->near_plane, rs->far_plane, rs->projection_matrix);
}

#endif