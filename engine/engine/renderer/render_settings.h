#ifndef RENDER_SETTINGS_H
#define RENDER_SETTINGS_H

#include "canvas.h"
#include "frustum_culling.h"

#include "maths/matrix4.h"
#include "maths/utils.h"

// TODO: RenderSettings is a bit off. We also want to store the view frustum
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

	// TODO: Should these go to the Renderer?
	M4 projection_matrix;
	ViewFrustum view_frustum; // TODO: Definitely should go in the renderer.

} RenderSettings;

inline void update_projection_m4(RenderSettings* rs, float aspect_ratio)
{
	// TODO: Get rid of this helper function.
	m4_projection(rs->fov, aspect_ratio, rs->near_plane, rs->far_plane, rs->projection_matrix);
}

#endif