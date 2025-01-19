#include "camera.h"

#include "maths/matrix4.h"
#include "maths/vector3.h"

void calculate_view_matrix(const Camera* camera, M4 out)
{
	// The view matrix transforms world space positions to view space. 
	// Rather than rotating and translating the camera around the 
	// world, the view matrix transforms vectors around the camera.
	// Therefore, we create a look at matrix that takes the inverse
	// of the camera's position and direction.
	look_at(v3_mul_f(camera->position, -1.f), v3_mul_f(camera->direction, -1.f), out);
}
