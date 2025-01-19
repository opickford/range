#include "matrix4.h"

#include "utils.h"

void m4_mul_m4(const M4 m0, const M4 m1, M4 out)
{
	// Post-multiplication, this means the combined out 
	// matrix will first apply m1, then m0.
	
	// For example, if m0 was a scale matrix, and m1 a 
	// translation matrix. m4_mul_m4(translation, scale, out)
	// would result in a combined transform matrix that
	// scales and then translates the point.

	// We use post-multiplication because we are using 
	// column major vectors/matrices and this is the 
	// industry standard as far as I know. Also, using
	// post-multiplication here aligns with the M * v
	// that we use for m4_mul_v4.

	for (int c = 0; c < 4; ++c)
	{
		for (int r = 0; r < 4; ++r)
		{
			int rowOffset = r * 4;

			out[rowOffset + c] = m1[rowOffset] * m0[c] +
								 m1[rowOffset + 1] * m0[c + 4] +
								 m1[rowOffset + 2] * m0[c + 8] +
								 m1[rowOffset + 3] * m0[c + 12];
		}
	}
}

// TODO: Rename out? Any performance critical code 
//		 can't be returning, takes too long.
void m4_mul_v4(const M4 m, V4 v, V4* out)
{
	out->x = m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12] * v.w;
	out->y = m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13] * v.w;
	out->z = m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * v.w;
	out->w = m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15] * v.w;
}

void m4_identity(M4 out)
{
	out[0] = 1;
	out[1] = 0;
	out[2] = 0;
	out[3] = 0;
	out[4] = 0;
	out[5] = 1;
	out[6] = 0;
	out[7] = 0;
	out[8] = 0;
	out[9] = 0;
	out[10] = 1;
	out[11] = 0;
	out[12] = 0;
	out[13] = 0;
	out[14] = 0;
	out[15] = 1;
}

// TODO: Order of in/out?
void m4_translation(V3 position, M4 out)
{
	m4_identity(out);
	out[12] = position.x;
	out[13] = position.y;
	out[14] = position.z;
}

void m4_rotation(const float pitch, const float yaw, const float roll, M4 out)
{
	// TODO: Look into quarternions. https://en.wikipedia.org/wiki/Quaternion

	// In a right handed coordinate system, the rotations are counter clockwise 
	// around the axis when looking at the axis. So for example, a positive yaw
	// would make whatever we rotate turn left.
	const float sinPitch = sinf(pitch);
	const float sinYaw = sinf(yaw);
	const float sinRoll = sinf(roll);

	const float cosPitch = cosf(pitch);
	const float cosYaw = cosf(yaw);
	const float cosRoll = cosf(roll);

	// TODO: I think labelling these as rotations around axis would make more sense tbf.

	// Rotation around the x axis.
	M4 pitch_rot;
	m4_identity(pitch_rot);

	pitch_rot[5] = cosPitch;
	pitch_rot[6] = sinPitch;
	pitch_rot[9] = -sinPitch;
	pitch_rot[10] = cosPitch;

	// Rotation around the y axis.
	M4 yaw_rot;
	m4_identity(yaw_rot);
	
	yaw_rot[0] = cosYaw;
	yaw_rot[2] = -sinYaw;
	yaw_rot[8] = sinYaw;
	yaw_rot[10] = cosYaw;

	// Rotation around the z axis.
	M4 roll_rot;
	m4_identity(roll_rot);

	roll_rot[0] = cosRoll;
	roll_rot[1] = sinRoll;
	roll_rot[4] = -sinRoll;
	roll_rot[5] = cosRoll;

	// TODO: The gimbal lock is actually quite noticeable. Maybe setting direction would be better?

	// Combine the rotations.
	m4_mul_m4(yaw_rot, pitch_rot, out); 
	// Seems to be the correct order, however, I believe when we come to use
	// quarternions, we won't be using pitch/yaw anymore. Which will be nice.

	// TODO: We want to obviously apply roll, however, we will get gimbal lock.
	//		 Switch to quaternions for this I believe. Also, gotta understand
	//		 order of matrix multplications a bit better.
}

void look_at(V3 position, V3 direction, M4 out)
{
	// TODO: This isn't quite right, most look_at matrices take to and from positions,
	//		 then the camera would get set to the from position. 
	//		 This is slightly misleading because we're really just making a rotation/translation
	//		 matrix from a position and direction.


	V3 world_up = { 0, 1, 0 };

	V3 x_axis;
	if (fabsf(direction.y) == fabsf(world_up.y))
	{
		// If direction.y == -1, the cross product will return 0,0,0. 
		// So hardcode the x axis to the world right?
		x_axis = (V3){ 1,0,0 };
	}
	else
	{
		// Calculate the other axis by using the z axis as the direction.
		x_axis = normalised(cross(world_up, direction));
	}

	V3 y_axis = cross(direction, x_axis);

	// Set the out matrix to the combined translation and rotation matrix.
	m4_identity(out);

	out[0] = x_axis.x; out[1] = y_axis.x; out[2] = direction.x;
	out[4] = x_axis.y; out[5] = y_axis.y; out[6] = direction.y;
	out[8] = x_axis.z; out[9] = y_axis.z; out[10] = direction.z;

	out[12] = dot(x_axis, position);
	out[13] = dot(y_axis, position);
	out[14] = dot(direction, position);
}

void m4_model_matrix(V3 position, V3 eulers, V3 scale, M4 out)
{
	// TODO: Eventually gonna switch to quarternions for rotations.

	// TODO: Look into avoiding the matrix multiplications? Can I set translation without multiplying?
	//		 Pretty sure I can.
	
	M4 translation_m4;
	m4_translation(position, translation_m4);

	M4 rotation_m4;
	m4_rotation(eulers.x, eulers.y, eulers.z, rotation_m4);
	
	M4 scale_m4;
	m4_identity(scale_m4);
	scale_m4[0] = scale.x;
	scale_m4[5] = scale.y;
	scale_m4[10] = scale.z;
	scale_m4[15] = 1;

	// We have to define an output matrix each time.
	// Although in my opinion this is fine it makes it more clear.
	M4 translation_rotation_m4;

	// We use post-matrix multiplication, so here, we end up
	// scaling, then rotating, then translating.
	m4_mul_m4(translation_m4, rotation_m4, translation_rotation_m4);
	m4_mul_m4(translation_rotation_m4, scale_m4, out);
}

void m4_normal_matrix(V3 eulers, V3 scale, M4 out)
{
	// Create a normal matrix from the given eulers and scale.
	// Essentially no translation, keep the rotation, and inverse scale.
	M4 rotation_m4;
	m4_rotation(eulers.x, eulers.y, eulers.z, rotation_m4);

	M4 scale_m4;
	m4_identity(scale_m4);
	scale_m4[0] = 1.f / scale.x;
	scale_m4[5] = 1.f / scale.y;
	scale_m4[10] = 1.f / scale.z;

	m4_mul_m4(rotation_m4, scale_m4, out);
}


void m4_transposed(const M4 in, M4 out)
{
	// Flip the matrix along the diagonal. Essentially column major to row major.
	out[0] = in[0];
	out[1] = in[4];
	out[2] = in[8];
	out[3] = in[12];
	out[4] = in[1];
	out[5] = in[5];
	out[6] = in[9];
	out[7] = in[13];
	out[8] = in[2];
	out[9] = in[6];
	out[10] = in[10];
	out[11] = in[14];
	out[12] = in[3];
	out[13] = in[7];
	out[14] = in[11];
	out[15] = in[15];
}

void m4_copy_m3(const M4 in, M4 out)
{
	// TODO: Rename better? Like specific it takes the rotation ??? 
	// or something about it overwritting the other values
	// Used for copying the top left 3x3 portion of the m4.

	
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
	out[3] = 0;

	out[4] = in[4];
	out[5] = in[5];
	out[6] = in[6];
	out[7] = 0;

	out[8] = in[8];
	out[9] = in[9];
	out[10] = in[10];
	out[11] = 0;

	out[12] = in[12];
	out[13] = in[13];
	out[14] = in[14];
	out[15] = 1; // TODO: Not sure about this
}

void m4_projection(float fov, float aspect_ratio, float near_plane, float far_plane, M4 out)
{
	// TODO: Fov is vertical fov here.
	// TODO: Comment all this properly to show I actually understand it all.

	// Currently the opengl perspective projection matrix.

	float y_scale = 1.f / tanf(radians(fov) / 2.f);
	float x_scale = y_scale / aspect_ratio;

	out[0] = x_scale;
	out[1] = 0;
	out[2] = 0;
	out[3] = 0;
	out[4] = 0;
	out[5] = y_scale;
	out[6] = 0;
	out[7] = 0;
	out[8] = 0;
	out[9] = 0;
	out[10] = -(far_plane + near_plane) / (far_plane - near_plane);
	out[11] = -1; // This negation our right handed coordinate system into a left handed coordinate system in ndc space?
	out[12] = 0;
	out[13] = 0;
	out[14] = -2 * far_plane * near_plane / (far_plane - near_plane);
	out[15] = 0;
}

char* m4_to_str(const M4 m)
{
	return format_str("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f",
		m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
}