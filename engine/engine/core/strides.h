#ifndef STRIDES_H
#define STRIDES_H

/*
TODO: Toplevel comments.

Defines the strides for models/lights/everything?


Here until I can think of a better place to put it.
// TODO: Should this go into common? or renderer? 
// TODO: File structure needs a bit of a restrucutre.
*/


// TODO: Are any of these unnecessary? Gotta think about it.
//		 A position could be 3 or 4, a colour could be 3 or 4 it's 
//		 a tricky one. Will have to take a look into the renderer basically.

#define STRIDE_FACE_VERTICES 3

#define STRIDE_POSITION 3 // TODO: Just remove this? STRIDE_V3 and STRIDE_V4 so it's not just 3s and 4s.
#define STRIDE_V4		4 // TODO: ^ Maybe temp.
#define STRIDE_COLOUR	3

#define STRIDE_NORMAL	3
#define STRIDE_UV		2

#define STRIDE_SPHERE	4				// Center (x,y,z), Radius	
#define STRIDE_POINT_LIGHT_ATTRIBUTES 4 // r,g,b,strength
#define STRIDE_MI_TRANSFORM 9			// Position, Eulers, Scale

// TODO: Not sure on the ENTIRE naming conventions. Could make this better.

/*
#define STRIDE_ENTIRE_VERTEX 11 // x, y, z, u, v, nx, ny, nz, r, g, b - TODO: Could think about this a bit more.


#define STRIDE_ENTIRE_FACE (STRIDE_ENTIRE_VERTEX * STRIDE_FACE_VERTICES)
*/
// TODO: Avoid copying uv when we don't need it? Might be difficult.

// New strides.


// TODO: Some of these strides are variable??? How do we handle this....

// Front faces go into clipping. they should have all information there ready.
// because lighting is performed before.


// These only define the constant part of the vertex?? obviously variable with
// if uv or shadow is needed.




/*

Plan is. front faces will read the light space positions and copy them 
to the buffer to clip.


*/

// Again these are all the static ones.

#define STRIDE_VISIBLE_VERTEX 14

// Front faces need a pos (V3), UV (V2), normal (V3), albedo (V3), diffuse (V3)
// as well as the space for light space coordinates.......
#define STRIDE_BASE_FRONT_VERTEX 14
#define STRIDE_BASE_FRONT_FACE (STRIDE_BASE_FRONT_VERTEX * STRIDE_FACE_VERTICES)

// TODO: We do NOT want the normal here. How do we remove it? Unless we write out 
//		 the data again after lighting.... but then the backface culling should not
//		 copy all the data. It would only need vertex pos, vertex normal and vertex albedo.
//		 it could then write out everything except the normal.....

// Same as a front face for now.
#define STRIDE_BASE_CLIPPED_VERTEX 14
#define STRIDE_BASE_CLIPPED_FACE (STRIDE_BASE_CLIPPED_VERTEX * STRIDE_FACE_VERTICES)

#endif