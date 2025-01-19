#ifndef COLOUR_H
#define COLOUR_H

// Define some standard colours.
// TODO: Rename these to show they're ints.
// TODO: More
#define COLOUR_WHITE 0x00FFFFFF
#define COLOUR_BLACK 0x00000000
#define COLOUR_RED 0x00FF0000
#define COLOUR_LIME 0x0000FF00
#define COLOUR_BLUE 0x000000FF

// TODO: Rename colour helpers or something?
inline int float_rgb_to_int(float r, float g, float b)
{
	return ((int)(r * 255) << 16
		| (int)(g * 255) << 8
		| (int)(b * 255));
}

inline int int_rgb_to_int(int r, int g, int b)
{
	return (r << 16 | g << 8 | b);
}

inline void unpack_int_rgb_to_floats(int colour, float* r, float* g, float* b)
{
	float n = 1.f / 255.f;

	*r = ((colour >> 16) & 0xFF) * n;
	*g = ((colour >> 8) & 0xFF) * n;
	*b = (colour & 0xFF) * n;
}

inline void unpack_int_rgb_to_ints(int colour, int* r, int* g, int* b)
{
	*r = (colour >> 16) & 0xFF;
	*g = (colour >> 8) & 0xFF;
	*b = colour & 0xFF;
}

#endif