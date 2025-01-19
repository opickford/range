#include "draw_2d.h"

#include <math.h>

void draw_line(Canvas* canvas, int x0, int y0, int x1, int y1, int colour)
{
	int dx = abs(x1 - x0);
	int sx = x0 < x1 ? 1 : -1;
	int dy = -abs(y1 - y0);
	int sy = y0 < y1 ? 1 : -1;
	int error = dx + dy;

	while (1)
	{
		if (x0 > -1 && x0 < canvas->width - 1 && y0 > -1 && y0 < canvas->height - 1)
		{
			int pos = y0 * canvas->width + x0;
			canvas->pixels[pos] = colour;
		}

		if (x0 == x1 && y0 == y1) break;

		int e2 = 2 * error;
		if (e2 >= dy)
		{
			error = error + dy;
			x0 = x0 + sx;
		}

		if (e2 <= dx)
		{
			error = error + dx;
			y0 = y0 + sy;
		}
	}
}

void draw_circle(Canvas* canvas, int cx, int cy, int r, int colour)
{

	// (x-a)^2 + (y-b)^2

	// All this could be better, doesn't work great.
	int rr = r * r;

	for (int y = -r; y < r; ++y)
	{
		if (cy + y < 0 || cy + y > canvas->height - 1)
		{
			continue;
		}

		int xx = abs(y * y - rr);
		float x = sqrtf((float)xx);

		int x0 = (int)(cx - x);
		int x1 = (int)(cx + x);

		if (x0 > -1 && x0 < canvas->width)
		{
			canvas->pixels[(int)((cy + y) * canvas->width + x0)] = colour;
		}

		if (x1 > -1 && x1 < canvas->width)
		{
			canvas->pixels[(int)((cy + y) * canvas->width + x1)] = colour;
		}
	}
}

void draw_rect(Canvas* canvas, int x0, int y0, int x1, int y1, int colour)
{
	// TODO: Can optimise.
	for (int y = y0; y < y1; ++y)
	{
		if (y < 0 || y >= canvas->height)
		{
			continue;
		}

		for (int x = x0; x < x1; ++x)
		{

			if (x < 0 || x >= canvas->width)
			{
				continue;
			}

			int i = y * canvas->width + x;
			canvas->pixels[i] = colour;
		}
	}
}