#ifndef DRAW_2D_H
#define DRAW_2D_H

#include "core/canvas.h"

// TODO: Documentation comments
void draw_line(canvas_t* canvas, int x0, int y0, int x1, int y1, int colour);
void draw_circle(canvas_t* canvas, int cx, int cy, int r, int colour);
void draw_rect(canvas_t* canvas, int x0, int y0, int x1, int y1, int colour);

#endif