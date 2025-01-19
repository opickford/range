#ifndef DRAW_2D_H
#define DRAW_2D_H

#include "canvas.h"

// TODO: Documentation comments
void draw_line(Canvas* canvas, int x0, int y0, int x1, int y1, int colour);
void draw_circle(Canvas* canvas, int cx, int cy, int r, int colour);
void draw_rect(Canvas* canvas, int x0, int y0, int x1, int y1, int colour);

#endif