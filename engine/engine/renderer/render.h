#ifndef RENDER_H
#define RENDER_H

#include "renderer.h"

// This stuff feels quite renderer specific, maybe except from resources.
// TODO: Think about global-stuff.
#include "scene.h"
#include "models.h"
#include "lights.h"
#include "resources.h"

#include "render_target.h"
#include "render_settings.h"
#include "canvas.h"
#include "depth_buffer.h"

#include "frustum_culling.h"

#include "maths/vector2.h"
#include "maths/vector3.h"
#include "maths/vector4.h"
#include "maths/matrix4.h"

// TODO: Organise this all.

// SECTION: Debug tools.
// TODO: Refactor and comments etc.
void debug_draw_point_lights(Canvas* canvas, const RenderSettings* settings, PointLights* point_lights);
void debug_draw_bounding_spheres(Canvas* canvas, const RenderSettings* settings, const Models* models, const M4 view_matrix);
void debug_draw_world_space_point(Canvas* canvas, const RenderSettings* settings, V3 point, const M4 view_matrix, int colour);
void debug_draw_view_space_point(Canvas* canvas, const RenderSettings* settings, V3 point, int colour);
void debug_draw_world_space_line(Canvas* canvas, const RenderSettings* settings, const M4 view_matrix, V3 v0, V3 v1, V3 colour);
void debug_draw_mi_normals(Canvas* canvas, const RenderSettings* settings, const Models* models, int mi_index);



// TODO: Split these functions into sections.

// TODO: Not sure where to put this?
float calculate_diffuse_factor(V3 v, V3 n, V3 light_pos, float a, float b);

// SECTION: Triangle rasterisation.
void draw_scanline(RenderTarget* rt, 
	RenderBuffers* rbs,
	int x0, int x1, 
	int y, 
	float z0, float z1, 
	float w0, float w1, 
	V3 ac0, V3 ac1, // Albedo
	V3 lc0, V3 lc1, // Light colour/contribution
	float* lsps, int lights_count, DepthBuffer* depth_maps);

void draw_flat_bottom_triangle(RenderTarget* rt, RenderBuffers* rbs, float* vc0, float* vc1, float* vc2, int vertex_stride, int lights_count, DepthBuffer* depth_maps);
void draw_flat_top_triangle(RenderTarget* rt, RenderBuffers* rbs, float* vc0, float* vc1, float* vc2, int vertex_stride, int lights_count, DepthBuffer* depth_maps);
void draw_triangle(RenderTarget* rt, RenderBuffers* rbs, float* vc0, float* vc1, float* vc2, float* vc3, int vertex_stride, int lights_count, DepthBuffer* depth_maps);

// TODO: Rename?
void draw_textured_scanline(RenderTarget* rt, int x0, int x1, int y, float z0, float z1, float w0, float w1, V3 c0, V3 c1, const V2 uv0, const V2 uv1, const Texture* texture);
void draw_textured_flat_bottom_triangle(RenderTarget* rt, V4 v0, V4 v1, V4 v2, V3 c0, V3 c1, V3 c2, V2 uv0, V2 uv1, V2 uv2, const Texture* texture);
void draw_textured_flat_top_triangle(RenderTarget* rt, V4 v0, V4 v1, V4 v2, V3 c0, V3 c1, V3 c2, V2 uv0, V2 uv1, V2 uv2, const Texture* texture);
void draw_textured_triangle(RenderTarget* rt, V4 v0, V4 v1, V4 v2, V3 c0, V3 c1, V3 c2, V2 uv0, V2 uv1, V2 uv2, const Texture* texture);

// TODO: Comments etc.
void draw_depth_scanline(DepthBuffer* db, int x0, int x1, int y, float z0, float z1);
void draw_depth_flat_bottom_triangle(DepthBuffer* db, V4 v0, V4 v1, V4 v2);
void draw_depth_flat_top_triangle(DepthBuffer* db, V4 v0, V4 v1, V4 v2);
void draw_depth_triangle(DepthBuffer* db, V4 v0, V4 v1, V4 v2);

// SECTION: Rendering Pipeline.
void project(const Canvas* canvas, const M4 projection_matrix, V4 v, V4* out);

void model_to_view_space(Models* models, const M4 view_matrix);

void lights_world_to_view_space(PointLights* point_lights, const M4 view_matrix);

void broad_phase_frustum_culling(Models* models, const ViewFrustum* view_frustum);

void cull_backfaces(Renderer* renderer, const Scene* scene);

void light_front_faces(Renderer* renderer, Scene* scene);

void clip_to_screen(Renderer* renderer, const M4 view_matrix, Scene* scene, const Resources* resources);

void project_and_draw_clipped(Renderer* renderer, Scene* scene, int mi_index, int clipped_face_count, const Resources* resources);

void render(Renderer* renderer, Scene* scene, const Resources* resources, const M4 view_matrix);

// REFACTORED PIPELINE

void calculate_model_matrices();




// TEMP


void update_depth_maps(Renderer* renderer, const Scene* scene);


#endif