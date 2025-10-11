#ifndef RENDER_H
#define RENDER_H

#include "renderer.h"

// This stuff feels quite renderer specific, maybe except from resources.
// TODO: Think about global-stuff.
#include "frame_data.h"
#include "render_target.h"
#include "render_settings.h"
#include "depth_buffer.h"
#include "frustum_culling.h"

#include "core/scene.h"
#include "core/mesh_base.h"
#include "core/mesh_instance.h"
#include "core/lights.h"
#include "core/resources.h"
#include "core/canvas.h"

#include "maths/vector2.h"
#include "maths/vector3.h"
#include "maths/vector4.h"
#include "maths/matrix4.h"

#include <cecs/ecs.h>

// TODO: Organise this all.

// SECTION: Debug tools.
// TODO: Refactor and comments etc.

void debug_draw_point_lights(
    canvas_t* canvas,
    const cecs_t* ecs,
    const cecs_view_id_t lighting_view,
    const frame_data_t* frame_data,
    const render_settings_t* settings
);

void debug_draw_view_space_point(canvas_t* canvas, const render_settings_t* settings, v3_t point, int colour);
void debug_draw_normals(canvas_t* canvas, const frame_data_t* frame_data, const render_settings_t* settings, const scene_t* scene);

/*void debug_draw_point_lights(canvas_t* canvas, const render_settings_t* settings, PointLights* point_lights);
void debug_draw_bounding_spheres(canvas_t* canvas, const render_settings_t* settings, const Models* models, const m4_t view_matrix);
*/
void debug_draw_world_space_point(canvas_t* canvas, const render_settings_t* settings, v3_t point, const m4_t view_matrix, int colour);
void debug_draw_world_space_line(canvas_t* canvas, const render_settings_t* settings, const m4_t view_matrix, v3_t v0, v3_t v1, v3_t colour);



// TODO: Split these functions into sections.

// TODO: Not sure where to put this?
float calculate_diffuse_factor(v3_t v, v3_t n, v3_t light_pos, float a, float b);

// SECTION: Triangle rasterisation.
void draw_scanline(render_target_t* rt, 
	int x0, int x1, 
	int y, 
	float z0, float z1, 
	float w0, float w1, 
	v3_t ac0, v3_t ac1);

void draw_flat_bottom_triangle(render_target_t* rt, float* vc0, float* vc1, float* vc2);
void draw_flat_top_triangle(render_target_t* rt, float* vc0, float* vc1, float* vc2);
void draw_triangle(render_target_t* rt, float* vc0, float* vc1, float* vc2);

// TODO: Rename?
void draw_textured_scanline(render_target_t* rt, int x0, int x1, int y, float z0, float z1, float w0, float w1, v3_t c0, v3_t c1, const v2_t uv0, const v2_t uv1, const texture_t* texture);
void draw_textured_flat_bottom_triangle(render_target_t* rt, float* vc0, float* vc1, float* vc2, const texture_t* texture);
void draw_textured_flat_top_triangle(render_target_t* rt, float* vc0, float* vc1, float* vc2, const texture_t* texture);
void draw_textured_triangle(render_target_t* rt, float* vc0, float* vc1, float* vc2, const texture_t* texture);

// TODO: Comments etc.
void draw_depth_scanline(depth_buffer_t* db, int x0, int x1, int y, float z0, float z1);
void draw_depth_flat_bottom_triangle(depth_buffer_t* db, v4_t v0, v4_t v1, v4_t v2);
void draw_depth_flat_top_triangle(depth_buffer_t* db, v4_t v0, v4_t v1, v4_t v2);
void draw_depth_triangle(depth_buffer_t* db, v4_t v0, v4_t v1, v4_t v2);

// SECTION: Rendering Pipeline.
void project(const canvas_t* canvas, const m4_t projection_matrix, v4_t v, v4_t* out);




// REFACTORED PIPELINE

void model_to_view_space(cecs_t* ecs, cecs_view_id_t render_view, frame_data_t* frame_data, scene_t* scene, const m4_t view_matrix);
void lights_world_to_view_space(cecs_t* ecs, cecs_view_id_t lighting_view, frame_data_t* frame_data, const scene_t* scene, const m4_t view_matrix);
void broad_phase_frustum_culling(cecs_t* ecs, cecs_view_id_t render_view, frame_data_t* frame_data, scene_t* scene, const view_frustum_t* view_frustum);
void cull_backfaces(cecs_t* ecs, cecs_view_id_t render_view, frame_data_t* frame_data, scene_t* scene);

void light_front_faces(
    cecs_t* ecs, 
    cecs_view_id_t render_view, 
    cecs_view_id_t lighting_view, 
    frame_data_t* frame_data, 
    scene_t* scene, 
    const v3_t ambient
);

void prepare_for_clipping(cecs_t* ecs, cecs_view_id_t render_view, frame_data_t* frame_data, scene_t* scene);
void clip_project_and_draw(
    cecs_view_id_t render_view,
	renderer_t* renderer,
	render_target_t* rt,
	frame_data_t* frame_data,
	scene_t* scene,
    const resources_t* resources);

// TODO: Hate this name.
void project_and_draw_clipped(
    renderer_t* renderer,
    render_target_t* rt,
    frame_data_t* frame_data,
    const mesh_instance_t* mi);

void project_and_draw_clipped_textured(
    renderer_t* renderer,
    render_target_t* rt,
    frame_data_t* frame_data,
    const mesh_instance_t* mi,
    const texture_t* texture);

void render(
    cecs_t* ecs,
    cecs_view_id_t render_view,
    cecs_view_id_t lighting_view,
    renderer_t* renderer,
    scene_t* scene,
    const resources_t* resources,
    const m4_t view_matrix);

// TEMP


void update_depth_maps(renderer_t* renderer, const scene_t* scene);


#endif