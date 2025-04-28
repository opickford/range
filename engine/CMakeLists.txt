# CMakeList.txt : CMake project for the range game engine.
#

# Compile the engine as a static .lib to be used by the game/tools.
add_library (engine STATIC

"engine/canvas.c"
"engine/engine.c"
"engine/lights.c"
 
"engine/window.c"
"engine/scene.c"
"engine/strides.c"

"engine/mesh_base.c" 
"engine/mesh_instance.c" 

"engine/maths/matrix4.c"
"engine/maths/plane.c"

"engine/renderer/camera.c"
"engine/renderer/frustum_culling.c"
"engine/renderer/render.c"
"engine/renderer/renderer.c"
"engine/renderer/draw_2d.c"
"engine/renderer/depth_buffer.c"
"engine/renderer/texture.c"
"engine/renderer/frame_data.c"

"engine/ui/font.c"
"engine/ui/text.c"
"engine/ui/ui.c"

"engine/utils/memory_utils.c"
"engine/utils/str_utils.c"
  )

target_include_directories (engine PUBLIC "${CMAKE_SOURCE_DIR}/engine/engine")

# TODO: Add tests and install targets if needed.
