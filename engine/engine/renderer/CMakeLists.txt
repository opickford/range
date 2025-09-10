add_library(renderer STATIC
    camera.h
    depth_buffer.h
    draw_2d.h
    frame_data.h
    frustum_culling.h
    render_settings.h
    render_target.h
    render.h
    renderer.h
    texture.h

    camera.c
    frustum_culling.c
    render.c
    renderer.c
    draw_2d.c
    depth_buffer.c
    texture.c
    frame_data.c
)

target_link_libraries(renderer 
    PUBLIC
        cecs::cecs core maths utils)

target_include_directories(renderer 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(renderer PROPERTIES FOLDER "Engine")