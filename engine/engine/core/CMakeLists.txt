# CMakeList.txt : CMake project for the range game engine.
#




# Compile the engine as a static .lib to be used by the game/tools.
add_library (core STATIC
    canvas.c
    engine.c
    lights.c
    window.c
    scene.c
    strides.c
    mesh_base.c
    mesh_instance.c
    canvas.h
    engine.h
    globals.h
    lights.h
    mesh_base.h
    mesh_instance.h
    resources.h
    scene.h
    strides.h
    window.h
	components.h
)

target_link_libraries(core 
    PUBLIC
        cecs::cecs common maths utils
)

target_include_directories(core 
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(core PROPERTIES FOLDER "Engine")

# TODO: Add tests and install targets if needed.

