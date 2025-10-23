add_library(maths STATIC
    bounding_sphere.h
    matrix4.h
    plane.h
    utils.h
    vector_maths.h
    vector2.h
    vector3.h
    vector4.h
    matrix4.c
    plane.c
)

target_include_directories(maths 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(maths 
    PUBLIC 
        utils)

set_target_properties(maths PROPERTIES FOLDER "engine")