add_library(physics STATIC
    physics.h
    physics.c
)

# TODO: This is the same for all libs right?
target_include_directories(physics 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(physics 
    PUBLIC 
        cecs::cecs core maths utils)

set_target_properties(physics PROPERTIES FOLDER "Engine")