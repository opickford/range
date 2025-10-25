add_library(rendering STATIC
    billboard.h
    billboard.c
)

target_include_directories(rendering
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(rendering  
    PRIVATE 
        utils maths common cecs::cecs chds::chds
)

set_target_properties(rendering PROPERTIES FOLDER "engine/scripts")
