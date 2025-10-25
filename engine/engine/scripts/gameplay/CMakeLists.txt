add_library(gameplay STATIC
    player_controller.h
    player_controller.c
)

target_include_directories(gameplay
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(gameplay  
    PRIVATE 
        utils maths common cecs::cecs chds::chds
)

set_target_properties(gameplay PROPERTIES FOLDER "engine/scripts")
