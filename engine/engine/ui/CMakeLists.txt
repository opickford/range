add_library(ui STATIC
    font.c
    font.h
    text.c
    text.h
    ui.c
    ui.h
)
target_link_libraries(ui 
    PUBLIC
        core
)

target_include_directories(ui 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)


target_include_directories(ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


set_target_properties(ui PROPERTIES FOLDER "engine_t")