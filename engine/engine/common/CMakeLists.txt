add_library(common INTERFACE
    colour.h
    status.h
)

target_include_directories(common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(common PROPERTIES FOLDER "Engine")