add_library(utils STATIC
    common.h
    logger.h
    memory_utils.h
    memory_utils.c
    str_utils.h
    str_utils.c
    timer.h
)

target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(utils PROPERTIES FOLDER "Engine")