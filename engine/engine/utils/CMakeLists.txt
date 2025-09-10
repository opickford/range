add_library(utils STATIC
    common.h
    logger.h
    str_utils.h
    str_utils.c
    timer.h
	vector.h
	vector.c
)

target_include_directories(utils 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(utils 
    PUBLIC 
        common)

set_target_properties(utils PROPERTIES FOLDER "Engine")