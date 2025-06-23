# CMakeList.txt : CMake project for the range game engine.
#

# TODO: Fix comments.

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/engine/libs/cecs")
find_package(cecs REQUIRED)

# Compile the engine as a static .lib to be used by the game/tools.
add_subdirectory(engine/common)
add_subdirectory(engine/utils)
add_subdirectory(engine/maths)
add_subdirectory(engine/core)
add_subdirectory(engine/renderer)
add_subdirectory(engine/ui)

add_library(engine INTERFACE)


target_link_libraries(engine 
    INTERFACE 
        core common maths renderer ui utils)

set_target_properties(engine PROPERTIES FOLDER "Engine")

get_target_property(EXISTS cecs::cecs NAME)
message(STATUS "Target exists: ${EXISTS}")


get_target_property(ENGINE_INCLUDES engine INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "engine includes: ${CECS_INCLUDES}")

get_target_property(CECS_INCLUDES cecs::cecs INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "cecs::cecs includes: ${CECS_INCLUDES}")

# TODO: Add tests and install targets if needed.

# Install libraries
# Tell CMake where to find cecsConfig.cmake

