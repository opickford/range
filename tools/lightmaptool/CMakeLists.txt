# CMakeList.txt : CMake project for the range lightmaptool.
#

# Compile the engine as a static .lib to be used by the game/tools.
add_executable (lightmaptool
"main.c"

"light_mapping.c"


)

 # Link the light baking tool to the engine
target_link_libraries(lightmaptool engine)

target_include_directories (lightmaptool PRIVATE "${CMAKE_SOURCE_DIR}/engine")


# TODO: Add tests and install targets if needed.
