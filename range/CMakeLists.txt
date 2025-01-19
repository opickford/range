# CMakeList.txt : CMake project for the range game.
#

# Compile the engine as a static .lib to be used by the game/tools.
add_executable (range
	"main.c"
 )

 # Link the light baking tool to the engine
target_link_libraries(range engine)

target_include_directories (range PRIVATE "${CMAKE_SOURCE_DIR}/engine")

# TODO: Add tests and install targets if needed.
