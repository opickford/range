# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("game" C)

# Compile options.
add_compile_options(/TC) # Compile to C.
add_compile_options(/W3) # Enable max warnings. TODO: Switch to W4.

# Include sub-projects.
add_subdirectory ("engine")
add_subdirectory ("game")

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT "game")
set_property(GLOBAL PROPERTY USE_FOLDERS ON)