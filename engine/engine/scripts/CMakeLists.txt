add_subdirectory(gameplay)
add_subdirectory(rendering)

add_library(scripts INTERFACE)

# target_include_directories(scripts 
#     PUBLIC
#         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
#         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
#         $<INSTALL_INTERFACE:include>
# )

target_link_libraries(scripts
    INTERFACE
        gameplay rendering
)

set_target_properties(scripts PROPERTIES FOLDER "engine/scripts")
