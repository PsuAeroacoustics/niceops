include_directories(source)
add_library(matplotlib-d
        source/matplotlibd/core/pycall.d
        source/matplotlibd/core/translate.d
        source/matplotlibd/pyplot.d
        )
target_link_libraries(matplotlib-d)

set(view_dir ${CMAKE_CURRENT_SOURCE_DIR}/views)
target_compile_options(matplotlib-d PRIVATE "-J ${view_dir}")
set_target_properties(matplotlib-d PROPERTIES TEXT_INCLUDE_DIRECTORIES view_dir)
set_target_properties(
    matplotlib-d PROPERTIES
    #RUNTIME_OUTPUT_DIRECTORY ./
    #LIBRARY_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_NAME matplotlib-d
    #LIBRARY_OUTPUT_NAME matplotlib-d
    #RUNTIME_OUTPUT_NAME matplotlib-d
    #COMPILE_FLAGS ""
    LINKER_LANGUAGE CXX
    #LINK_FLAGS ""
    # Use a custom .props file to set up Visual D (import paths, predefined versions...).
    #VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
)
