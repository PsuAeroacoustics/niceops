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
