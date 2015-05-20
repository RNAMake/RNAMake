add_executable(design_primers ../../design_primers/main.cpp)
target_link_libraries(design_primers ${link_libraries})
add_custom_target(dp_symlink ALL)
add_custom_command(TARGET dp_symlink POST_BUILD COMMAND cmake -E create_symlink ../../bin/design_primers design_primers)
