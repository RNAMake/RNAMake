add_executable(simulate_tectos ../../apps/simulate_tectos/simulate_tectos.cc)
target_link_libraries(simulate_tectos ${link_libraries})
add_custom_target(st_symlink ALL)
add_custom_command(TARGET st_symlink POST_BUILD COMMAND cmake -E create_symlink ../../bin/simulate_tectos simulate_tectos)
