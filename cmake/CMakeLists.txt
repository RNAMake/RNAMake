####################################################################################################
# Project Level Info
####################################################################################################
cmake_minimum_required(VERSION 3.0)
project(RNAMake)

set(CMAKE_BUILD_TYPE Release)
set(RNAMAKE /Users/cjurich/projects/RNAMake/)

if (NOT RNAMAKE)
    message(FATAL_ERROR "Variable \"RNAMAKE\" is not set in the core CMakeLists.txt file. Please set it. Exiting")
endif()

include(compiler.cmake)
include(sqlite.cmake)

include(libs.cmake)

include_directories(${RNAMAKE}/src/)
include_directories(${RNAMAKE}/src/external/)
include_directories(${RNAMAKE}/apps/)

if(TESTS)
    include_directories(${RNAMAKE}/unittests/)
    include(unittests.cmake)
endif()
# library imports are found in here

# deciding which apps will be made
foreach(app_name
        ddg_tecto;
        #structure_info;
        #design_rna;
        #simulate_tectos;
        simulate_tectos_devel;
        sequence_optimizer_app;
        #sample_helix;
        thermo_simulation;
        #general_helix_sampler;
        #opt_tecto_cutoff;
        #apt_stablization;
        #dock_motif;
        #apt_new_interface;
        #build_flex_helices:
        sequence_optimization_benchmarks;
        #hash_motifs:
        get_best_solutions;
        #mini_ttr;
        #exhustive_eternabot;
        eternabot;
        #follow_path;
        #follow_path_process;
        #path_builder;
        #path_builder_new;
        #apt_new_interface;
        #apt_stablization;
        #build_flex_helices;
        #ddg_tecto;
        #design_rna:
        design_rna_scaffold;
        #dock_motif;
        #eternabot;
        #exhustive_eternabot;
        #follow_path;
        #follow_path_process;
        #general_helix_sampler;
        #get_best_solutions;
        #hash_motifs;
        #mini_ttr;
        #opt_tecto_cutoff;
        #p5abc;
        #path_builder;
        #path_builder_new;
        #sample_helix;
        #sequence_optimization_benchmarks;
        #sequence_optimizer_app;
        #simulate_tectos:
        #simulate_tectos_devel:
        #structure_info:
        #thermo_simulation;
        )
    add_executable(${app_name} ${RNAMAKE}/apps/${app_name}/${app_name}.cc)
    target_link_libraries(${app_name} all_lib)
endforeach()


