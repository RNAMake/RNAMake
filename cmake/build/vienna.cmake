
include_directories(${RNAMAKE}/src/external/ViennaRNA-2.1.9/lib/)
include_directories(${RNAMAKE}/src/external/ViennaRNA-2.1.9/H/)
include_directories(${RNAMAKE}/src/external/ViennaRNA-2.1.9/)
include_directories(${RNAMAKE}/src/external/ViennaRNA-2.1.9/libsvm-2.91/)

set(vienna_files
    ${RNAMAKE}/src/external/ViennaRNA-2.1.9/lib/energy_par.c
    ${RNAMAKE}/src/external/ViennaRNA-2.1.9/lib/fold.c
    ${RNAMAKE}/src/external/ViennaRNA-2.1.9/lib/fold_vars.c
    ${RNAMAKE}/src/external/ViennaRNA-2.1.9/lib/gquad.c
    ${RNAMAKE}/src/external/ViennaRNA-2.1.9/lib/params.c
    ${RNAMAKE}/src/external/ViennaRNA-2.1.9/lib/utils.c
)

set_source_files_properties( ${vienna_files}
        PROPERTIES LANGUAGE C
        )

add_library(vienna_RNA ${vienna_files})



target_link_libraries( vienna_RNA )