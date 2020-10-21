

set(x3dna_files
        ${RNAMAKE}/src/external/x3dna/alc2img.c
        ${RNAMAKE}/src/external/x3dna/ana_fncs.c
        ${RNAMAKE}/src/external/x3dna/analyze.c
        ${RNAMAKE}/src/external/x3dna/anyhelix.c
        ${RNAMAKE}/src/external/x3dna/app_fncs.c
        ${RNAMAKE}/src/external/x3dna/cehs.c
        ${RNAMAKE}/src/external/x3dna/changelog.md
        ${RNAMAKE}/src/external/x3dna/cmn_fncs.c
        ${RNAMAKE}/src/external/x3dna/comb_str.c
        ${RNAMAKE}/src/external/x3dna/ex_str.c
        ${RNAMAKE}/src/external/x3dna/fiber.c
        ${RNAMAKE}/src/external/x3dna/file_list.rb
        ${RNAMAKE}/src/external/x3dna/find_pair.c
        ${RNAMAKE}/src/external/x3dna/find_platform.c
        ${RNAMAKE}/src/external/x3dna/fncs_slre.c
        ${RNAMAKE}/src/external/x3dna/frame_mol.c
        ${RNAMAKE}/src/external/x3dna/get_part.c
        ${RNAMAKE}/src/external/x3dna/mutate_bases.c
        ${RNAMAKE}/src/external/x3dna/nrutil.c
        ${RNAMAKE}/src/external/x3dna/o1p_o2p.c
        ${RNAMAKE}/src/external/x3dna/pdb2img.c
        ${RNAMAKE}/src/external/x3dna/r3d_atom.c
        ${RNAMAKE}/src/external/x3dna/reb_fncs.c
        ${RNAMAKE}/src/external/x3dna/rebuild.c
        ${RNAMAKE}/src/external/x3dna/regular_dna.c
        ${RNAMAKE}/src/external/x3dna/rotate_mol.c
        ${RNAMAKE}/src/external/x3dna/stack2img.c
        ${RNAMAKE}/src/external/x3dna/std_base.c
        ${RNAMAKE}/src/external/x3dna/step_hel.c
        ${RNAMAKE}/src/external/x3dna/x3dna_fncs.h
        ${RNAMAKE}/src/external/x3dna/x3dna.h
        ${RNAMAKE}/src/external/x3dna/functions.h
        ${RNAMAKE}/src/external/x3dna/x3dna_wrapper.cpp


        )

add_library(x3dna ${x3dna_files})

target_link_libraries(x3dna)