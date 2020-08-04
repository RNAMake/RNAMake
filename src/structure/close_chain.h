//
// Created by Joseph Yesselman on 4/4/18.
//

#ifndef TEST_CLOSE_CHAIN_H
#define TEST_CLOSE_CHAIN_H

#include <math/xyz_matrix.h>
#include <structure/chain.h>

namespace structure {

math::Matrix
create_coord_system(
        AtomOPs const &);

float
to_radians(
        float);

AtomOP
virtual_atom(
        String const &,
        float,
        float,
        float,
        AtomOPs const &);

math::Vector
get_projection(
        math::Point const &,
        math::Point const &,
        math::Vector const &);

math::Matrix
axis_angle_to_rot_matrix(
        float,
        math::Vector const &);

void
close_torsion(
        int,
        AtomOPs const &,
        AtomOPs,
        AtomOPs const &,
        AtomOPs const &);

math::Matrix
get_res_ref_frame(
        ResidueOP);

void
replace_missing_phosphate_backbone(
        ResidueOP,
        ResidueOP);


void
close_chain(
        ChainOP);

}


#endif //TEST_CLOSE_CHAIN_H
