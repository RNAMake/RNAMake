//
// Created by Joseph Yesselman on 4/4/18.
//

#ifndef TEST_CLOSE_CHAIN_H
#define TEST_CLOSE_CHAIN_H

#include <math/matrix_3x3.hpp>
#include <structure/chain.h>

namespace structure {

math::Matrix3x3 create_coord_system(AtomOPs const &);

float to_radians(float);

AtomOP virtual_atom(String const &, float, float, float, AtomOPs const &);

math::Vector3 get_projection(math::Vector3 const &, math::Vector3 const &,
                            math::Vector3 const &);

math::Matrix3x3 axis_angle_to_rot_matrix(float, math::Vector3 const &);

void close_torsion(int, AtomOPs const &, AtomOPs, AtomOPs const &,
                   AtomOPs const &);

math::Matrix3x3 get_res_ref_frame(ResidueOP);

void replace_missing_phosphate_backbone(ResidueOP, ResidueOP);

void close_chain(ChainOP);

} // namespace structure

#endif // TEST_CLOSE_CHAIN_H
