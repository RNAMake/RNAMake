
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "math/numerical.h"
#include "util/x3dna.h"
#include "structure/structure.h"
#include "structure/basepair.h"

state::BasepairOP
_get_bp_from_str(String const & s) {
    auto spl = split_str_by_delimiter(s, ";");
    auto d = vector_from_str(spl[0]);
    auto r = matrix_from_str(spl[1]);
    auto sugars = vectors_from_str(spl[2]);
    auto bp_name_str = String("test");
    auto bp_name = Chars(4);
    int i = 0;
    for(auto const & e : bp_name_str) { bp_name[i] = e; i++; }
    return std::make_shared<state::Basepair>(Uuid(), Uuid(), r, d, sugars[0], sugars[1], bp_name,
                                             X3dna::X3dnaBPType::cDDD,
                                             primitives::Basepair::BasepairType::WC, Uuid());
}

TEST_CASE( "Test Basepairs for Structure", "[Basepair]" ) {
    auto path = py_unittest_path() + "/resources/motifs/p4p6/p4p6.pdb";
    auto rts = ResidueTypeSet();
    auto s = structure_from_pdb(path, rts);
    auto x = X3dna();
    auto x_bps = x.get_basepairs(path);
    auto x_bp = x_bps[0];

    auto res1 = s->get_residue(x_bp.res1.num, x_bp.res1.chain_id, x_bp.res1.i_code);
    auto res2 = s->get_residue(x_bp.res2.num, x_bp.res2.chain_id, x_bp.res2.i_code);
    auto res = ResidueOPs { res1, res2};
    auto center = _calc_center(res);
    auto bp_name_str = primitives::calc_bp_name<Residue>(res);
    auto sugars = Points();
    sugars.push_back(res1->get_atom("C1'")->coords());
    sugars.push_back(res2->get_atom("C1'")->coords());
    auto bp_name = Chars(bp_name_str.length());
    int i = 0;
    for(auto const & e : bp_name_str) { bp_name[i] = e; i++; }

    auto bp = std::make_shared<Basepair>(res1->uuid(), res2->uuid(), x_bp.r, center, sugars,
                                         bp_name, x_bp.bp_type,
                                         primitives::Basepair::BasepairType::WC, Uuid());

    SECTION("test whether basepair transformations reflect residue transformations") {
        auto path = unittest_resource_dir() + "/math/random_transformations.dat";
        auto lines = get_lines_from_file(path);
        for(auto const & l : lines) {
            auto spl = split_str_by_delimiter(l, "|");
            if(spl.size() < 2) { break; }
            auto r = matrix_from_str(spl[0]);
            auto trans = vector_to_str(spl[1]);

            auto t = Transform(r, trans);
            bp->transform(t);

            for(auto & r : res) { r->transform(t); }

        }

        REQUIRE(are_xyzVector_equal(bp->d(), _calc_center(res)));
    }

    SECTION("test copying") {
        auto bp_copy = std::make_shared<Basepair>(*bp);
        REQUIRE(are_basepairs_equal(bp, bp_copy));

        bp_copy = std::make_shared<Basepair>(*bp, Uuid(), Uuid(), Uuid());
        REQUIRE(are_basepairs_equal(bp, bp_copy) == 0);
        REQUIRE(are_basepairs_equal(bp, bp_copy, 0));
    }

    SECTION("test to str") {
        auto s = bp->to_str();
        auto bp_copy = std::make_shared<Basepair>(s, bp->res1_uuid(), bp->res2_uuid());

        REQUIRE(are_basepairs_equal(bp, bp_copy, 0));
    }

    SECTION("test get state") {
        auto bp_state = bp->get_state();
        bp->move(Point(10, 10, 10));
        REQUIRE(are_xyzVector_equal(bp->d(), bp_state->d()) == 0);

        auto bp_state_2 = std::make_shared<state::Basepair>(*bp_state);
        REQUIRE(state::are_basepairs_equal(bp_state, bp_state_2));

        bp_state_2 = std::make_shared<state::Basepair>(bp_state->to_str(), bp_state->res1_uuid(),
                                                       bp_state->res2_uuid(), bp->uuid());
        REQUIRE(state::are_basepairs_equal(bp_state, bp_state_2));


    }

    SECTION("test get_transforming_r_and_t") {
        auto path = py_unittest_path() + "resources/motif_state/get_transforming_r_and_t_test.dat";
        auto lines = get_lines_from_file(path);
        for(auto const & l : lines) {
            if(l.length() < 10) { break; }
            auto spl = split_str_by_delimiter(l, "|");
            auto bp1 = _get_bp_from_str(spl[0]);
            auto bp2 = _get_bp_from_str(spl[1]);
            auto t = vector_from_str(spl[2]);
            auto r = matrix_from_str(spl[3]);
            auto transform_info = state::TransformInfo();

            bp1->get_transforming_r_and_t(*bp2, transform_info);
            transform_info.r.transpose();
            auto result_t = t + bp1->d();
            REQUIRE(are_xyzMatrix_equal(r, transform_info.r));
            REQUIRE(are_xyzVector_equal(result_t, transform_info.d));

        }
    }

    SECTION("test random tranformation between basepairs") {
        auto path = unittest_resource_dir() + "/structure/test_get_transformed_state.dat";
        auto lines = get_lines_from_file(path);
        auto transform_info = state::TransformInfo();

        for(auto const & l : lines) {
            if(l.length() < 5) { break;}
            auto spl = split_str_by_delimiter(l, "|");
            auto bp1 = _get_bp_from_str(spl[0]);
            auto bp2 = _get_bp_from_str(spl[1]);

            bp1->get_transforming_r_and_t(*bp2, transform_info);
            bp2->fast_transform(transform_info);

            REQUIRE(are_xyzMatrix_equal(bp1->r(), bp2->r()));
            REQUIRE(are_xyzVector_equal(bp1->d(), bp2->d()));
            //REQUIRE(are_xyzVector_equal(bp1->sugars()[0], bp2->sugars()[0]));
            //REQUIRE(are_xyzVector_equal(bp1->sugars()[1], bp2->sugars()[1]));

        }


    }

}











































