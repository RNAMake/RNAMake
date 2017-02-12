
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "math/numerical.h"
#include "util/x3dna.h"
#include "structure/structure.h"
#include "structure/basepair.h"
#include "structure/is_equal.hpp"

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
    auto bp_name = primitives::calc_bp_name<Residue>(res);
    auto sugars = Points();
    sugars.push_back(res1->get_atom("C1'")->coords());
    sugars.push_back(res2->get_atom("C1'")->coords());

    //auto bp = std::make_shared<Basepair>(res1->uuid(), res2->uuid(), x_bp.r, center, sugars,
    //                                     std::make_shared<String>(bp_name),
    //                                     std::make_shared<String>(x_bp.bp_type),
    //                                     primitives::Basepair::BasepairType::WC,
    //                                     std::make_shared<Uuid>());

    /*SECTION("test whether basepair transformations reflect residue transformations") {
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

        bp_copy = std::make_shared<Basepair>(*bp, std::make_shared<Uuid>(), std::make_shared<Uuid>());
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
        bp->move(Point(10, 10, 01));
        REQUIRE(are_xyzVector_equal(bp->d(), bp_state->d()) == 0);

    }*/

}











































