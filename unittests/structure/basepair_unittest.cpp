

#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/structure.h"
#include "structure/basepair.h"
#include "structure/pdb_parser.h"
#include "structure/is_equal.h"
#include "util/random_number_generator.h"


TEST_CASE ("Test base pair ") {

    auto rts = structure::ResidueTypeSet();
    auto parser = structure::PDBParser(rts);

    auto path = base::unittest_resources_path() + "/structure/p4p6.pdb";
    auto residues = parser.parse(path);
    auto s = structure::get_structure_from_residues(residues->RNA_residues);

    auto x = util::X3dna();
    auto x3dna_basepairs = x.get_basepairs(path);

    auto bps = structure::get_basepairs_from_x3dna(x3dna_basepairs, *s);

    SUBCASE("test sync between basepair and residues") {

        auto bp = bps->at(0);
        auto bp_res = structure::Residues();
        bp_res.push_back(s->get_residue(bp.get_res1_uuid()));
        bp_res.push_back(s->get_residue(bp.get_res2_uuid()));

        auto rng = util::RandomNumberGenerator();
        for (int i = 0; i < 100; i++) {
            auto rot = math::get_random_rotation_matrix();
            auto dist = util::get_random_point(rng, 10);

            bp.transform(rot, dist);
            bp_res[0].transform(rot, dist);
            bp_res[1].transform(rot, dist);
        }

// calculate center of base pair
        auto center = math::Point();
        auto count = 0;
        for (auto const &a : bp_res[0]) {
            center += a.get_coords();
            count += 1;
        }
        for (auto const &a : bp_res[1]) {
            center += a.get_coords();
            count += 1;
        }
        center /= count;

        CHECK(math::are_points_equal(bp.get_center(), center));

        CHECK(math::are_points_equal(bp.get_res1_c1_prime_coord(),
                                     bp_res[0].get_coords("C1'")));

        CHECK(math::are_points_equal(bp.get_res2_c1_prime_coord(),
                                     bp_res[1].get_coords("C1'")));

    }

}

