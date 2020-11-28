#ifndef RNAMAKE_SS_MOTIF_UNITTESTS
#define RNAMAKE_SS_MOTIF_UNITTESTS

#include <algorithm>

#include "../common.hpp"

#include <base/types.h>

#include <rnamake2d/ss_motif.h>
#include <secondary_structure/secondary_structure_parser.h>

TEST_CASE("Test rnamake2d::Motif() Classes", "[rnamake2d::Motif()]") {
    SECTION("Small motif test case") {
        auto parser = secondary_structure::Parser();
        const auto seq = "GGGGAAAACCCC";
        const auto db = "((((....))))";
        const auto pose = parser.parse_to_pose(seq, db);
        auto [hps, helices, jncs, ss] = rnamake2d::parse_to_motif2ds(pose, db, seq);

        auto strands1 = std::vector<Ints>{Ints{0,1,2,3}, Ints{8,9,10,11}};
        auto strands2 = std::vector<Ints>{Ints{3,4,5,6,7,8}};
        auto test_hel = rnamake2d::Helix(strands1, seq);
        auto test_hp = rnamake2d::Hairpin(strands2, seq);
        test_hp.buffer(4);

        // checking proper output of helix vectors
        REQUIRE(hps.size() == 1);
        REQUIRE(helices.size() == 1);
        REQUIRE(jncs.empty());
        REQUIRE(ss.empty());
        // check specifics of the motifs
        REQUIRE(hps[0]->size_ == 4);
        // check the motifs are the same
        REQUIRE(*helices[0] == test_hel);
        REQUIRE(*(hps[0]) == test_hp);
    }

    SECTION("Testing Single Strand Edge Case") {
        auto parser = secondary_structure::Parser();
        const auto seq = "AAAAGGGGAAAACCCCAAA";
        const auto db = "....((((....))))...";
        const auto pose = parser.parse_to_pose(seq, db);
        auto [hps, helices, jncs, ss] = rnamake2d::parse_to_motif2ds(pose, db, seq);

        auto strands1 = std::vector<Ints>{Ints{0,1,2,3}};
        auto strands2 = std::vector<Ints>{Ints{16,17,18}};

        auto test_ss1 = rnamake2d::SingleStrand(strands1, seq);
        auto test_ss2 = rnamake2d::SingleStrand(strands2 ,seq);

        std::sort(ss.begin(), ss.end(), [] (const Motif2DOP& m1, const Motif2DOP& m2) {
            return m1->sequence.size() < m2->sequence.size();
        }) ;

        REQUIRE(*ss[0] == test_ss2);
        REQUIRE(*ss[1] == test_ss1);

    }

    SECTION("Medium size structure with a variety of Motifs") {

    }
}



#endif // RNAMAKE_SS_MOTIF_UNITTESTS
