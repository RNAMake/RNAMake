

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include <math/quaternion.h>
#include "motif/motif_state.h"
#include "motif/motif.h"
#include <motif/motif_state_aligner.h>

math::Matrix
rotation_about_x_axis(
        float degrees) {
    float a = 3.14/180*degrees;
    return math::Matrix(1.0, 0.0, 0.0,
                  0.0, cos(a), sin(a),
                  0.0, -sin(a), cos(a));
}

TEST_CASE( "Test Motif states, motifs that dont have coordinates" ) {
    auto path = base::motif_dirs() + "base.motif";
    auto m = motif::file_to_motif(path);
    path = base::motif_dirs() + "ref.motif";
    auto ref_m = motif::file_to_motif(path);

    auto ms = m->get_state();

    CHECK(ms->name() == m->name());
    CHECK(ms->score() == m->score());
    CHECK(ms->size() == m->residues().size());
    CHECK(ms->end_ids() == m->end_ids());
    CHECK(ms->block_end_add() == m->block_end_add());

    int i = -1;
    for (auto const & end : m->ends()) {
        i++;
        CHECK(end->name() == ms->end_names()[i]);
        CHECK(end->d() == ms->end_states()[i]->d());
        CHECK(math::are_xyzMatrix_equal(end->r(), ms->end_states()[i]->r()));
    }

    SUBCASE("test copy constructor") {
        auto ms_copy = std::make_shared<motif::MotifState>(*ms);

        int i = -1;
        for (auto const & end : ms->end_states()) {
            i++;
            CHECK(ms_copy->end_names()[i] == ms->end_names()[i]);
            CHECK(ms_copy->end_states()[i]->d() == ms->end_states()[i]->d());
            CHECK(math::are_xyzMatrix_equal(ms_copy->end_states()[i]->r(),
                                              ms->end_states()[i]->r()));
        }

    }

    SUBCASE("test stringifying motif state") {
        auto s = ms->to_str();
        auto ms_copy = std::make_shared<motif::MotifState>(s);

        int i = -1;
        for (auto const & end : ms->end_states()) {
            i++;
            CHECK(ms_copy->end_names()[i] == ms->end_names()[i]);
            CHECK(math::are_xyzVector_equal(ms_copy->end_states()[i]->d(),
                                              ms->end_states()[i]->d()));

            CHECK(math::are_xyzMatrix_equal(ms_copy->end_states()[i]->r(),
                                              ms->end_states()[i]->r()));
        }


    }

    SUBCASE("compare motif state align to motif align") {
        auto m2 = std::make_shared<motif::Motif>(*m);
        auto m_aligned = get_aligned_motif(m->ends()[1], m2->ends()[0], m2);

        auto ms2 = std::make_shared<motif::MotifState>(*ms);
        auto ms3 = std::make_shared<motif::MotifState>(*ms);
        get_aligned_motif_state(ms->end_states()[1], ms2, ms3);


        CHECK(math::are_xyzVector_equal(m_aligned->ends()[1]->d(),
                                          ms2->end_states()[1]->d()));

        CHECK(math::are_xyzMatrix_equal(m_aligned->ends()[1]->r(),
                                          ms2->end_states()[1]->r()));

    }

    SUBCASE("test getting motif end information") {
        REQUIRE_NOTHROW(ms->get_end_index("A5-A6"));
        REQUIRE_THROWS_AS(ms->get_end_index("FAKE"), motif::MotifStateException);

        REQUIRE_NOTHROW(ms->get_end_state("A5-A6"));
        REQUIRE_THROWS_AS(ms->get_end_state("FAKE"), motif::MotifStateException);

        auto end = ms->end_states()[0];

        CHECK(ms->get_end_index(end) == 0);

    }

    SUBCASE("Test behavior between basepair and basepair state") {
        auto m_copy = std::make_shared<motif::Motif>(*m);

        auto bp = m_copy->ends()[0];
        auto bp_target = m_copy->ends()[1];

        // make sure its a deep copy
        auto bp_state = bp->state();
        auto bp_state_target = bp_target->state();
        auto s = bp_state->to_str();
        auto dummy = structure::BasepairState();
        bp_state = std::make_shared<structure::BasepairState>(s);

        bp_state_target->get_transforming_r_and_t(*bp_state, dummy);

        bp_state->transform(dummy.r().transposed(), dummy.d());
        bp->transform(dummy.r().transposed(), dummy.d());

        CHECK(math::are_xyzVector_equal(bp_state->d(), bp->d()));

    }

    SUBCASE("test transformation of motif states") {
        auto ms = m->get_state();
        auto ms_copy = std::make_shared<motif::MotifState>(*ms);
        auto m_copy = std::make_shared<motif::Motif>(*m);

        auto t = math::Point(0, 0, 0);
        auto r = rotation_about_x_axis(60);

        for (int i = 0; i < 5; i++) {
            //align_motif(ms->end_states()[0], m_copy->ends()[0], m_copy);
            //m_copy->to_pdb("frame."+std::to_string(i)+".pdb");
            ms->transform(r, t);
        }

        // should be at a 300 degree rotation thus not in the same place
        CHECK(!math::are_xyzVector_equal(ms->end_states()[1]->d(), ms_copy->end_states()[1]->d()));

        ms->transform(r, t);
        auto dist = ms->end_states()[1]->d().distance(ms_copy->end_states()[1]->d());
        CHECK(dist < 0.1);


    }

    SUBCASE("test transformation with aligner ") {
        auto aligner = motif::MotifStateAligner();

        auto ms = m->get_state();
        auto ms1 = std::make_shared<motif::MotifState>(*ms);
        auto ms2 = std::make_shared<motif::MotifState>(*ms);
        auto ms3 = std::make_shared<motif::MotifState>(*ms);
        auto ms4 = std::make_shared<motif::MotifState>(*ms);

        aligner.get_aligned_motif_state(ms1->end_states()[1], ms2, ms3);
        aligner.get_aligned_motif_state(ms1->end_states()[1], ms4);

        CHECK(math::are_xyzVector_equal(ms2->end_states()[1]->d(), ms4->end_states()[1]->d()));
        CHECK(math::are_xyzMatrix_equal(ms2->end_states()[1]->r(), ms4->end_states()[1]->r()));

    }

}




























