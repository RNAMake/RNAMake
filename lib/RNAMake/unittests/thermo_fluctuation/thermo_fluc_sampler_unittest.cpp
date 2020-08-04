
//headers for testing
#include <thermo_fluctuation/graph/sampler.h>
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/graph/sampler.h"

#include <math/quaternion.h>
#include <math/euler.h>
#include <resources/resource_manager.h>

math::Matrix
rotation_about_x_axis(
        float degrees) {
    float a = 3.14/180*degrees;
    return math::Matrix(1.0, 0.0, 0.0,
                        0.0, cos(a), sin(a),
                        0.0, -sin(a), cos(a));
}

math::Matrix
rotation_about_y_axis(
        float degrees) {
    float b = 3.14/180*degrees;
    return math::Matrix(cos(b), 0.0, -sin(b),
                        0.0, 1.0, 0.0,
                        sin(b), 0.0, cos(b));
}

math::Matrix
rotation_about_z_axis(
        float degrees) {
    float g = 3.14/180*degrees;
    return math::Matrix(cos(g), sin(g), 0.0,
                        -sin(g), cos(g), 0.0,
                        0.0, 0.0, 1.0);
}

math::Matrix
rotation_between_frames(
        math::Matrix const & r1,
        math::Matrix const & r2) {

    auto ref_T = math::Matrix();
    auto r = math::Matrix();
    transpose(r1, ref_T);
    dot(ref_T, r2, r);
    r.unitarize();
    return r;
}

float
rmsd_between_basepairs(
        structure::Basepair const & bp1,
        structure::Basepair const & bp2) {
    auto sum = 0.0f;
    for(int i = 0; i < bp1.atoms().size(); i++) {
        sum += (bp1.atoms()[i]->coords().x() - bp2.atoms()[i]->coords().x())*(bp1.atoms()[i]->coords().x() - bp2.atoms()[i]->coords().x());
        sum += (bp1.atoms()[i]->coords().y() - bp2.atoms()[i]->coords().y())*(bp1.atoms()[i]->coords().y() - bp2.atoms()[i]->coords().y());
        sum += (bp1.atoms()[i]->coords().z() - bp2.atoms()[i]->coords().z())*(bp1.atoms()[i]->coords().z() - bp2.atoms()[i]->coords().z());
    }
    sum /= bp1.atoms().size();
    sum = sqrt(sum);

    return sum;


}


TEST_CASE( "Test Thermo Flucuation Sampler ", "[thermo_fluctuation::ThermoFluctuationSampler]" ) {

    auto & rm = resources::Manager::instance();

    SECTION("test moving one frame") {
        auto sampler = thermo_fluctuation::ThermoFlucSampler();
        auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>();
        mset->add_ensemble(resources::Manager::instance().motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(resources::Manager::instance().motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(resources::Manager::instance().motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(resources::Manager::instance().motif_state_ensemble("GG_LL_CC_RR"));
        sampler.setup(mset);

        auto names = Strings();
        for(auto const & n : *sampler.mst()) { names.push_back(n->data()->name()); }
        
        while(sampler.next() == 0) {}
        auto new_names = Strings();

        for(auto const & n : *sampler.mst()) { new_names.push_back(n->data()->name()); }

        int diff = 0;
        for(int i = 0; i < names.size(); i++) {
            if(names[i] != new_names[i]) { diff = 1;}
        }
        
        //check to make sure is actually a different motif state present
        REQUIRE(diff == 1);
        
    }

    SECTION("test graph sampler") {

        SECTION("test random number generator") {
            auto rng = util::RandomNumberGenerator();
            int max = 1;
            int fail = 0;
            for(int i = 0; i < 1000000; i++) {
                if(rng.randrange(max) == max) {
                    fail = 1;
                }
            }
            REQUIRE(fail == 0);
        }

        auto mseg = std::make_shared<motif_data_structure::MotifStateEnsembleGraph>();
        auto & rm = resources::Manager::instance();
        mseg->add_ensemble(*rm.motif_state_ensemble("GG_LL_CC_RR"));
        mseg->add_ensemble(*rm.motif_state_ensemble("GG_LL_CC_RR"));

        auto sampler = thermo_fluctuation::graph::Sampler(*mseg);
        auto msg = sampler.get_initial_state();

        REQUIRE(msg->size() == mseg->size());

        auto names = Strings();
        for(auto const & n : *msg) { names.push_back(n->data()->name()); }

        for(int i = 0; i < 1000; i++) { sampler.next(msg); }
        auto new_names = Strings();
        for(auto const & n : *msg) { new_names.push_back(n->data()->name()); }

        int diff = 0;
        for(int i = 0; i < names.size(); i++) {
            if(names[i] != new_names[i]) { diff = 1;}
        }

        //check to make sure is actually a different motif state present
        REQUIRE(diff == 1);

    }

    SECTION("test rotation invariance") {
        auto m = rm.motif("HELIX.IDEAL.1");
        auto ms = m->get_state();
        auto t = math::Point(0, 0, 0);
        auto dummy = structure::BasepairState();

        ms->end_states()[0]->get_transformed_state(*ms->end_states()[1], dummy);
        auto org_dist = ms->end_states()[0]->d().distance(ms->end_states()[1]->d());
        auto org_dummy = structure::BasepairState(dummy);
        auto org_euler = math::Vector(3);
        //math::calc_euler(org_dummy.r(), )

        for(int i = 0; i < 1000; i++) {
            auto r = math::get_random_quaternion().get_rotation_matrix();

            ms->transform(r.transposed(), t);
            ms->end_states()[0]->get_transformed_state(*ms->end_states()[1], dummy);
            auto dist = ms->end_states()[0]->d().distance(ms->end_states()[1]->d());

            auto r_diff = org_dummy.r().difference(dummy.r());
            if(r_diff > 0.1) { FAIL("orientation has changed"); }
            if(abs(org_dist - dist) > 0.1) { FAIL("distance has changed"); }
        }

        ms->end_states()[0]->get_transformed_state(*ms->end_states()[1], dummy);
        auto new_dist = ms->end_states()[0]->d().distance(ms->end_states()[1]->d());

        REQUIRE(abs(org_dist - new_dist) < 0.1);

        auto new_m = rm.get_motif_from_state(ms);
        auto new_ms = new_m->get_state();
        new_ms->end_states()[0]->get_transformed_state(*new_ms->end_states()[1], dummy);

    }

    SECTION("test axis angle scoring") {
        auto m = rm.motif("HELIX.IDEAL.1");
        auto m1 = rm.motif("HELIX.IDEAL.1");

        auto bp = m->ends()[0];
        auto bp_org = m1->ends()[0];
        //auto org_bp = std::make_shared<structure::Basepair>(*bp);
        auto rx = rotation_about_x_axis(10);
        auto ry = rotation_about_y_axis(10);
        auto rz = rotation_about_z_axis(10);
        auto t = math::Point();

        auto aa = math::AxisAngle();

        std::ofstream out;
        out.open("rot_diff_vs_rmsd.csv");
        out << "rot_diff,rmsd" << std::endl;

        for (int i = 0; i < 36; i++) {
            bp->transform(rx, t);
            for(int j = 0; j < 36; j++) {
                bp->transform(ry, t);
                for(int k = 0; k < 36; k++) {
                    bp->transform(rz, t);
                    auto diff = rotation_between_frames(bp_org->r(), bp->r());
                    //math::axis_angle_from_matrix(diff, aa);
                    //std::cout << math::degrees(aa.angle) << " " << rmsd_between_basepairs(*bp, *bp_org) << std::endl;
                    out << bp_org->r().difference(bp->r()) << "," << rmsd_between_basepairs(*bp, *bp_org) << std::endl;
                }
            }
        }
        out.close();
    }
    
}






























