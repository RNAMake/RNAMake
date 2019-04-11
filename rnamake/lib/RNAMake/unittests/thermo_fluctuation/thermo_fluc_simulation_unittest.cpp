
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "thermo_fluctuation/thermo_fluc_simulation.h"
#include <thermo_fluctuation/graph/simulation.h>

// copy old code scoring
class OldFrameScorer : public thermo_fluctuation::graph::Scorer {
public:
    OldFrameScorer() : thermo_fluctuation::graph::Scorer() {}

    ~OldFrameScorer() {}

    thermo_fluctuation::graph::Scorer *
    clone() const { return new OldFrameScorer(*this); };

public:

    inline
    float
    score(
            structure::BasepairState const & state_1,
            structure::BasepairState const & state_2) {
        score_ = state_1.d().distance(state_2.d());
        r_diff_ = state_1.r().difference(state_2.r());
        flipped_ = state_2.r().get_flip_orientation();
        r_diff_flip_ = state_1.r().difference(flipped_);
        if (r_diff_ > r_diff_flip_) { score_ += r_diff_flip_; }
        else { score_ += r_diff_; }
        return score_;
    }

private:
    float r_diff_, r_diff_flip_,  score_;
    math::Matrix flipped_;
};

TEST_CASE( "Test Thermo Flucuation Simulation ", "[thermo_fluctuation::ThermoFluctuationSimulation" ) {

    auto & rm = resources::Manager::instance();

    SECTION("test running simulation") {
        auto tfs = thermo_fluctuation::ThermoFlucSimulation();
        tfs.set_option_value("cutoff", 20);
        auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>();
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        tfs.setup(mset, 0, 3, 0, 1);
        auto count = tfs.run();
        REQUIRE(count > 0);
    }
    
    SECTION("catch common error in setting up simulation") {
        auto tfs = thermo_fluctuation::ThermoFlucSimulation();
        // cant execute run until setup has been called
        REQUIRE_THROWS_AS(tfs.run(), thermo_fluctuation::ThermoFlucSimulationException);
        
        auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>();
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));

        // try using node indexes and end indexes that dont exist
        REQUIRE_THROWS_AS(tfs.setup(mset, 5, 2, 0, 1), thermo_fluctuation::ThermoFlucSimulationException);
        REQUIRE_THROWS_AS(tfs.setup(mset, 0, 5, 0, 1), thermo_fluctuation::ThermoFlucSimulationException);
        REQUIRE_THROWS_AS(tfs.setup(mset, 0, 2, 5, 1), thermo_fluctuation::ThermoFlucSimulationException);
        REQUIRE_THROWS_AS(tfs.setup(mset, 0, 2, 0, 5), thermo_fluctuation::ThermoFlucSimulationException);

    }

    SECTION("Test graph simulation") {

        SECTION("Test scorer") {
            auto scorer = thermo_fluctuation::graph::FrameScorer();
            auto ms = rm.motif_state("HELIX.IDEAL.2");

            auto score = scorer.score(*ms->end_states()[0], *ms->end_states()[0]);
            REQUIRE(score == 0.0);
        }

        SECTION("Test simulation") {
            auto scorer = std::make_shared<thermo_fluctuation::graph::FrameScorer>();
            auto sim = thermo_fluctuation::graph::Simulation(scorer);

            auto mseg = std::make_shared<motif_data_structure::MotifStateEnsembleGraph>();
            mseg->add_ensemble(*rm.motif_state_ensemble("GG_LL_CC_RR"));
            mseg->add_ensemble(*rm.motif_state_ensemble("GG_LL_CC_RR"), data_structure::NodeIndexandEdge{0, 1});

            sim.setup(*mseg, data_structure::NodeIndexandEdge{1, 1}, data_structure::NodeIndexandEdge{0, 0});
            sim.set_option_value("cutoff", 100);
            REQUIRE(sim.next() == true);
        }

        /*SECTION("compare to old code") {
            auto scorer = std::make_shared<OldFrameScorer>();
            auto sim = thermo_fluctuation::graph::Simulation(scorer);

            auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>();
            auto mseg = std::make_shared<motif_data_structure::MotifStateEnsembleGraph>();

            for(int i = 0; i < 10; i++) {
                mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
                if(i == 0) {
                    mseg->add_ensemble(*rm.motif_state_ensemble("GG_LL_CC_RR"));
                }
                else {
                    mseg->add_ensemble(*rm.motif_state_ensemble("GG_LL_CC_RR"), data_structure::NodeIndexandEdge{i-1, 1});
                }
            }

            auto tfs = thermo_fluctuation::ThermoFlucSimulation();
            tfs.set_option_value("cutoff", 25);
            tfs.setup(mset, 0, 9, 0, 1);
            auto avg_1 = 0;

            for(int i = 0; i < 500; i++) {
                auto count = tfs.run();
                avg_1 += count;
            }
            avg_1 /= 500;

            sim.set_option_value("cutoff", 25);
            sim.setup(*mseg, data_structure::NodeIndexandEdge{0, 0}, data_structure::NodeIndexandEdge{9, 1});
            int avg_2 = 0, new_count = 0;
            for(int j = 0; j < 500; j++) {
                new_count = 0;
                for (int i = 0; i < 100000; i++) {
                    new_count += sim.next();
                }
                avg_2 += new_count;
            }
            avg_2 /= 500;

            std::cout << avg_1 << " " << avg_2 << std::endl;
        }*/
    }
    
}