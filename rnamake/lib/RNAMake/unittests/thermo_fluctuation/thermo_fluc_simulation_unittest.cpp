
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "thermo_fluctuation/thermo_fluc_simulation.h"

TEST_CASE( "Test Thermo Flucuation Simulation ", "[ThermoFluctuationSimulation" ) {

    auto & rm = resources::Manager::instance();

    SECTION("test running simulation") {
        auto tfs = ThermoFlucSimulation();
        tfs.set_option_value("cutoff", 20);
        auto mset = std::make_shared<MotifStateEnsembleTree>();
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        tfs.setup(mset, 0, 3, 0, 1);
        auto count = tfs.run();
        REQUIRE(count > 0);
      
    }
    
    SECTION("catch common error in setting up simulation") {
        auto tfs = ThermoFlucSimulation();
        // cant execute run until setup has been called
        REQUIRE_THROWS_AS(tfs.run(), ThermoFlucSimulationException);
        
        auto mset = std::make_shared<MotifStateEnsembleTree>();
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));
        mset->add_ensemble(rm.motif_state_ensemble("GG_LL_CC_RR"));

        // try using node indexes and end indexes that dont exist
        REQUIRE_THROWS_AS(tfs.setup(mset, 5, 2, 0, 1), ThermoFlucSimulationException);
        REQUIRE_THROWS_AS(tfs.setup(mset, 0, 5, 0, 1), ThermoFlucSimulationException);
        REQUIRE_THROWS_AS(tfs.setup(mset, 0, 2, 5, 1), ThermoFlucSimulationException);
        REQUIRE_THROWS_AS(tfs.setup(mset, 0, 2, 0, 5), ThermoFlucSimulationException);

    }
    
    
}