
//headers for testing
#include <thermo_fluctuation/graph/sampler.h>
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/graph/sampler.h"

TEST_CASE( "Test Thermo Flucuation Sampler ", "[thermo_fluctuation::ThermoFluctuationSampler]" ) {
    
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

    
}