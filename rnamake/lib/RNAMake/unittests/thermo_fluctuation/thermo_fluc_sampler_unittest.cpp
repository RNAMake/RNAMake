
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/thermo_fluc_graph_sampler.h"

TEST_CASE( "Test Thermo Flucuation Sampler ", "[ThermoFluctuationSampler]" ) {
    
    SECTION("test moving one frame") {
        auto sampler = ThermoFlucSampler();
        auto mset = std::make_shared<MotifStateEnsembleTree>();
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
        auto sampler = ThermoFlucGraphSampler();
        /*auto mset = std::make_shared<MotifStateEnsembleTree>();
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
        REQUIRE(diff == 1);*/

    }

    
}