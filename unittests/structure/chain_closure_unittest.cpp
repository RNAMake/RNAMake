//
// Created by Joseph Yesselman on 4/4/18.
//


#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/pdb_parser.h"
#include "structure/close_chain.h"


TEST_CASE( "Test chain closure", "[ChainClosure]" ) {

    SECTION("test most basic fix") {
        auto m_path = base::unittest_resource_dir() + "/motifs/HELIX.IDEAL/HELIX.IDEAL.pdb";
        auto parser = structure::PDBParser();
        auto residues = parser.parse(m_path);
        auto chains = structure::ChainOPs();
        connect_residues_into_chains(residues, chains);

        close_chain(chains[0]);
    }

    SECTION("test fixing missing phosphates") {
         
        auto m_path = base::base_dir() + "//unittests/unittest_resources/motifs/BP.0.22.pdb";
        auto parser = structure::PDBParser();
        auto residues = parser.parse(m_path);
        auto chains = structure::ChainOPs();
        connect_residues_into_chains(residues, chains);
        //chains[1]->to_pdb("org_test.pdb");
        if(chains.size() > 1) { 
            close_chain(chains[1]);
        }
        //chains[1]->to_pdb("test.pdb");
        

    }

}
