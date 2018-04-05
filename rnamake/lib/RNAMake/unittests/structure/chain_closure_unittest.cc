//
// Created by Joseph Yesselman on 4/4/18.
//


#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/pdb_parser.h"
#include "structure/close_chain.h"


TEST_CASE( "Test chain closure", "[ChainClosure]" ) {
    auto m_path = base_dir() + "/rnamake/unittests/resources/motifs/HELIX.IDEAL/HELIX.IDEAL.pdb";
    auto parser = PDBParser();
    auto residues = parser.parse(m_path);
    auto chains = ChainOPs();
    connect_residues_into_chains(residues, chains);

    close_chain(chains[0]);

}