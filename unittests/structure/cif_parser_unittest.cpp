#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/cif_parser.h"
#include "structure/pdb_parser.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.h"


// what do I want to test?
// 1. doesn't return anything for a .cif with nothing in it
// 2. returns the same for the pdb parser
// 3. generic variation (i.e. small, medium, multi-chain, etc
// 4. raises for the errors

TEST_CASE( "Test CIF Parser", "[CIFParser]" ) {
    
    const auto cif_path = base::unittest_resource_dir() + "/cifs/";
    const auto pdb_path = base::unittest_resource_dir() + "/pdbs/";
    auto m_path = String{"../../data/255D.cif"};
    auto parser = structure::CIFParser();

    SECTION("Simple comparison with PDBParser") {
        auto cifparser = structure::CIFParser();   
        auto pdbparser = structure::PDBParser();
        
        auto pdb_residues = pdbparser.parse(pdb_path +"/255D.pdb");
        auto cif_residues = cifparser.parse(cif_path + "/255D.cif");
        
        REQUIRE(pdb_residues.size() == cif_residues.size());
        
        auto pdb_chains = structure::ChainOPs{};
        auto cif_chains = structure::ChainOPs{};
        
        connect_residues_into_chains(pdb_residues,pdb_chains);
        connect_residues_into_chains(cif_residues,cif_chains);
        
        REQUIRE(pdb_chains.size() == cif_chains.size());
        
        REQUIRE(are_chains_equal(pdb_chains[0],cif_chains[0],0));
    }

    SECTION("Medium difficulty comparison with PDBParser") {
        auto cifparser = structure::CIFParser();   
        auto pdbparser = structure::PDBParser();
        
        auto pdb_residues = pdbparser.parse(pdb_path + "/4OJI.pdb");
        auto cif_residues = cifparser.parse(cif_path+"/4OJI.cif");
        
        REQUIRE(pdb_residues.size() == cif_residues.size());
        
        auto pdb_chains = structure::ChainOPs{};
        auto cif_chains = structure::ChainOPs{};
        
        connect_residues_into_chains(pdb_residues,pdb_chains);
        connect_residues_into_chains(cif_residues,cif_chains);
        
        REQUIRE(pdb_chains.size() == cif_chains.size());
        
        auto pdb_it = pdb_chains.cbegin();
        auto cif_it = cif_chains.cbegin();
        
        const auto cif_last = cif_chains.cend();

        for( ; cif_it != cif_last; ++cif_it,++pdb_it) {
            REQUIRE(are_chains_equal(*pdb_it,*cif_it,0)); 
        }
 
    }

    SECTION("Mixed RNA/non-rna") {
        auto cifparser = structure::CIFParser();   
        auto pdbparser = structure::PDBParser();
        
        auto pdb_residues = pdbparser.parse(pdb_path+"/4K4W.pdb");
        auto cif_residues = cifparser.parse(cif_path+"/4K4W.cif");
        
        REQUIRE(pdb_residues.size() == cif_residues.size());
        
        auto pdb_chains = structure::ChainOPs{};
        auto cif_chains = structure::ChainOPs{};
        
        connect_residues_into_chains(pdb_residues,pdb_chains);
        connect_residues_into_chains(cif_residues,cif_chains);
        
        REQUIRE(pdb_chains.size() == cif_chains.size());
        
        auto pdb_it = pdb_chains.cbegin();
        auto cif_it = cif_chains.cbegin();
        
        const auto cif_last = cif_chains.cend();

        for( ; cif_it != cif_last; ++cif_it,++pdb_it) {
            REQUIRE(are_chains_equal(*pdb_it,*cif_it,0)); 
        }

    }

    SECTION("Testing for Parsing errors") {
        auto cifparser = structure::CIFParser{};
        // incomplete file
        REQUIRE_THROWS(cifparser.parse(cif_path + "/incomplete.cif"));
        // missing file
        REQUIRE_THROWS(cifparser.parse(cif_path + "/nontexistent.cif"));
        // bad headers 
        REQUIRE_THROWS(cifparser.parse(cif_path + "/bad_headers.cif"));
        }

}

