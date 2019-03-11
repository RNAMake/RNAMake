
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/pdb_parser.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.hpp"


TEST_CASE( "Test PDB Parser", "[PDBParser]" ) {
    auto m_path = base::base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
    auto parser = PDBParser();
    auto residues = parser.parse(m_path);
    auto chains = ChainOPs();
    connect_residues_into_chains(residues, chains);
    auto s = std::make_shared<Structure>(chains);

    SECTION("compare parsed structure to stringifed structure") {
    
        auto rts = ResidueTypeSet();
        String path = base::unittest_resource_dir() + "/structure/test_str_to_structure.dat";
        auto lines =base::get_lines_from_file(path);
        auto s_org = std::make_shared<Structure>(lines[0], rts);
        
        REQUIRE(are_structures_equal(s, s_org, 0));
    }
}
