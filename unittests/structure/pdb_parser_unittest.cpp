
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/pdb_parser.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.h"


TEST_CASE( "Test PDB Parser", "[PDBParser]" ) {
    auto m_path = base::base_dir() + "/unittests/unittest_resources/motifs/p4p6/p4p6.pdb";
    auto parser = structure::PDBParser();
    auto residues = parser.parse(m_path);
    auto chains = structure::ChainOPs();
    connect_residues_into_chains(residues, chains);
    auto s = std::make_shared<structure::Structure>(chains);

    SECTION("compare parsed structure to stringifed structure") {
    
        auto rts = structure::ResidueTypeSet();
        String path = base::unittest_resource_dir() + "/structure/test_str_to_structure.dat";
        auto lines =base::get_lines_from_file(path);
        auto s_org = std::make_shared<structure::Structure>(lines[0], rts);
         
        REQUIRE(are_structures_equal(s, s_org, 0));
    }
}
