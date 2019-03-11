
#include "../../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/pdb_parser.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.hpp"


TEST_CASE( "Load all PDBs", "[PDBParser]" ) {
    auto parser = PDBParser();
    
    SECTION("compare all parsed structures to python counterpart") {
        
        auto path = unittest_resource_dir() + "/structure/seqs_in_structures.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto base = String("/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas/");
        
        auto fail = 0;
        for(auto const & l : lines) {
            auto spl = base::split_str_by_delimiter(l, " ");
            if(spl.size() < 2) { continue; }
            auto s_path = base + "/" + spl[0] + "/" + spl[0] + ".pdb";
            auto res = parser.parse(s_path);
            
            auto chains = ChainOPs();
            connect_residues_into_chains(res, chains);
            auto seq = String("");
            for(auto const & c : chains) {
                for(auto const & r : c->residues()) {
                    seq += r->name();
                }
            }
            
            if(spl[1].length() != seq.length()) {
                std::cout << spl[0] << " " << spl[1] << " " << seq << std::endl;
                fail = 1;
            }
        }
        
        REQUIRE(fail == 0);
        
    }
}


