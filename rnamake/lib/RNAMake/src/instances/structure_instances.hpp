//
//  structure_instances.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/13/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef structure_instances_hpp
#define structure_instances_hpp

#include <stdio.h>

//RNAMake Headers
#include "base/settings.h"
#include "base/file_io.h"
#include "structure/residue.h"

namespace instances {

    

inline
structure::ResidueOP
residue() {
    auto unittest_path = base::base_dir() + "/rnamake/lib/RNAMake/unittests/resources/";
    auto path = unittest_path + "residue/test_str_to_residue.dat";
    Strings lines =base::get_lines_from_file(path);
    structure::ResidueTypeSet rts;
    auto r = std::make_shared<structure::Residue>(lines[0], rts);
    return r;
}
    
    
}

#endif /* structure_instances_hpp */
