//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/residue.h"

namespace secondary_structure {

ResType
convert_res_name_to_type(
        char c) {
    if (c == 'A')      { return ResType::ADE;  }
    else if (c == 'C') { return ResType::CYT;  }
    else if (c == 'G') { return ResType::GUA;  }
    else if (c == 'U') { return ResType::URA;  }
    else if (c == 'T') { return ResType::URA;  }
    else if (c == 'N') { return ResType::NONE; }
    else {
        throw secondary_structure::Exception("incorrect character for secondary string");
    }
}

}
