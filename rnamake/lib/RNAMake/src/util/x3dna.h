//
//  x3dna.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__x3dna__
#define __RNAMake__x3dna__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//RNAMake Headers
#include "base/string.h"
#include "util/settings.h"

class X3dna {
public:
    
    X3dna();
    
    ~X3dna() {}

public:
    void
    generate_ref_frame(String const &);
    
    void
    generate_dssr_file(String const &);
    
private:
    
    String
    _get_ref_frame_path(String const &);
    
    String
    _get_dssr_file_path(String const &);
    
    
    
private:
    String bin_path_;
    
};


#endif /* defined(__RNAMake__x3dna__) */
