//
//  motif_type.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_type__
#define __RNAMake__motif_type__

#include <stdio.h>
#include "types.h"

enum MotifType {
    TWOWAY            = 0,
    NWAY              = 1,
    HAIRPIN           = 2,
    TCONTACT_HP_HP    = 3,
    TCONTACT_H_HP     = 4,
    TCONTACT_H_H      = 5,
    T_T               = 6,
    T_T_T             = 7,
    TWOWAY_SEGMENTS   = 8,
    HELIX             = 9,
    SSTRAND           = 10,
    TCONTACT          = 11,
    UNKNOWN           = 99,
    ALL               = 999
};

String const
type_to_str(MotifType const);

MotifType const
str_to_type(String const);


#endif /* defined(__RNAMake__motif_type__) */