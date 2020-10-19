//
//  direction_of_gc.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__direction_of_gc__
#define __RNAMake__direction_of_gc__

#include <stdio.h>

#include "eternabot/strategy.h"

namespace eternabot {
    
class DirectionofGCPairsinMultiLoops : public Strategy {
public:
    DirectionofGCPairsinMultiLoops() {
        mean_ = 85.2869664088;
        stdev_ = 26.9535204308;
        params_ = Floats(2);
        params_[0] = -5.71;
        params_[1] = -6.63;
        name_ = "DirectionofGCPairsinMultiLoops";

    }
    
    ~DirectionofGCPairsinMultiLoops() = default;
    
    float
    score(FeaturesOP const & features) {
        float penalty = 0;
        int i = 0;
        for(auto const & m : features->multi_loops) {
            i = 0;
            for(auto const & end : m->ends()) {
                if(i == 0) {
                    if(secondary_structure::get_bp_type(*end) == secondary_structure::BPType::CG) {
                        continue;
                    }
                    else if(secondary_structure::get_bp_type(*end) == secondary_structure::BPType::GC) {
                        penalty += params_[0];
                    }
                    else {
                        penalty += params_[1];
                    }
                }
                else {
                    if(secondary_structure::get_bp_type(*end) == secondary_structure::BPType::GC) {
                        continue;
                    }
                    else if(secondary_structure::get_bp_type(*end) == secondary_structure::BPType::CG) {
                        penalty += params_[0];
                    }
                    else {
                        penalty += params_[1];
                    }
                }
                i++;
            }
        }
        return 100 + penalty;
    }

private:

};


    
}

#endif /* defined(__RNAMake__direction_of_gc__) */
