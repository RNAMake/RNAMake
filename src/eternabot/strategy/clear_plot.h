//
//  clear_plot.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__clear_plot__
#define __RNAMake__clear_plot__

#include <stdio.h>

#include <eternabot/strategy.h>

namespace eternabot {

class CleanPlotStackCapsandSafeGC : public Strategy {
public:
    CleanPlotStackCapsandSafeGC() {
        params_ = std::vector<float>(4);
        params_[0] = 0.100135909783;
        params_[1] = 1.76372839803;
        params_[2] = 3.11085515568;
        params_[3] = 0.966424875922;
        mean_ = 82.3365692703;
        stdev_ = 12.050647236;
    }
    
    ~CleanPlotStackCapsandSafeGC() {}
    
public:
    
    float
    score(FeaturesOP const & features) {
        
        float penalty = 0.0;
        int n = features->length;
        float npairs = features->gc + features->gu + features->ua;
        
        int i_index, j_index;
        int fail = 0;
        for(int i = 0; i < n*n; i++) {
            if(features->dotplot[i].p < 0.0001) { continue; }
            fail = 0;
            i_index = features->dotplot[i].i;
            j_index = features->dotplot[i].j;
            
            if     (features->pairmap.find(i_index) == features->pairmap.end()) { fail = 1; }
            else if(features->pairmap.at(i_index) != j_index)                   { fail = 1; }

            if(fail) { penalty += features->dotplot[i].p; }
            
        }
        
        float plotscore = 0;
        float gc_penalty = 0;
        if(npairs > 0) {
            plotscore = (1.0 - (penalty / npairs));
            if(features->gc / npairs > params_[3]) {
                gc_penalty = 1;
            }
        }
        
        float cap_score = 0.0f, stack_count = 0.0f;
        int length;
        for(auto const & helix : features->helices) {
            stack_count++;
            if     (helix->basepairs().size() == 1) {
                if(secondary_structure::is_gc_pair(helix->basepairs()[0])) { cap_score += 1;   }
            }
            else if(helix->basepairs().size() == 2) {
                if(secondary_structure::is_gc_pair(helix->basepairs()[0])) { cap_score += 0.5; }
                if(secondary_structure::is_gc_pair(helix->basepairs()[1])) { cap_score += 0.5; }
                
            }
            
            else if(helix->basepairs().size() == 3) {
                if(secondary_structure::is_gc_pair(helix->basepairs()[0])) { cap_score += 0.4; }
                if(secondary_structure::is_gc_pair(helix->basepairs()[1])) { cap_score += 0.4; }
                if(secondary_structure::is_gc_pair(helix->basepairs()[2])) { cap_score += 0.4; }
            }
            
            else {
                length = (int)helix->basepairs().size();
                if(secondary_structure::is_gc_pair(helix->basepairs()[0])) { cap_score += 1.0 / 3.0; }
                if(secondary_structure::is_gc_pair(helix->basepairs()[1])) { cap_score += 1.0 / 6.0; }
                if(secondary_structure::is_gc_pair(helix->basepairs()[length-2])) { cap_score += 1.0 / 6.0; }
                if(secondary_structure::is_gc_pair(helix->basepairs()[length-1])) { cap_score += 1.0 / 3.0; }
            }
        }
        
        if(stack_count > 0) {
            cap_score = cap_score / stack_count;
        }
        

        
        float score =  (2.0 + cap_score * params_[1] + plotscore * params_[0] - gc_penalty * params_[2]) * 25;
        return score;
    }
    
};

}

#endif /* defined(__RNAMake__clear_plot__) */
