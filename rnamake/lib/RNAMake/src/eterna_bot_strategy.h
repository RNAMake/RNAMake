//
//  eterna_bot_strategy.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/22/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__eterna_bot_strategy__
#define __RNAMake__eterna_bot_strategy__

#include <stdio.h>
#include <vector>
#include <memory>
#include <math.h>
#include <iostream>
#include <map>
#include "vienna.h"
#include "secondary_structure_tree.h"

struct EternabotScorerData {
public:
    
    EternabotScorerData():
    gu(0), gc(0), ua(0),
    meltpoint(97), fe(0),
    dotplot(NULL),
    length ( 0 ),
    pairmap ( std::map<int, int>() ),
    seq ( String() )
    {}
    
    
public:
    float gu, gc, ua;
    float meltpoint, fe;
    plists dotplot;
    int length;
    std::map<int, int> pairmap;
    SecondaryStructureTree sstree;
    String seq;
    
};


class EternabotStrategy {
public:
    
    virtual
    float
    score(EternabotScorerData const &) = 0;
    
public:
    inline
    float
    mean() const { return mean_; }
    
    inline
    float
    stdev() const { return stdev_; }
    
protected:
    float mean_, stdev_;
    
};

class ABasicTest : public EternabotStrategy {
public:
    ABasicTest() {
        params_ = std::vector<float>(7);
        params_[0] = 0.303497971269;
        params_[1] = 92.9893755247;
        params_[2] = -1.37878787864;
        params_[3] = 0.512804062262;
        params_[4] = 0.477932936507;
        params_[5] = 84.4793979751;
        params_[6] = 124.345433009;
        mean_ = 83.5007560083;
        stdev_ = 10.5290224709;

    }
    
    ~ABasicTest() {}
    
    float
    score(EternabotScorerData const & data) {
        float total_pairs = data.gc + data.gu + data.ua;
        float score = 100;
        if(total_pairs > 0) {
            score -= fabs(data.ua / total_pairs - params_[0]) * params_[1];
        }
        float target_fe = params_[2] * total_pairs;
        score -= fabs(target_fe - data.fe) * params_[3];

        if(data.meltpoint < params_[5]) {
            score -= fabs(data.meltpoint - params_[5]) * params_[4];
        }
        else if(data.meltpoint > params_[6]) {
            score -= fabs(data.meltpoint - params_[6]) * params_[4];
        }
        return score;
    }

    
private:
    std::vector<float> params_;
};

class ClearPlotStackCapsandSafeGC : public EternabotStrategy {
public:
    ClearPlotStackCapsandSafeGC() {
        params_ = std::vector<float>(4);
        params_[0] = 0.100135909783;
        params_[1] = 1.76372839803;
        params_[2] = 3.11085515568;
        params_[3] = 0.966424875922;
        mean_ = 82.3365692703;
        stdev_ = 12.050647236;
    }
    
    float
    score(EternabotScorerData const & data) {
        
        float penalty = 0.0;
        int n = data.length;
        float npairs = data.gc + data.gu + data.ua;
        
        int i_index, j_index;
        int fail = 0;
        for(int i = 0; i < n*n; i++) {
            if(data.dotplot[i].p < 0.0001) { continue; }
            fail = 0;
            i_index = data.dotplot[i].i;
            j_index = data.dotplot[i].j;
            
            if     (data.pairmap.find(i_index) == data.pairmap.end()) {
                fail = 1;
            }
            else if(data.pairmap.at(i_index) != j_index) {
                fail = 1;
            }
            
            if(fail) {
                penalty += data.dotplot[i].p;
            }
            
        }

        float plotscore = 0;
        float gc_penalty = 0;
        if(npairs > 0) {
            plotscore = (1.0 - (penalty / npairs));
            if(data.gc / npairs > params_[3]) {
                gc_penalty = 1;
            }
        }
        
        float cap_score = 0.0f, stack_count = 0.0f;
        int length;
        for(auto const & helix : data.sstree.helices()) {
            stack_count++;
            if     (helix.size() == 1) {
                if(helix[0]->bp_type_num() == 0 || helix[0]->bp_type_num() == 1) {
                    cap_score += 1;
                }
            }
            else if(helix.size() == 2) {
                if(helix[0]->bp_type_num() == 0 || helix[0]->bp_type_num() == 1) {
                    cap_score += 0.5;
                }
                if(helix[1]->bp_type_num() == 0 || helix[1]->bp_type_num() == 1) {
                    cap_score += 0.5;
                }

            }
            
            else if(helix.size() == 3) {
                if(helix[0]->bp_type_num() == 0 || helix[0]->bp_type_num() == 1) {
                    cap_score += 0.4;
                }
                if(helix[1]->bp_type_num() == 0 || helix[1]->bp_type_num() == 1) {
                    cap_score += 0.4;
                }
                if(helix[2]->bp_type_num() == 0 || helix[2]->bp_type_num() == 1) {
                    cap_score += 0.4;
                }
            }
            
            else {
                length = (int)helix.size();
                if(helix[0]->bp_type_num() == 0 || helix[0]->bp_type_num() == 1) {
                    cap_score += 1.0 / 3.0;
                }
                if(helix[1]->bp_type_num() == 0 || helix[1]->bp_type_num() == 1) {
                    cap_score += 1.0 / 6.0;
                }
                if(helix[length-2]->bp_type_num() == 0 || helix[length-2]->bp_type_num() == 1) {
                    cap_score += 1.0 / 6.0;
                }
                if(helix[length-1]->bp_type_num() == 0 || helix[length-1]->bp_type_num() == 1) {
                    cap_score += 1.0 / 3.0;
                }
            }
        }
        
        if(stack_count > 0) {
            cap_score = cap_score / stack_count;
        }
        
        float score =  (2.0 + cap_score * params_[1] + plotscore * params_[0] - gc_penalty * params_[2]) * 25;
        return score;
    }
    
private:
    std::vector<float> params_;

};

class BerexTest: public EternabotStrategy {
public:
    BerexTest() {
        params_ = std::vector<float>(10);
        params_[0] = 0.233331320365;
        params_[1] = 67.1928215814;
        params_[2] = 0.0846125555161;
        params_[3] = 97.8926617326;
        params_[4] = 0.124731583984;
        params_[5] = 76.3152913746;
        params_[6] = -64.8000139266;
        params_[7] = -25.1005947191;
        params_[8] = 1.2097925221;
        params_[9] = 33.3729652592;
        params_[10] = 171.705337159;
        params_[11] = 1.28554296532;
        mean_ = 84.0125821249;
        stdev_ = 8.91633847502;
    }
    
    ~BerexTest() {}
    
    float
    score(EternabotScorerData const & data) {
        g_count_ = 0; c_count_ = 0; a_count_ = 0; u_count_ = 0;
        for (auto const & s : data.seq) {
            if      (s == 'G') { g_count_++; }
            else if (s == 'C') { c_count_++; }
            else if (s == 'A') { a_count_++; }
            else               { u_count_++; }
        }
        
        float score = 100;
        score -= fabsf(g_count_ / data.length - params_[0]) * params_[1];
        score -= fabsf(u_count_ / data.length - params_[2]) * params_[3];
        score -= fabsf(c_count_ / data.length - params_[4]) * params_[5];

        if     (data.fe < params_[6]) {
            score -= fabsf(data.fe - params_[6]) * params_[8];
        }
        else if(data.fe > params_[7]) {
            score -= fabsf(data.fe - params_[7]) * params_[8];
        }
        
        if     (data.meltpoint < params_[9]) {
            score -= fabsf(data.meltpoint - params_[9]) * params_[11];
        }
        else if(data.meltpoint > params_[10]) {
            score -= fabsf(data.meltpoint - params_[10]) * params_[11];
        }
        
        return score;
    }
    
private:
    std::vector<float> params_;
    float g_count_, c_count_, a_count_, u_count_;
};

class NumofYellowNucleotidesperLengthofString : public EternabotStrategy {
public:
    NumofYellowNucleotidesperLengthofString() {
        params_ = std::vector<float>(1);
        params_[0] = 791.641998291;
        upper_length_ = Ints(10);
        upper_length_[0] = 0; upper_length_[1] = 1; upper_length_[2] = 2;
        upper_length_[3] = 2; upper_length_[4] = 2; upper_length_[5] = 3;
        upper_length_[6] = 3; upper_length_[7] = 4; upper_length_[8] = 5;
        upper_length_[9] = 4;
        lower_length_ = Ints(10);
        lower_length_[0] = 0; lower_length_[1] = 0; lower_length_[2] = 0;
        lower_length_[3] = 1; lower_length_[4] = 1; lower_length_[5] = 1;
        lower_length_[6] = 2; lower_length_[7] = 3; lower_length_[8] = 2;
        lower_length_[9] = 1;
        mean_ = 91.2420911348;
        stdev_ = 12.5663926344;
    }
    
    ~NumofYellowNucleotidesperLengthofString() {}
    
    float
    score(EternabotScorerData const & data) {
        
        float penalty = 0, count = 0;
        int stack_length;
        float yellow_count = 0;
        for(auto const & helix : data.sstree.helices()) {
            if(helix.size() < 2 || helix.size() > 10) { continue; }
            stack_length = (int)helix.size();
            count ++;
            for(auto const & n : helix) {
                //is a bp of AU or UA
                if(n->bp_type_num() == 2 || n->bp_type_num() == 3) { yellow_count ++; }
            }
            if     (upper_length_[stack_length] < yellow_count) {
                penalty += yellow_count - upper_length_[stack_length];
            }
            else if(lower_length_[stack_length] > yellow_count) {
                penalty += lower_length_[stack_length] - yellow_count;
            }

        }
        if(count == 0) { return 0; }
        
        return 100 - params_[0] * penalty/float(data.length);
    
    }
    
    
    
private:
    Ints upper_length_, lower_length_;
    std::vector<float> params_;
};

//implement at some point
class DirectionofGCPairsinMultiLoops : public EternabotStrategy {
public:
    DirectionofGCPairsinMultiLoops() {
        mean_ = 85.2869664088;
        stdev_ = 26.9535204308;
    }
    
    ~DirectionofGCPairsinMultiLoops() {}
    
    float
    score(EternabotScorerData const & data) {
        return 0;
    }
    
};


typedef std::shared_ptr<EternabotStrategy> EternabotStrategyOP;
typedef std::vector<EternabotStrategyOP>   EternabotStrategyOPs;


#endif /* defined(__RNAMake__eterna_bot_strategy__) */
