//
//  motif_state_search_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search_scorer__
#define __RNAMake__motif_state_search_scorer__

#include <stdio.h>

//RNAMAke Headers
#include "structure/basepair_state.h"
#include "motif_state_search/motif_state_search_node.h"

class MotifStateSearchScorer {
public:
    MotifStateSearchScorer() {}
    
    ~MotifStateSearchScorer() {}
    
public:
    
    void
    set_target(BasepairStateOP const &);
    
    virtual
    inline
    float
    score(
        MotifStateSearchNodeOP const & node ) {
        best_score_ = 1000;
        int i = -1;
        for(auto const & state : node->cur_state()->end_states() ) {
            i++;
            if(i == 0) { continue; }
            
            score_ = new_score_function_new(state, target_, target_flip_);
            
            if(score_ < best_score_) {
                best_score_ = score_;
            }
            
        }
        return best_score_;
    
    }
    
    virtual
    float
    accept_score(MotifStateSearchNodeOP const &);
    
public:
    virtual
    inline
    void
    level_weight(float const & nlevel_weight) {}
    
    virtual
    inline
    void
    ss_score_weight(float const & nss_score_weight) {}
    
protected:
    BasepairStateOP target_, target_flip_;
    float best_score_, score_, r_diff_, r_diff_flip_;
    
};

class MSS_GreedyBestFirstSearch : public MotifStateSearchScorer {
public:
    MSS_GreedyBestFirstSearch() {}
    
public:
    
    inline
    float
    score(MotifStateSearchNodeOP const & node) {
        return MotifStateSearchScorer::score(node);
    }
    
public:
    
};

class MSS_Astar: public MotifStateSearchScorer {
public:
    MSS_Astar() {
        g_ = 0; h_ = 0;
        ss_score_weight_ = 0.25;
        level_weight_ = 2.0;
    }
    
    inline
    float
    score(MotifStateSearchNodeOP const & node) {
        g_ = node->ss_score()*ss_score_weight_;
        if(node->level() > 2) {
            g_ += node->level()*level_weight_;
        }
        h_ = MotifStateSearchScorer::score(node);
        return g_ + h_;
    }
    
public:
    
    inline
    void
    level_weight(float const & nlevel_weight) { level_weight_ = nlevel_weight; }
    
    inline
    void
    ss_score_weight(float const & nss_score_weight) { ss_score_weight_ = nss_score_weight; }
    
private:
    float g_, h_;
    float ss_score_weight_, level_weight_;
};

class MSS_PathFollow: public MotifStateSearchScorer {
public:
    MSS_PathFollow(Points const & path) {
        path_ = path;
        beads_ = Points(10000);
        bead_size_ = 0;
        weight_ = 2.0;
        
    }
    
    inline
    float
    score(MotifStateSearchNodeOP const & node) {
        current_ = node;
        bead_size_ = 0;
        weight_ = 2.0;
        length_ = 0;
        while(1) {
            for(auto const & b : current_->cur_state()->beads()) {
                beads_[bead_size_] = b;
                bead_size_++;
            }
            length_ += current_->cur_state()->end_states()[1]->d().distance(current_->cur_state()->end_states()[0]->d());
            
            current_ = current_->parent();
            if(current_ == nullptr) { break; }
        }
        
        score_ = 0;
        int sum = 0;
        int i = 0;
        for(auto const & b1 : path_) {
            best_ = 1000000;
            for(int i = 0; i < bead_size_; i++) {
                dist_ = b1.distance(beads_[i]);
                if(best_ > dist_) {
                    best_ = dist_;
                }
                if(best_ < 10) { break; }
            }
            if(best_ > 10) {
                score_ += 1;
            }
            else {
                if(sum < i) {
                    score_ += (float)(i - sum)/2;
                }
                sum += 1;
            }
            i++;
        }
        
        //return length_*0.01 + node->level()*0.50 + score_;
        return length_*0.02 + score_;
  
    }
    
    inline
    float
    _score(MotifStateSearchNodeOP const & node) {
        current_ = node;
        bead_size_ = 0;
        weight_ = 2.0;
        while(1) {
            //for(auto const & b : current_->cur_state()->beads()) {
            beads_[bead_size_] = current_->center();
            bead_size_++;
            //}
            current_ = current_->parent();
            if(current_ == nullptr) { break; }
        }
        
        score_ = 0;
        int i = 0;

        for(auto const & b1 : path_) {
            best_ = 1000000;
            for(int i = 0; i < bead_size_; i++) {
                dist_ = b1.distance(beads_[i]);
                if(best_ > dist_) {
                    best_ = dist_;
                }
                //if(best_ < 5) { break; }
            }
            score_ += best_*weight_;
            weight_ *= 0.99;
        }
        
        return score_;
        
    }
    
    float
    accept_score(MotifStateSearchNodeOP const & node) {
        return node->score();
    }
    

private:
    Points path_;
    MotifStateSearchNodeOP current_;
    Points beads_;
    int bead_size_;
    float score_, dist_, best_;
    float weight_;
    float length_;
};





typedef std::shared_ptr<MotifStateSearchScorer> MotifStateSearchScorerOP;

#endif /* defined(__RNAMake__motif_state_search_scorer__) */



























