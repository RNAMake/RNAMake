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
        lengths_ = Floats(1000);
        seen_ = Ints(path_.size());
        bead_arrays_ = std::vector<Points>(1000);
        bead_lengths_ = Floats(path_.size());
        bead_size_ = 0;
        weight_ = 2.0;
        path_length_ = 0;
        end_direction_ =  path[path_.size()-2] - path_.back();
        end_direction_ = end_direction_.normalize();
        current_direction_ = Vector();
        bead_lengths_[0] = 0;
        for(int i = 1; i < path_.size(); i++) {
            path_length_ += path_[i-1].distance(path_[i]);
            bead_lengths_[i] = path_length_;
            seen_[i] = 0;
        }
        
        
    }
    
    inline
    float
    _score(MotifStateSearchNodeOP const & node) {
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
                if(best_ < 7) { break; }
            }
            if(best_ > 7) {
                //score_ += 1 - (1000 - best_*5)/1000;
                score_ += best_-7;
            }
            else {
                if(sum < i) {
                    score_ += (float)(i - sum)*10;
                }
                sum += 1;
            }
            i++;
        }
        
        //return length_*0.01 + node->level()*0.50 + score_;
        //score_ += length_*0.01 + score_;
        if(length_ > path_length_) {
            score_ += (length_ - path_length_);
        }
        score_ += path_.back().distance(beads_[bead_size_-1]);
    

        return score_;
        
    }
    
    inline
    float
    _score_2(MotifStateSearchNodeOP const & node) {
        current_ = node;
        bead_size_ = 0;
        weight_ = 1.0;
        length_ = 0;
        int bead_count = 0;
        while(1) {
            /*for(auto const & b : current_->cur_state()->beads()) {
             beads_[bead_size_] = b;
             bead_size_++;
             }*/
            bead_arrays_[bead_size_] = current_->cur_state()->beads();
            
            length_ += current_->cur_state()->end_states()[1]->d().distance(current_->cur_state()->end_states()[0]->d());
            lengths_[bead_size_] = current_->cur_state()->end_states()[1]->d().distance(current_->cur_state()->end_states()[0]->d());
            bead_size_++;
            
            current_ = current_->parent();
            if(current_ == nullptr) { break; }
        }
        
        if(length_ > path_length_*2) {
            return 1000000;
        }
        
        score_ = 0;
        int sum = 0;
        int i = 0, j = -1;
        
        int pos = bead_size_-2;
        int b_pos = 0;
        int b_end = 0;
        float diff, best_diff = 10000;
        int best_b_pos = 0;
        float avg_diff = 0;
        float current_length = lengths_[pos];
        
        for(auto const & p : path_) {
            j++;
            
            for(i = bead_arrays_.size(); i > 0; i--) {
                
            }
            
           
            b_end = best_b_pos;
            //std::cout << pos << " " << b_pos << " " << b_end << " " << current_length << " " << bead_lengths_[best_b_pos] << " " << path_.size() << std::endl;
            
            for(i = b_pos; i < b_end; i++) {
                best_ = 10000;
                for(auto const & b : bead_arrays_[pos]) {
                    dist_ = path_[i].distance(b);
                    if(best_ > dist_) {
                        best_ = dist_;
                    }
                    // if(best_ < 10) { break; }
                }
                /*if(best_ > 10) {
                 //score_ += 1 - (1000 - best_*5)/1000;
                 score_ += best_-10;
                 }*/
                score_ += best_;
                
            }
            pos -= 1;
            if(pos > -1) {
                current_length += lengths_[pos];
            }
            weight_ += 1;
            
        }
        
        //avg_diff /= path_.size();
        score_ += (path_.size() - b_end)*100;
        
        return score_;
        
    }

    
    inline
    float
    _score_3(MotifStateSearchNodeOP const & node) {
        current_ = node;
        bead_size_ = 0;
        weight_ = 1.0;
        length_ = 0;
        int bead_count = 0;
        while(1) {
            /*for(auto const & b : current_->cur_state()->beads()) {
             beads_[bead_size_] = b;
             bead_size_++;
             }*/
            bead_arrays_[bead_size_] = current_->cur_state()->beads();
            
            length_ += current_->cur_state()->end_states()[1]->d().distance(current_->cur_state()->end_states()[0]->d());
            lengths_[bead_size_] = current_->cur_state()->end_states()[1]->d().distance(current_->cur_state()->end_states()[0]->d());
            bead_size_++;
            
            current_ = current_->parent();
            if(current_ == nullptr) { break; }
        }
        
        score_ = 0;
        int sum = 0;
        int i = 0, j = -1;
        
        int pos = bead_size_-2;
        int b_pos = 0;
        int b_end = 0;
        float diff, best_diff = 10000;
        int best_b_pos = 0;
        float avg_diff = 0;
        float current_length = lengths_[pos];
        
        for(i = 0; i < seen_.size(); i++) {
            seen_[i] = 500;
        }
        
        int last_seen_pos_ = 0;
        int total = 0;
        int k = 0;
        
        for(i = pos; i >= 0; i--) {
            
            for(auto const & b : bead_arrays_[i]) {
                best_ = 10000;
                best_b_pos = 0;
                for(j = 0; j < path_.size(); j++ ) {
                    dist_ = path_[j].distance(b);
                    if(dist_ < best_) {
                        best_b_pos = j;
                        best_ = dist_;
                    }
                    if(seen_[j] > dist_) {
                        seen_[j] = dist_;
                    }
                    
                }
                
                total = 0;
                for(k = 0; k < best_b_pos; k++){
                    if(seen_[k] > 10) { total += 1; }
                }
                
                //if(best_ < 10) {
                //    score_ += (best_b_pos - total);
                //}
                //else {
                score_ += (total /(float)seen_.size() )*10;
                //}
                
            }
            
        }
        
        int count = 0;
        
        for(auto const & spos : seen_) {
            count = spos;
            /*if(count < -15) {
                count = -15;
            }*/
            score_+= count;
        }
        
        score_ += fabsf(length_ - bead_lengths_.back()) * 5;
        
        //score_ += seen_.back()*10;
        
        return score_;
        
    }

    
    inline
    float
    score(MotifStateSearchNodeOP const & node) {
        current_ = node;
        bead_size_ = 0;
        weight_ = 1.0;
        length_ = 0;
        int bead_count = 0;
        while(1) {
            /*for(auto const & b : current_->cur_state()->beads()) {
             beads_[bead_size_] = b;
             bead_size_++;
             }*/
            bead_arrays_[bead_size_] = current_->cur_state()->beads();
            
            length_ += current_->cur_state()->end_states()[1]->d().distance(current_->cur_state()->end_states()[0]->d());
            lengths_[bead_size_] = current_->cur_state()->end_states()[1]->d().distance(current_->cur_state()->end_states()[0]->d());
            bead_size_++;
            
            current_ = current_->parent();
            if(current_ == nullptr) { break; }
        }
        
        score_ = 0;
        int sum = 0;
        int i = 0, j = -1;
        
        int pos = bead_size_-2;
        int b_pos = 0;
        int b_end = 0;
        float diff, best_diff = 100000000;
        int best_b_pos = 0;
        float avg_diff = 0;
        float current_length = lengths_[pos];
        
        for(i = 0; i < seen_.size(); i++) {
            seen_[i] = 500;
        }
        
        int last_seen_pos_ = 0;
        int total = 0;
        int k = 0;
        
        for(i = 0; i < pos; i++) {
            
            for(auto const & b : bead_arrays_[i]) {
                best_ = 10000000;
                best_b_pos = 0;
                for(j = 0; j < path_.size(); j++ ) {
                    dist_ = path_[j].distance(b);
                    if(dist_ < best_) {
                        best_b_pos = j;
                        best_ = dist_;
                    }
                    if(seen_[j] > dist_) {
                        seen_[j] = dist_;
                    }
                    
                }
                
            }
            
            
        }
        
        int count = 0;
        
        for(auto const & spos : seen_) {
            count = spos;
            /*if(count < -15) {
             count = -15;
             }*/
            score_+= count;
        }
        
        current_direction_.x(node->cur_state()->end_states()[1]->r().zx());
        current_direction_.y(node->cur_state()->end_states()[1]->r().zy());
        current_direction_.z(node->cur_state()->end_states()[1]->r().zz());

        
        score_ += path_.back().distance(node->cur_state()->end_states()[1]->d())*10;
        score_ += end_direction_.distance(current_direction_.normalize())*50;
        
        if(length_ > bead_lengths_.back()) {
            score_ += (length_ - bead_lengths_.back()) * 2;
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
    Vector end_direction_, current_direction_;
    std::vector<Points> bead_arrays_;
    Floats bead_lengths_, lengths_;
    Ints seen_;
    int bead_size_;
    float score_, dist_, best_;
    float weight_;
    float length_;
    float path_length_;
};





typedef std::shared_ptr<MotifStateSearchScorer> MotifStateSearchScorerOP;

#endif /* defined(__RNAMake__motif_state_search_scorer__) */



























