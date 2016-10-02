//
//  sequence_optimizer_3d.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/30/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef sequence_optimizer_3d_hpp
#define sequence_optimizer_3d_hpp

#include <stdio.h>

#include "base/types.h"
#include "base/option.h"
#include "util/random_number_generator.h"
#include "eternabot/sequence_designer.h"
#include "motif_data_structures/motif_tree.h"

class SequenceOptimizer3D : public OptionClass {
public:
    
    SequenceOptimizer3D();
    
    ~SequenceOptimizer3D() {}
    
private:
    
    struct OptimizedSequence {
        String sequence;
        float dist_score, eterna_score;
    };
    
    typedef std::shared_ptr<OptimizedSequence> OptimizedSequenceOP;
    typedef std::vector<OptimizedSequenceOP> OptimizedSequenceOPs;
    
    struct DesignableBP {
        sstruct::BasepairOP bp;
        String last_state;
        std::shared_ptr<Uuid> m_id_bot, m_id_top;
    };
    
    typedef std::shared_ptr<DesignableBP> DesignableBPOP;
    typedef std::vector<DesignableBPOP> DesignableBPOPs;
    
public:
    
    OptimizedSequenceOPs
    get_optimized_sequences(
        MotifTreeOP const &,
        BasepairOP const &,
        int,
        int);
    
public: //option wrappers
    
    inline
    Options &
    options() { return options_; }
    
    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
    void
    update_var_options();
    
protected:
    void
    setup_options();
    
private:
    Options options_;
    eternabot::SequenceDesigner sequence_designer_;
    RandomNumberGenerator rng_;
    
    
};

#endif /* sequence_optimizer_3d_hpp */
