//
//  sequence_optimization.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef sequence_optimizer_app_hpp
#define sequence_optimizer_app_hpp

#include <stdio.h>

#include "base/application.hpp"
#include "sequence_optimization/sequence_optimizer_3d.hpp"

class SequenceOptimizerAppException : public std::runtime_error {
public:
    SequenceOptimizerAppException(
            String const & message):
            std::runtime_error(message)
    {}
};


struct NodeIndexandEdge {
    int ni; //node index
    int ei; //end index

    inline
    bool
    operator == (
            NodeIndexandEdge const & nie) const {
        return(nie.ni == ni && nie.ei == ei);
    }
};

struct ConnectionTemplate {
    NodeIndexandEdge start;
    NodeIndexandEdge end;
    String type;

    inline
    bool
    operator == (
            ConnectionTemplate const & c) const {
        return (c.start == start && c.end == end && type == c.type);
    }
};

class ScoreFileWriter {
public:


    ScoreFileWriter(
            String const & score_file_path) {
        out_ = std::make_shared<std::ofstream>();
        out_->open(score_file_path);
        *out_ << "design_num,opt_num,opt_score,sequence" << std::endl;
    }

    ~ScoreFileWriter() {
        out_->close();
    }

public:
    void
    write(
            int design_num,
            int opt_num,
            float opt_score,
            String sequence) {
        *out_ << design_num << "," << opt_num << "," << opt_score << "," << sequence << std::endl;
    }

private:
    std::shared_ptr<std::ofstream> out_;

};


class SequenceOptimizerApp : public base::Application {
public:
    SequenceOptimizerApp() : base::Application(),
    optimizer_(sequence_optimization::SequenceOptimizer3D()) {}
    
    ~SequenceOptimizerApp() {}
    
public:
    
    void
    setup_options();
    
    void
    parse_command_line(
        int,
        const char **);
    
    void
    run();

private:
    std::vector<ConnectionTemplate>
    _parse_connections_from_str(
            String const &,
            motif_data_structure::MotifGraphOP);

    std::vector<ConnectionTemplate>
    _guess_connections(
            motif_data_structure::MotifGraphOP);

    sequence_optimization::SequenceOptimizerScorerOP
    _setup_optimizer_scorer(
            std::vector<ConnectionTemplate> const &,
            motif_data_structure::MotifGraphOP);

    void
    _fix_flex_helices_mtype(
            motif_data_structure::MotifGraphOP);

    
private:
    std::vector<ConnectionTemplate> connections_;
    sequence_optimization::SequenceOptimizer3D optimizer_;
    
    
};



#endif /* sequence_optimizer_hpp */
