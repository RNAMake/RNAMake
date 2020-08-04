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
        //*out_ << "design_num,opt_num,opt_score,sequence" << std::endl;

        *out_ << "design_num,design_score,design_sequence,design_structure,motifs_uses,opt_num,";
        *out_ << "opt_sequence,opt_score,eterna_score" << std::endl;
    }

    ~ScoreFileWriter() {
        out_->close();
    }

public:
    void
    write(
            motif_data_structure::MotifGraphOP mg,
            sequence_optimization::OptimizedSequenceOP sequence_opt_sol,
            int design_num,
            int sequence_opt_num) {
        _get_motif_names(mg);
        *out_ << design_num << "," << sequence_opt_sol->dist_score << "," << mg->designable_sequence() << ",";
        *out_ << mg->dot_bracket() << "," << motif_names_ << "," << sequence_opt_num << "," << sequence_opt_sol->sequence;
        *out_ << "," << sequence_opt_sol->dist_score << "," << sequence_opt_sol->eterna_score << std::endl;
    }

private:
    String const &
    _get_motif_names(
            motif_data_structure::MotifGraphOP mg) {
        motif_names_ = "";
        for (auto const & n : *mg) {
            motif_names_ += n->data()->name() + ";";
        }
        return motif_names_;
    }

private:
    std::shared_ptr<std::ofstream> out_;
    String motif_names_;

};


class SequenceOptimizerApp : public base::Application {
public:
    struct Parameters {
        String design_file, log_level, out_file, score_file, join_points;
        int n, start_design;
        float cutoff;
        bool return_lowest, pdbs;
    };


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
    Parameters parameters_;
    
};



#endif /* sequence_optimizer_hpp */
