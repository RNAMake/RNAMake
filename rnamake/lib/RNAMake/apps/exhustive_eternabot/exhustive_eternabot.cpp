//
//  exhustive_eternabot.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/5/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <fstream>
#include "exhustive_eternabot.h"

#include "util/cartesian_product.h"
#include "eternabot/scorer.h"
#include "secondary_structure/secondary_structure_parser.h"

CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    CommandLineOptions cl_opts;
    cl_opts.add_option("seq", String("NNNNAAAANNNN"), OptionType::STRING, false);
    cl_opts.add_option("ss",  String("((((....))))"), OptionType::STRING, false);
    cl_opts.add_option("out", String("4bp.out"), OptionType::STRING, false);

    
    cl_opts.parse_command_line(argc, argv);
    return cl_opts;
}

void
ExhustiveEternabot::setup(CommandLineOptions const & opts) {
    
    out_name_ = opts.get_string("out");
    auto parser = sstruct::SecondaryStructureParser();
    p_ = parser.parse_to_pose(opts.get_string("seq"), opts.get_string("ss"));
    
    for(auto const & bp : p_->basepairs()) {
        if     (bp->res1()->res_type() == -1 && bp->res2()->res_type() == -1) {
            enumerated_bps_.push_back(bp); }
        else if(bp->res1()->res_type() == -1 || bp->res2()->res_type() == -1) {
            throw std::runtime_error("only one residue in a basepair is marked as N, this is not allowed");
        }
    }
    
    std::cout << "num of designed bps " << enumerated_bps_.size() << std::endl;
    
    pairs_ = std::vector<Strings>();
    pairs_.push_back(Strings{"A", "U"});
    pairs_.push_back(Strings{"U", "A"});
    pairs_.push_back(Strings{"C", "G"});
    pairs_.push_back(Strings{"G", "C"});

}


void
ExhustiveEternabot::run() {
    
    auto all_pairs = std::vector<std::vector<Strings>>(enumerated_bps_.size());
    for(int i = 0; i < all_pairs.size(); i++) {
        all_pairs[i] = pairs_;
    }
    
    auto pair_iterator = CartesianProduct<Strings>(all_pairs);
    auto current = std::vector<Strings>();
    
    auto scorer = eternabot::Scorer();
    scorer.setup(p_);
    float score = 0;
    float best = 0;
    auto best_seq = p_->sequence();
    auto best_helix = String("");
    
    auto out = std::ofstream(out_name_);
    out << "sequence,score,a_basic_test,clean_plot,berex_test,num_of_yellow,direction_of_gc," << std::endl;
    
    int i = 0;
    auto str = String("");
    while (!pair_iterator.end()) {
        current = pair_iterator.next();
        i = 0;
        str = "";
        for(auto const & p : current) {
            str += p[0];
            enumerated_bps_[i]->res1()->name(p[0]);
            enumerated_bps_[i]->res2()->name(p[1]);
            i++;
        }
        score =  scorer.score_secondary_structure(p_);
        if(score > best) {
            best = score;
            best_seq = p_->sequence();
            best_helix = str;
        }
        out << str << "," << score << ",";
        for(auto const & v : scorer.scores()) {
            out << v << ",";
        }
        out << std::endl;
        
    }
    
    std::cout << best << " " << best_helix << " " << best_seq << std::endl;
    out.close();
}



int main(int argc, const char * argv[]) {
    auto cmd_opts = parse_command_line(argc, argv);
    
    auto app = ExhustiveEternabot();
    app.setup(cmd_opts);
    app.run();
    
}


