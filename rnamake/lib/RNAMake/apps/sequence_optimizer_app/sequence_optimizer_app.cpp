//
//  sequence_optimization.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <stdexcept>

#include "base/backtrace.hpp"
#include <base/log.h>
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_topology.h"
#include "motif_data_structure/motif_graph.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"
#include "sequence_optimizer_app.hpp"
#include "util/steric_lookup.hpp"


void
SequenceOptimizerApp::setup_options() {
    add_option("design_file", String(""), base::OptionType::STRING, true);
    add_option("join_points", String(""), base::OptionType::STRING, false);
    
    
    add_option("log_level", "info", base::OptionType::STRING, false);
    add_option("out_file", "default.out", base::OptionType::STRING);
    add_option("score_file", "default.scores", base::OptionType::STRING);
    add_option("n", 1, base::OptionType::INT);
    add_option("opt", String("Internal"), base::OptionType::STRING);
    add_option("pdbs", false, base::OptionType::BOOL);
    add_option("start_design", -1, base::OptionType::INT);


}

void
SequenceOptimizerApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    base::Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, optimizer_.options(), "optimizer");
}

void
SequenceOptimizerApp::run() {
    auto log_level = base::log_level_from_str(get_string_option("log_level"));
    base::init_logging(log_level);
    LOG_INFO << "log level set to: " << get_string_option("log_level");

    auto connections = std::vector<ConnectionTemplate>();

    auto lines =base::get_lines_from_file(get_string_option("design_file"));
    LOG_INFO << lines.size() << " designs loaded from " << get_string_option("design_file");

    auto writer = ScoreFileWriter("seq_opt.scores");

    int start_design = get_int_option("start_design");

    auto design_num = -1;
    auto mg = motif_data_structure::MotifGraphOP(nullptr);
    for(auto const & l : lines) {
        design_num += 1;
        if(design_num < start_design) { continue; }

        mg = std::make_shared<motif_data_structure::MotifGraph>(l, motif_data_structure::MotifGraphStringType::MG);
        _fix_flex_helices_mtype(mg);
        mg->replace_ideal_helices();

        if(get_string_option("join_points") != "") {
            connections = _parse_connections_from_str(get_string_option("join_points"), mg);
        }
        else {
            connections = _guess_connections(mg);
        }
        auto scorer = _setup_optimizer_scorer(connections, mg);
        auto sols = optimizer_.get_optimized_sequences(mg, scorer);

        LOG_INFO << "solution found for design_num: " << design_num << " with score: " << sols[0]->dist_score;

        writer.write(design_num, 0, sols[0]->dist_score, sols[0]->sequence);
    }
}

std::vector<ConnectionTemplate>
SequenceOptimizerApp::_parse_connections_from_str(
        String const & connection_str,
        motif_data_structure::MotifGraphOP mg) {
    auto connections = std::vector<ConnectionTemplate>();
    auto connection_strs = base::split_str_by_delimiter(get_string_option("join_points"), ";");
    for (auto const & connection_str : connection_strs) {
        auto spl = base::split_str_by_delimiter(connection_str, " ");
        auto ni = std::stoi(spl[0]);
        auto ei = mg->get_node(ni)->data()->get_end_index(spl[1]);
        auto end_node = mg->get_node(ni);
        auto partner = end_node->connections()[ei]->partner(end_node->index());


        auto start = NodeIndexandEdge{ni, ei};
        auto end = NodeIndexandEdge{partner->index(), 1};
        connections_.push_back(ConnectionTemplate{start, end, "Internal"});
    }

    if (connections_.size() == 0) {
        throw SequenceOptimizerAppException("connections were supplied but were not processed correctly");
    }
    return connections;

}

std::vector<ConnectionTemplate>
SequenceOptimizerApp::_guess_connections(
        motif_data_structure::MotifGraphOP mg) {
    auto connections = std::vector<ConnectionTemplate>();
    for(auto const & n : *mg) {
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }

            auto ei1 = c->end_index(n->index());
            auto partner = c->partner(n->index());
            auto ei2 = c->end_index(partner->index());

            auto end1 = n->data()->ends()[ei1];
            auto end2 = partner->data()->ends()[ei2];
            auto diff = end1->state()->diff(end2->state());

            if(diff < 3) { continue; }

            auto nie1 = NodeIndexandEdge{n->index(), ei1};
            auto nie2 = NodeIndexandEdge{partner->index(), ei2};

            if(n->index() > partner->index()) { std::swap(nie1, nie2); }

            auto ct = ConnectionTemplate{nie1, nie2, "Internal"};
            bool found = false;
            for(auto const & ct2 : connections) {
                if(ct == ct2) { found = true; }
            }
            if(!found) {
                connections.push_back(ct);
            }
        }
    }
    return connections;
}


sequence_optimization::SequenceOptimizerScorerOP
SequenceOptimizerApp::_setup_optimizer_scorer(
        std::vector<ConnectionTemplate> const & connections,
        motif_data_structure::MotifGraphOP mg) {
    if(connections.size() == 1) {
        auto c = connections[0];
        return std::make_shared<sequence_optimization::InternalTargetScorer>(
                c.start.ni, c.start.ei, c.end.ni, c.end.ei, false);
    }
    else {
        auto sub_scorers = std::vector<sequence_optimization::SequenceOptimizerScorerOP>();
        for(auto const & c : connections) {
            bool target_an_aligned_end = false;
            if(mg->get_node(c.end.ni)->data()->block_end_add() == c.end.ei) {
                target_an_aligned_end = true;
            }
            auto scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(
                    c.start.ni, c.start.ei, c.end.ni, c.end.ei, target_an_aligned_end);
            sub_scorers.push_back(scorer);
        }
        return std::make_shared<sequence_optimization::MultiTargetScorer>(sub_scorers);
    }
}


void
SequenceOptimizerApp::_fix_flex_helices_mtype(
        motif_data_structure::MotifGraphOP mg) {
    for(auto & n : *mg) {
        if(n->data()->name().substr(0, 5) == "HELIX") {
            n->data()->mtype(util::MotifType::HELIX);
        }
    }
}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    auto app = SequenceOptimizerApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}
