//
//  design_rna.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/26/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//


#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_tools/segmenter.h"
#include "motif_data_structures/motif_topology.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "design_rna.hpp"

DesignRNAApp::DesignRNAApp() : Application(),
search_(MotifStateSearch()),
mg_(std::make_shared<MotifGraph>()),
lookup_(StericLookup())
{}

void
DesignRNAApp::setup_options() {
    
    add_option("designs", 1, OptionType::INT, false);
    add_option("seqs_per_design", 1, OptionType::INT, false);
    add_option("sol_pdbs", false, OptionType::BOOL, false);
    add_option("full_pdbs", false, OptionType::BOOL, false);
    add_option("out_file", "default.out", OptionType::STRING, false);
    add_option("score_file", "default.scores", OptionType::STRING, false);
    add_option("verbose", false, OptionType::BOOL, false);
    add_option("pdbs", false, OptionType::BOOL, false);
    
    // start from pdb
    add_option("pdb", String(""), OptionType::STRING, false);
    add_option("start_bp", String(""), OptionType::STRING, false);
    add_option("end_bp", String(""), OptionType::STRING, false);

    add_cl_options(search_.options(), "search");

}

void
DesignRNAApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, search_.options(), "search");
    search_.update_var_options();

}

void
DesignRNAApp::_setup_sterics() {
    auto beads = Points();
    for(auto & n : *mg_) {
        for(auto const & b : n->data()->beads()) {
            if(b.btype() == BeadType::PHOS) { continue; }
            beads.push_back(b.center());
            
        }
        
        for(auto const & b : n->data()->protein_beads()) { beads.push_back(b.center()); }
    }
    lookup_.add_points(beads);
}

void
DesignRNAApp::_setup_from_pdb() {
    auto struc = RM::instance().get_structure(get_string_option("pdb"), "scaffold");
    
    std::cout << "DESIGN RNA: loaded pdb from file: " << get_string_option("pdb") << std::endl;

    auto start_bp_name = get_string_option("start_bp");
    auto end_bp_name = get_string_option("end_bp");
    
    if(start_bp_name == "" || end_bp_name == "") {
        throw DesignRNAAppException(
            "must supply the name of the start_bp and end_bp using option -start_bp and "
            "-end_bp respectively when using -pdb option");
    }
    
    auto start_bps = struc->get_basepair(start_bp_name);
    auto end_bps = struc->get_basepair(end_bp_name);
    
    if(start_bps.size() == 0) {
        throw DesignRNAAppException("cannot find start basepair: " + start_bp_name);
    }
    
    if(start_bps.size() > 1) {
        throw DesignRNAAppException(
            "start basepair name: " + start_bp_name + " is not unique name, multiple basepairs "
            "have this name please pick another basepair or reformat your pdb");
    }
    
    if(end_bps.size() == 0) {
        throw DesignRNAAppException("cannot find end basepair: " + end_bp_name);
    }
    
    if(end_bps.size() > 1) {
        throw DesignRNAAppException(
            "end basepair name: " + end_bp_name + " is not unique name, multiple basepairs "
            "have this name please pick another basepair or reformat your pdb");
    }

    if(start_bps[0]->bp_type() != "cW-W") {
        throw DesignRNAAppException(
            "start basepair is not a watson and crick basepair building from non-canonical "
            "leads to bizzare effects");
    }
    
    if(end_bps[0]->bp_type() != "cW-W") {
        throw DesignRNAAppException(
            "ends basepair is not a watson and crick basepair building from non-canonical "
            "leads to bizzare effects");
    }
    
    
    auto bps = BasepairOPs{start_bps[0], end_bps[0]};
    
    auto segmenter = Segmenter();
    auto segments = segmenter.apply(struc, bps);
    
    std::cout << "DESIGN RNA: segmentation success" << std::endl;
    std::cout << "DESIGN RNA: removed rna segment=removed.pdb" << std::endl;
    std::cout << "DESIGN RNA: remaining rna segment=remaining.pdb" << std::endl;

    segments->remaining->to_pdb("remaining.pdb");
    segments->removed->to_pdb("removed.pdb");
    segments->remaining->mtype(MotifType::TWOWAY);
    
    RM::instance().register_motif(segments->remaining);
    
    //auto new_struc = std::make_shared<RNAStructure>(*struc);
    segments->remaining->block_end_add(-1);
    mg_->add_motif(segments->remaining);
    start_ = EndStateInfo{start_bp_name, 0};
    end_ = EndStateInfo{end_bp_name, 0};
}

void
DesignRNAApp::run() {
    if(get_string_option("pdb") != "") { _setup_from_pdb(); }
    
    auto start = mg_->get_node(start_.n_pos)->data()->get_basepair(start_.name)[0]->state();
    
    auto end_bp = mg_->get_node(end_.n_pos)->data()->get_basepair(end_.name)[0];
    auto end = end_bp->state();

    _setup_sterics();
    
    search_.setup(start, end);
    search_.lookup(lookup_);
    search_.set_option_value("max_solutions", 10000);
    search_.set_option_value("verbose", get_bool_option("verbose"));
    
    auto end_n_uuid = mg_->get_node(end_.n_pos)->data()->id();
    mg_->increase_level();
    
    auto so = SequenceOptimizer3D();
    so.set_option_value("verbose", get_bool_option("verbose"));
    
    std::ofstream out, sf_out;
    sf_out.open(get_string_option("score_file"));
    sf_out << "design_num,design_score,design_sequence,design_structure,opt_num,";
    sf_out << "opt_sequence,opt_score,eterna_score" << std::endl;
    
    out.open(get_string_option("out_file"));
    
    
    int design_num = 0;
    int solution_count = 0;
    while(! search_.finished()) {
        auto sol = search_.next();
        if(sol == nullptr) {
            std::cout << "DESIGN RNA: cannot find anymore solutions!" << std::endl;
            std::cout << "DESIGN RNA: generated " << design_num << " designs!" << std::endl;
            sf_out.close();
            out.close();
            exit(0);
        }
        
        auto mt = sol->to_motif_tree();

        mg_->add_motif_tree(mt, start_.n_pos, start_.name);
        auto last_node = mg_->last_node();
        mg_->add_connection(end_.n_pos, mg_->last_node()->index(), end_.name, "");
        mg_->replace_ideal_helices();
        mg_->write_pdbs();
        mg_->to_pdb("test.pdb", 1);
        
        auto c = GraphtoTree();
        auto d_mt = c.convert(mg_, nullptr, -1, last_node);
        d_mt->set_option_value("sterics", false);
        
        auto sols = so.get_optimized_sequences(d_mt, end_bp, d_mt->last_node()->index(), 1);

        if(sols.size() > 0) {
            int opt_num = 0;
            for(auto const & s : sols) {
                sf_out << design_num << "," << sol->score() << "," <<  mg_->designable_sequence() << ",";
                sf_out << opt_num << "," << s->sequence << "," << s->dist_score << "," << s->eterna_score;
                sf_out << std::endl;
                
                opt_num++;
                
                auto copy_mg = std::make_shared<MotifGraph>(*mg_);
                auto dss = copy_mg->designable_secondary_structure();
                dss->replace_sequence(s->sequence);
                copy_mg->replace_helical_sequence(dss);
                
                out << copy_mg->to_str() << std::endl;
                if(get_bool_option("pdbs")) {
                    copy_mg->to_pdb("design." + std::to_string(solution_count) + ".pdb", 1);
                    solution_count++;
                }
                
                break;
            }
            
            design_num++;
            
            if(design_num >= get_int_option("designs")) {
                std::cout << "DESIGN RNA: generated " << get_int_option("designs") << " ";
                std::cout << "design(s)! if you would like more please specify how many you ";
                std::cout << "would like with -designs #Num" << std::endl;
                sf_out.close();
                out.close();
                
                exit(0);
            }
        }
        
        mg_->remove_level(1);
        
        
        
    }
    
    sf_out.close();
    out.close();
    
    
}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);
    
    auto app = DesignRNAApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;
    
}