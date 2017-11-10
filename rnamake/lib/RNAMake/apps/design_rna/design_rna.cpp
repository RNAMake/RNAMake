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
        lookup_(StericLookup()),
        optimizer_(SequenceOptimizer3D()) {}

void
DesignRNAApp::setup_options() {

    //general options
    add_option("designs", 1, OptionType::INT, false);
    add_option("seqs_per_design", 1, OptionType::INT, false);
    add_option("out_file", "default.out", OptionType::STRING, false);
    add_option("score_file", "default.scores", OptionType::STRING, false);
    add_option("verbose", false, OptionType::BOOL, false);
    add_option("pdbs", false, OptionType::BOOL, false);
    add_option("design_pdbs", false, OptionType::BOOL, false);
    add_option("show_sections", false, OptionType::BOOL, false);
    add_option("defined_motif_path", "", OptionType::STRING, false);

    //no sequence opt
    add_option("only_ideal", false, OptionType::BOOL, false);


    // start from pdb
    add_option("pdb", String(""), OptionType::STRING, false);
    add_option("start_bp", String(""), OptionType::STRING, false);
    add_option("end_bp", String(""), OptionType::STRING, false);
    add_option("no_segment", false, OptionType::BOOL, false);

    // start from motif graph
    add_option("mg", String(""), OptionType::STRING, false);

    add_cl_options(search_.options(), "search");
    add_cl_options(optimizer_.options(), "optimizer");

}

void
DesignRNAApp::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, search_.options(), "search");
    cl_parser_.assign_options(cl_options_, optimizer_.options(), "optimizer");
    search_.update_var_options();
    optimizer_.update_var_options();

}

void
DesignRNAApp::_setup_sterics() {
    auto beads = Points();
    for (auto & n : *mg_) {
        n->data()->get_beads(n->data()->ends());
        //std::cout << n->data()->beads().size() << std::endl;
        for (auto const & b : n->data()->beads()) {
            if (b.btype() == BeadType::PHOS) { continue; }
            beads.push_back(b.center());
        }

        for (auto const & b : n->data()->protein_beads()) { beads.push_back(b.center()); }
    }
    lookup_.add_points(beads);
}

void
DesignRNAApp::_setup_from_pdb() {
    auto struc = RM::instance().get_structure(get_string_option("pdb"), "scaffold");
    std::cout << "DESIGN RNA: loaded pdb from file: " << get_string_option("pdb") << std::endl;

    auto start_bp_name = get_string_option("start_bp");
    auto end_bp_name = get_string_option("end_bp");

    if (start_bp_name == "" || end_bp_name == "") {
        throw DesignRNAAppException(
                "must supply the name of the start_bp and end_bp using option -start_bp and "
                        "-end_bp respectively when using -pdb option");
    }

    auto start_bps = struc->get_basepair(start_bp_name);
    auto end_bps = struc->get_basepair(end_bp_name);

    if (start_bps.size() == 0) {
        throw DesignRNAAppException("cannot find start basepair: " + start_bp_name);
    }

    if (start_bps.size() > 1) {
        throw DesignRNAAppException(
                "start basepair name: " + start_bp_name + " is not unique name, multiple basepairs "
                        "have this name please pick another basepair or reformat your pdb");
    }

    if (end_bps.size() == 0) {
        throw DesignRNAAppException("cannot find end basepair: " + end_bp_name);
    }

    if (end_bps.size() > 1) {
        throw DesignRNAAppException(
                "end basepair name: " + end_bp_name + " is not unique name, multiple basepairs "
                        "have this name please pick another basepair or reformat your pdb");
    }

    if (start_bps[0]->bp_type() != "cW-W") {
        throw DesignRNAAppException(
                "start basepair is not a watson and crick basepair building from non-canonical "
                        "leads to bizzare effects");
    }

    if (end_bps[0]->bp_type() != "cW-W") {
        throw DesignRNAAppException(
                "ends basepair is not a watson and crick basepair building from non-canonical "
                        "leads to bizzare effects");
    }


    auto bps = BasepairOPs{start_bps[0], end_bps[0]};

    if (get_bool_option("no_segment")) {
        RM::instance().add_motif(get_string_option("pdb"), "scaffold");
        auto m = RM::instance().motif("scaffold", "", end_bp_name);
        //RM::instance().register_motif(m);
        mg_->add_motif(m);
        mg_->write_pdbs();
        start_ = EndStateInfo{start_bp_name, 0};
        end_ = EndStateInfo{end_bp_name, 0};
        return;
    }
    auto segmenter = Segmenter();
    auto segments = segmenter.apply(struc, bps);

    std::cout << "DESIGN RNA: segmentation success" << std::endl;
    std::cout << "DESIGN RNA: removed rna segment=removed.pdb" << std::endl;
    std::cout << "DESIGN RNA: remaining rna segment=remaining.pdb" << std::endl;

    if (get_bool_option("show_sections")) {
        segments->remaining->to_pdb("remaining.pdb");
        segments->removed->to_pdb("removed.pdb");
    }
    segments->remaining->mtype(MotifType::TWOWAY);

    RM::instance().register_motif(segments->remaining);

    //auto new_struc = std::make_shared<RNAStructure>(*struc);
    segments->remaining->block_end_add(-1);
    mg_->add_motif(segments->remaining);
    start_ = EndStateInfo{start_bp_name, 0};
    end_ = EndStateInfo{end_bp_name, 0};
}

void
DesignRNAApp::_setup_from_mg() {

    auto lines = get_lines_from_file(get_string_option("mg"));
    mg_ = std::make_shared<MotifGraph>(lines[0], MotifGraphStringType::MG);
    mg_->set_option_value("sterics", false);
    auto spl = split_str_by_delimiter(lines[1], " ");
    start_ = EndStateInfo{spl[1], std::stoi(spl[0])};
    spl = split_str_by_delimiter(lines[2], " ");
    end_ = EndStateInfo{spl[1], std::stoi(spl[0])};

    try {mg_->get_node(start_.n_pos); }
    catch(...) {
        throw DesignRNAAppException("node: " + std::to_string(start_.n_pos) + " does not exist!");
    }

    try {mg_->get_node(end_.n_pos); }
    catch(...) {
        throw DesignRNAAppException("node: " + std::to_string(end_.n_pos) + " does not exist!");
    }


}

std::shared_ptr<MSS_Path>
DesignRNAApp::_setup_path() {
    auto spl = split_str_by_delimiter(get_string_option("defined_motif_path"), ",");
    auto selector = std::make_shared<MSS_Path>();
    auto i = 0;
    for(auto const & name : spl) {
        if(name.length() < 2) { continue; }
        if(name == "ideal_helices_min" || name == "unique_twoway" || name == "tcontact" ||
                name == "twoway") {
            selector->add(name);
        }
        else {
            auto m = RM::instance().motif(name);
            auto motif_states = MotifStateOPs();
            auto scores = Floats();

            for(auto const & end : m->ends()) {
                auto m_new = RM::instance().motif(name, "", end->name());
                m_new->new_res_uuids();
                //m_new->ends()[1]->flip();
                //std::cout << m_new->name() << std::endl;
                motif_states.push_back(m_new->get_state());
                scores.push_back(1);
            }
            auto mse = std::make_shared<MotifStateEnsemble>(motif_states, scores);
            selector->add("", mse);
        }
        if(i > 0) {
            selector->connect(i-1, i);
        }
        i++;
    }
    return selector;
}

void
DesignRNAApp::run() {
    if (get_string_option("pdb") != "") { _setup_from_pdb(); }
    if (get_string_option("mg") != "") { _setup_from_mg(); }

    //mg_->write_pdbs();
    auto start = mg_->get_node(start_.n_pos)->data()->get_basepair(start_.name)[0]->state();
    auto end_bp = mg_->get_node(end_.n_pos)->data()->get_basepair(end_.name)[0];
    auto end = end_bp->state();

    _setup_sterics();

    if(get_string_option("defined_motif_path") != "") {
        std::shared_ptr<MSS_Path> custom_selector = _setup_path();
        search_.selector(custom_selector);
    }
    /*auto selector = std::make_shared<MSS_Path>();
    //selector->add("ideal_helices_min");
    selector->add("unique_twoway");
    selector->add("ideal_helices_min");
    selector->add("tcontact");
    selector->add("ideal_helices_min");
    //selector->add("unique_twoway");
    //selector->add("ideal_helices_min");
    selector->connect(0, 1);
    selector->connect(1, 2);
    selector->connect(2, 3);
    //selector->connect(3, 4);
    //selector->connect(4, 5);
    //selector->connect(5, 6);
    */
    search_.setup(start, end);
    search_.lookup(lookup_);
    search_.set_option_value("max_solutions", 10000000);
    search_.set_option_value("verbose", get_bool_option("verbose"));
    search_.set_option_value("accept_score", 10);

    auto end_n_uuid = mg_->get_node(end_.n_pos)->data()->id();
    mg_->increase_level();

    optimizer_.set_option_value("verbose", get_bool_option("verbose"));
    //optimizer_.set_option_value("cutoff", 5.0f);

    std::ofstream out, sf_out;
    sf_out.open(get_string_option("score_file"));
    sf_out << "design_num,design_score,design_sequence,design_structure,motifs_uses,opt_num,";
    sf_out << "opt_sequence,opt_score,eterna_score" << std::endl;

    out.open(get_string_option("out_file"));

    int design_num = 0;
    int solution_count = 0;
    while (!search_.finished()) {
        auto sol = search_.next();
        if (sol == nullptr || sol->path().size() == 0) {
            std::cout << "DESIGN RNA: cannot find anymore solutions!" << std::endl;
            std::cout << "DESIGN RNA: generated " << design_num << " designs!" << std::endl;
            sf_out.close();
            out.close();
            exit(0);
        }


        auto mt = sol->to_motif_tree();
        mg_->add_motif_tree(mt, start_.n_pos, start_.name);
        for(auto const & n : *mt) {
            std::cout << n->data()->name() << std::endl;
        }
        //mg_->write_pdbs();

        mg_->add_connection(end_.n_pos, mg_->last_node()->index(), end_.name, "");
        mg_->replace_ideal_helices();

        auto end_node = mg_->get_node(end_.n_pos);
        auto end_i = end_node->data()->get_end_index(end_.name);
        auto partner = end_node->connections()[end_i]->partner(end_node->index());

        if (get_bool_option("only_ideal")) {
            auto motif_names = String("");
            for (auto const & n : *mt) {
                motif_names += n->data()->name() + ";";
            }

            sf_out << design_num << "," << sol->score() << ","; //<< mg_->designable_sequence();
            //sf_out << "," << mg_->secondary_structure()->dot_bracket() << ",";
            sf_out << motif_names << std::endl;

            out << mg_->to_str() << std::endl;

            design_num++;

            if (get_bool_option("pdbs")) {
                mg_->to_pdb("design." + std::to_string(solution_count) + ".pdb", 1);
                solution_count++;
                //mg_->write_pdbs();
                //exit(0);
            }

            if(get_bool_option("design_pdbs")) {
                mt->to_pdb("design." + std::to_string(solution_count) + ".pdb", 1);
                solution_count++;
            }

            if (design_num >= get_int_option("designs")) {
                std::cout << "DESIGN RNA: generated " << get_int_option("designs") << " ";
                std::cout << "design(s)! if you would like more please specify how many you ";
                std::cout << "would like with -designs #Num" << std::endl;
                sf_out.close();
                out.close();

                exit(0);
            }

            mg_->remove_level(1);
            continue;
        }

        auto scorer = std::make_shared<ExternalTargetScorer>(end_bp->state(), partner->index(), 1);
        auto sols = optimizer_.get_optimized_sequences(mg_, scorer);

        if (sols.size() > 0) {
            int opt_num = 0;
            for (auto const & s : sols) {
                sf_out << design_num << "," << sol->score() << "," << mg_->designable_sequence() << ",";
                sf_out << opt_num << "," << s->sequence << "," << s->dist_score << "," << s->eterna_score;
                sf_out << std::endl;

                opt_num++;

                try {
                    auto copy_mg = std::make_shared<MotifGraph>(*mg_);
                    auto dss = copy_mg->designable_secondary_structure();
                    dss->replace_sequence(s->sequence);
                    copy_mg->replace_helical_sequence(dss);
                    out << copy_mg->to_str() << std::endl;

                    if (get_bool_option("pdbs")) {
                        copy_mg->to_pdb("design." + std::to_string(solution_count) + ".pdb", 1);
                        solution_count++;
                    }
                } catch (std::runtime_error const & e) {
                    std::cout << e.what() << std::endl;
                }

                break;
            }

            design_num++;

            if (design_num >= get_int_option("designs")) {
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


int main(int argc, const char *argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);

    //load extra motifs being used
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    //RM::instance().add_motif(base_path+"GAAA_tetraloop");
    //RM::instance().add_motif(base_path+"GGAA_tetraloop");
    RM::instance().add_motif(base_path+"ATP_apt.pdb");

    auto app = DesignRNAApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}
