//
//  simulate_tectos.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include "base/cl_option.h"
#include "base/settings.h"
#include "math/hashing.h"
#include "structure/residue_type_set_manager.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_tree.h"
#include "simulate_tectos_devel.h"
#include "math/euler.h"

// loggers               ////////////////////////////////////////////////////////////////////////////

class SimulateTectosRecordAllLogger : public thermo_fluctuation::ThermoFlucSimulationLogger {
public:
    SimulateTectosRecordAllLogger(
            String const & fname,
            Strings const & motif_names,
            Strings const & ggaa_ttr_end_names,
            Strings const & gaaa_ttr_end_names):
            thermo_fluctuation::ThermoFlucSimulationLogger(fname),
            motif_names_(motif_names),
            ggaa_ttr_end_names_(ggaa_ttr_end_names),
            gaaa_ttr_end_names_(gaaa_ttr_end_names) {
    }

public:

    void
    log(
            motif_data_structure::MotifStateTreeOP const & mst,
            float score) {
        int i = 0;
        for(auto const & n : *mst) {
            i++;
            if(i == 1) { continue; } // ignore first state
            if(n->data()->cur_state->end_states().size() == 2) {
                out_ << vector_to_str(n->data()->get_end_state(0)->d()) << ",";
                out_ << matrix_to_str(n->data()->get_end_state(0)->r()) << ",";
                out_ << n->data()->cur_state->name() << ",";
                out_ << n->data()->cur_state->end_ids()[0];
            }
            else {
                int i = 0;
                int size = n->data()->cur_state->end_states().size();
                for(auto const & end_state : n->data()->cur_state->end_states()) {
                    i++;
                    out_ << matrix_to_str(end_state->r()) << ",";
                    out_ << vector_to_str(end_state->d());
                    if(i != size) { out_ << ","; }
                }

            }
            out_ << ",";
        }
        out_ << score << std::endl;

    }

protected:
    void
    _output_header(
            motif_data_structure::MotifStateTreeOP const & mst) {
        for(auto const & name : motif_names_) {
            if(name.length() < 4) {
                this->out_ << name + "_d," <<  name + "_r," << name + "_bp," << name << "_ei";
                this->out_ << ",";
            }
            else {
                if(name == "gaaa_ttr") { _output_ttr_header_section(name, gaaa_ttr_end_names_); }
                if(name == "ggaa_ttr") { _output_ttr_header_section(name, ggaa_ttr_end_names_); }
                this->out_ << ",";

            }

        }
        this->out_ << "score\n";

        outputed_header_ = true;
    }

    void
    _output_ttr_header_section(
            String const & name,
            Strings const & end_names) {
        int i = 0;
        for(auto const & end_name : end_names) {
            this->out_ << name + "_" + end_name + "_r," <<  name + "_" + end_name + "_d";
            i++;
            if(i != end_names.size()) { this->out_ << ","; }
        }

    }

private:
    Strings motif_names_;
    Strings ggaa_ttr_end_names_;
    Strings gaaa_ttr_end_names_;

};

class SimulateTectosRecord6D : public thermo_fluctuation::ThermoFlucSimulationLogger {
public:
    SimulateTectosRecord6D(
            String const & fname,
            String const & record_constraints,
            structure::BasepairOP ref_bp ):
            thermo_fluctuation::ThermoFlucSimulationLogger(fname),
            ref_bp_(ref_bp) {
        _setup_constraints();
        _parse_constraints(record_constraints);
        ref_r_t_ = ref_bp_->r().transposed();
    }

public:
    void
    log(
            motif_data_structure::MotifStateTreeOP const & mst,
            float score) {

        auto end_state_1 = mst->get_node(ni1_)->data()->get_end_state(ei1_);
        auto end_state_2 = mst->get_node(ni2_)->data()->get_end_state(ei2_);

        r1_ = end_state_1->r();
        d1_ = end_state_1->d();

        r2_ = end_state_2->r();
        d2_ = end_state_2->d();

        dot(ref_r_t_, r1_, rot_);
        rot_t_.unitarize();
        rot_t_ = rot_.transposed();


        dot(r1_, rot_t_, r1_trans_);
        dot(r2_, rot_t_, r2_trans_);

        dot(r1_trans_.transposed(), r2_trans_, r_);
        r_.unitarize();

        auto euler = math::Vector();
        math::calc_euler(r_, euler);

        d_ = d2_ - d1_;

        for(int i = 0; i < 3; i++) {
            euler[i] = euler[i]*180/M_PI;
            if(euler[i] > 180) {
                euler[i] -= 360;
            }
        }

        value_[0] = d_[0];
        value_[1] = d_[1];
        value_[2] = d_[2];
        value_[3] = euler[0];
        value_[4] = euler[1];
        value_[5] = euler[2];

        for(int i = 0; i < 6; i++) {
            if(constraints_[i][0] > value_[i] || value_[i] > constraints_[i][1]) { return; }
        }

        dist_ = sqrt(value_[0]*value_[0] + value_[1]*value_[1] + value_[2]*value_[2]);
        if(dist_ > 7) { return; }


        out_ << d_[0] << "," << d_[1] << "," << d_[2] << "," << euler[0] << "," << euler[1] << "," << euler[2] << ",";
        out_ << score << std::endl;


    }

protected:
    void
    _output_header(
            motif_data_structure::MotifStateTreeOP const & mst) {
        this->out_ << "x,y,z,a,b,g,score" << std::endl;
        outputed_header_ = true;

    }

private:
    void
    _parse_constraints(
            String const & constraints) {

        auto spl =  base::split_str_by_delimiter(constraints, ";");
        for(auto const & s : spl) {
            auto spl2 = base::split_str_by_delimiter(s, ",");
            if(spl2.size() != 3) { throw SimulateTectosAppException("invalid record constraint: " + s); }
            auto pos = _parse_constraint_position(spl2[0]);
            auto lower = std::stod(spl2[1]);
            auto upper = std::stod(spl2[2]);
            constraints_[pos] = math::Real2{lower, upper};
        }
    }

    int
    _parse_constraint_position(
            String const & constraint_name) {
        if     (constraint_name == "x") { return 0; }
        else if(constraint_name == "y") { return 1; }
        else if(constraint_name == "z") { return 2; }
        else if(constraint_name == "a") { return 3; }
        else if(constraint_name == "b") { return 4; }
        else if(constraint_name == "g") { return 5; }
        else {
            throw SimulateTectosAppException("unknown constraint name: " + constraint_name);
        }

    }

    void
    _setup_constraints() {
        for(int i = 0; i < 3; i++) { constraints_[i] = math::Real2{ -100, 100}; }
        for(int i = 3; i < 6; i++) { constraints_[i] = math::Real2{ -180, 180}; }
    }

private:
    structure::BasepairOP ref_bp_;
    math::Matrix rot_, rot_t_, ref_r_t_, r_, r1_, r2_, r1_trans_, r2_trans_;
    math::Point d1_, d2_, d_;
    math::Real6 value_;
    std::array<math::Real2, 6> constraints_;
    double dist_;


};

class SimulateTectosRecord6DHistogram : public thermo_fluctuation::ThermoFlucSimulationLogger {
public:
    SimulateTectosRecord6DHistogram(
            String const & fname,
            String const & constraints,
            structure::BasepairOP ref_bp ):
            thermo_fluctuation::ThermoFlucSimulationLogger("out.out"),
            ref_bp_(ref_bp),
            histo_(math::SixDHistogram(math::BoundingBox(), math::Real6{0.1, 0.1, 0.1, 0.1, 0.1, 0.1})),
            file_name_(fname) {
        ref_r_t_ = ref_bp_->r().transposed();
        _setup_constraints();
        _parse_constraints(constraints);
        auto lower = math::Point(constraints_[0][0], constraints_[1][0], constraints_[2][0]);
        auto upper = math::Point(constraints_[0][1], constraints_[1][1], constraints_[2][1]);
        auto bb = math::BoundingBox(lower, upper);
        auto bin_widths = math::Real6{0.25, 0.25, 0.25, 5.0, 5.0, 5.0};
        histo_ = math::SixDHistogram(bb, bin_widths);
    }

public:
    void
    log(
            motif_data_structure::MotifStateTreeOP const & mst,
            float score) {

        auto end_state_1 = mst->get_node(ni1_)->data()->get_end_state(ei1_);
        auto end_state_2 = mst->get_node(ni2_)->data()->get_end_state(ei2_);

        r1_ = end_state_1->r();
        d1_ = end_state_1->d();

        r2_ = end_state_2->r();
        d2_ = end_state_2->d();

        dot(ref_r_t_, r1_, rot_);
        rot_t_.unitarize();
        rot_t_ = rot_.transposed();


        dot(r1_, rot_t_, r1_trans_);
        dot(r2_, rot_t_, r2_trans_);

        dot(r1_trans_.transposed(), r2_trans_, r_);
        r_.unitarize();

        auto euler = math::Vector();
        math::calc_euler(r_, euler);

        d_ = d2_ - d1_;
        for(int i = 0; i < 3; i++) {
            euler[i] = euler[i]*180/M_PI;
            euler[i] += 180;
        }
        values_[0] = d_[0];
        values_[1] = d_[1];
        values_[2] = d_[2];
        values_[3] = euler[0];
        values_[4] = euler[1];
        values_[5] = euler[2];

        for(int i = 0; i < 6; i++) {
            if(constraints_[i][0] > values_[i] || values_[i] > constraints_[i][1]) { return; }
        }

        dist_ = sqrt(values_[0]*values_[0] + values_[1]*values_[1] + values_[2]*values_[2]);
        if(dist_ > 7) { return; }

        histo_.add(values_);
    }

    void
    finalize() {
        histo_.to_text_file("test.txt");
        histo_.to_binary_file(file_name_);
    }

protected:
    void
    _output_header(
            motif_data_structure::MotifStateTreeOP const & mst) {
        outputed_header_ = true;

    }


private:
    void
    _parse_constraints(
            String const & constraints) {

        auto spl =  base::split_str_by_delimiter(constraints, ";");
        for(auto const & s : spl) {
            auto spl2 = base::split_str_by_delimiter(s, ",");
            if(spl2.size() != 3) { throw SimulateTectosAppException("invalid record constraint: " + s); }
            auto pos = _parse_constraint_position(spl2[0]);
            auto lower = std::stod(spl2[1]);
            auto upper = std::stod(spl2[2]);
            constraints_[pos] = math::Real2{lower, upper};
        }

    }

    int
    _parse_constraint_position(
            String const & constraint_name) {
        if     (constraint_name == "x") { return 0; }
        else if(constraint_name == "y") { return 1; }
        else if(constraint_name == "z") { return 2; }
        else if(constraint_name == "a") { return 3; }
        else if(constraint_name == "b") { return 4; }
        else if(constraint_name == "g") { return 5; }
        else {
            throw SimulateTectosAppException("unknown constraint name: " + constraint_name);
        }

    }

    void
    _setup_constraints() {
        for(int i = 0; i < 3; i++) { constraints_[i] = math::Real2{ -100, 100}; }
        for(int i = 3; i < 6; i++) { constraints_[i] = math::Real2{ 0, 360}; }
    }

private:
    structure::BasepairOP ref_bp_;
    math::Matrix rot_, rot_t_, ref_r_t_, r_, r1_, r2_, r1_trans_, r2_trans_;
    math::Point d1_, d2_, d_;
    String file_name_;
    math::SixDHistogram histo_;
    math::Real6 values_;
    std::array<math::Real2, 6> constraints_;
    double dist_;



};

SimulateTectosApp::SimulateTectosApp() : base::Application(),
        tfs_(thermo_fluctuation::ThermoFlucSimulationDevel()),
        motif_names_(Strings()),
        ggaa_ttr_end_names_(Strings()),
        gaaa_ttr_end_names_(Strings()) {}


// application setups functions ////////////////////////////////////////////////////////////////////

void
SimulateTectosApp::setup_options() {
    add_option("fseq", "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG", base::OptionType::STRING);
    add_option("fss",  "((((((....((((((((((((....))))))))))))....))))))", base::OptionType::STRING);
    add_option("cseq", "CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG",  base::OptionType::STRING);
    add_option("css",  "(((((((..((((((((((((....))))))))))))...)))))))",  base::OptionType::STRING);
    add_option("s", 1000000, base::OptionType::INT);
    add_option("start_pose", false, base::OptionType::BOOL);
    add_option("start_pdbs", false, base::OptionType::BOOL);
    add_option("coorigin", false, base::OptionType::BOOL);
    add_option("gaaa_coorigin", false, base::OptionType::BOOL);

    //extra ensembles
    add_option("extra_me", "", base::OptionType::STRING);

    //new ggaa loop
    add_option("new_ggaa_model", false, base::OptionType::BOOL);
    add_option("ggaa_model", "", base::OptionType::STRING);
    
    add_option("extra_motifs", "", base::OptionType::STRING);

    // recording info from simulation
    add_option("record", false, base::OptionType::BOOL);
    add_option("record_constraints", "", base::OptionType::STRING);
    add_option("record_file", "test.out", base::OptionType::STRING);
    add_option("record_only_bound", false, base::OptionType::BOOL);
    add_option("record_only_unbound", false, base::OptionType::BOOL);
    add_option("record_file_type", "", base::OptionType::STRING);
    add_option("dump_state", false, base::OptionType::BOOL);
    add_option("dump_pdbs", false, base::OptionType::BOOL);

    add_option("scorer", "", base::OptionType::STRING);
    add_option("constraints", "", base::OptionType::STRING);


    //add_cl_options(tfs_.options(), "simulation");
    
}

void
SimulateTectosApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    
    base::Application::parse_command_line(argc, argv);
    
    //cl_parser_.assign_options(cl_options_, tfs_.options(), "simulation");
    tfs_.update_var_options();
}

void
SimulateTectosApp::run() {

    /*auto lines =base::get_lines_from_file("state.out");
    auto mt = std::make_shared<motif_data_structure::MotifTree>(lines[0], MotifTreeStringType::MT_STR);

    auto path = base::motif_dirs() + "ref.motif";
    auto ref_motif = motif::file_to_motif(path);
    auto logger = std::make_shared<SimulateTectosRecord6D>(
            "new.out",
            ref_motif->basepairs()[0]);
    auto mst = std::make_shared<motif_data_structure::MotifStateTree>(mt);
    logger->setup(mst, 1, 23, 1, 1);
    logger->log(mst, 0);

    exit(0);*/

    // load extra motifs in resource manager
    if(get_string_option("extra_motifs") != "") {
        auto spl = base::split_str_by_delimiter(get_string_option("extra_motifs"), ",");
        for(auto const & path : spl) {
            std::cout << path << std::endl;
            auto m = motif::file_to_motif(path);
            resources::Manager::instance().add_motif(m);
            
        }
    }
    
    auto fseq = get_string_option("fseq");
    auto fss  = get_string_option("fss");
    auto cseq = get_string_option("cseq");
    auto css  = get_string_option("css");
    auto mset = motif_data_structure::MotifStateEnsembleTreeOP(nullptr);
    
    if(get_string_option("extra_me") != "") {
        std::cout << "SIMULATE_TECTOS: registered extra motif ensembles from file: ";
        std::cout << get_string_option("extra_me") << std::endl;
        resources::Manager::instance().register_extra_motif_ensembles(get_string_option("extra_me"));
    }

    if     (get_bool_option("new_ggaa_model")) {
        mset = get_mset_new_receptor(remove_Us(fseq), fss, remove_Us(cseq), css);
    }
    else if(get_bool_option("coorigin")) {
        mset = get_mset_old_coorigin(remove_Us(fseq), fss, remove_Us(cseq), css);
    }
    else if(get_bool_option("gaaa_coorigin")) {
        mset = get_mset_old_gaaa_coorigin(remove_Us(fseq), fss, remove_Us(cseq), css);
    }

    else {
        mset = get_mset_old(remove_Us(fseq), fss, remove_Us(cseq), css);
    }

    if(get_bool_option("start_pose")) {
        auto mt = mset->to_mst()->to_motif_tree();
        mt->to_pdb("start_pose.pdb", 1, 1, 1);
        std::cout << "SIMULATE_TECTOS: outputing starting pose: start_pose.pdb" << std::endl;
    }
    
    if(get_bool_option("start_pdbs")) {
        auto mt = mset->to_mst()->to_motif_tree();
        mt->write_pdbs();
        std::cout << "SIMULATE_TECTOS: outputing each motif as nodes.*.pdb" << std::endl;

    }

    auto target_node_index = 1;
    auto target_node_end_index = -1;

    auto final_node_index = mset->last_node()->index();
    auto final_node_end_index = 1;

    auto steric_node_str = String("");
    steric_node_str  = std::to_string(final_node_index) + "," + std::to_string(final_node_index-1);

    if(get_bool_option("coorigin")) {
        for(auto const & n : *mset) {
            if(n->data()->get_member(0)->motif_state->name() == "GAAA_tetraloop") {
                target_node_index = n->index();
                break;
            }
        }
        if(target_node_index == -1) {
            throw SimulateTectosAppException("cannot find GAAA_tetraloop");
        }

        target_node_end_index = mset->get_node(target_node_index)->data()->get_member(0)->motif_state->get_end_index("A222-A251");

        steric_node_str += ":" + std::to_string(target_node_index);
    }
    else if(get_bool_option("gaaa_coorigin")) {
        for(auto const & n : *mset) {
            if(n->data()->get_member(0)->motif_state->name() == "GGAA_tetraloop") {
                target_node_index = n->index();
                break;
            }
        }
        if(target_node_index == -1) {
            throw SimulateTectosAppException("cannot find GGAA_tetraloop");
        }
        target_node_end_index = mset->get_node(target_node_index)->data()->get_member(0)->motif_state->get_end_index("A1-A6");

        steric_node_str += ":" + std::to_string(target_node_index);

    }

    else {
        steric_node_str += ":1";

        target_node_end_index = mset->get_node(1)->data()->get_member(0)->motif_state->get_end_index("A1-A6");
    }

    tfs_.set_option_value("steps", get_int_option("s"));
    tfs_.set_option_value("steric_nodes", steric_node_str);
    tfs_.setup(mset, target_node_index, final_node_index, target_node_end_index, final_node_end_index);

    tfs_.set_option_value("record_only_bound", get_bool_option("record_only_bound"));
    tfs_.set_option_value("record_only_unbound", get_bool_option("record_only_unbound"));
    tfs_.set_option_value("dump_state", get_bool_option("dump_state"));
    tfs_.set_option_value("dump_pdbs", get_bool_option("dump_pdbs"));

    if(get_bool_option("record")) {
        tfs_.set_option_value("record", true);
        auto logger = _get_logger(get_string_option("record_file_type"));
        tfs_.set_logger(logger);
    }

    auto scorer = _get_scorer(get_string_option("scorer"));
    tfs_.set_scorer(scorer);

    auto count = tfs_.run();
    std::cout << count << std::endl;
    
}


motif_data_structure::MotifStateEnsembleTreeOP
SimulateTectosApp::get_mset_old(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    auto ggaa_ttr = resources::Manager::instance().motif("GGAA_tetraloop", "", "A14-A15");
    auto gaaa_ttr = resources::Manager::instance().motif("GAAA_tetraloop", "", "A149-A154");

    for(auto const & end : gaaa_ttr->ends()) { gaaa_ttr_end_names_.push_back(end->name()); }
    for(auto const & end : ggaa_ttr->ends()) { ggaa_ttr_end_names_.push_back(end->name()); }
    
    auto flow_motifs = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motifs = get_motifs_from_seq_and_ss(cseq, css);
    
    auto mt = std::make_shared<motif_data_structure::MotifTree>();
    mt->set_option_value("sterics", false);
    auto m = resources::Manager::instance().bp_step("GG_LL_CC_RR");
    mt->add_motif(m);
    mt->add_motif(ggaa_ttr);
    motif_names_.push_back("ggaa_ttr");
    motif_names_.push_back("f1");
    mt->add_motif(flow_motifs[1], 1, "A7-A22");
    for(int i = 2; i < flow_motifs.size(); i++) {
        motif_names_.push_back("f"+std::to_string(i));
        mt->add_motif(flow_motifs[i]);
    }

    mt->add_motif(gaaa_ttr);
    mt->add_motif(chip_motifs[1], -1, "A222-A251");
    motif_names_.push_back("gaaa_ttr");
    motif_names_.push_back("c1");
    for(int i = 2; i < chip_motifs.size(); i++) {
        motif_names_.push_back("c"+std::to_string(i));
        mt->add_motif(chip_motifs[i]);
    }
    mt->add_connection(1, mt->last_node()->index(), "A1-A6", mt->last_node()->data()->end_name(1));

    //mt->write_pdbs();
    auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>(mt);
    return mset;
}

motif_data_structure::MotifStateEnsembleTreeOP
SimulateTectosApp::get_mset_new_receptor(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    if(get_string_option("ggaa_model") == "") {
        String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
        auto lines =base::get_lines_from_file(base_path+"new_ggaa_tetraloop.motif");
        std::cout << "SIMULATE_TECTOS: using custom ggaa_model: ";
        std::cout << base_path+"new_ggaa_tetraloop.motif" << std::endl;

        auto new_ggaa_tetraloop = std::make_shared<motif::Motif>(
                lines[0],
                structure::ResidueTypeSetManager::getInstance().residue_type_set());
        resources::Manager::instance().add_motif(new_ggaa_tetraloop);
    }
    else {
        auto lines =base::get_lines_from_file(get_string_option("ggaa_model"));
        std::cout << "SIMULATE_TECTOS: using custom ggaa_model: ";
        std::cout << get_string_option("ggaa_model") << std::endl;
        auto new_ggaa_tetraloop = std::make_shared<motif::Motif>(
            lines[0],
            structure::ResidueTypeSetManager::getInstance().residue_type_set());
        resources::Manager::instance().add_motif(new_ggaa_tetraloop);

    }
    
    auto ggaa_ttr = resources::Manager::instance().motif("new_ggaa_tetraloop", "", "A13-A16");
    auto gaaa_ttr = resources::Manager::instance().motif("GAAA_tetraloop", "", "A149-A154");

    for(auto const & end : gaaa_ttr->ends()) { gaaa_ttr_end_names_.push_back(end->name()); }
    for(auto const & end : ggaa_ttr->ends()) { ggaa_ttr_end_names_.push_back(end->name()); }

    auto flow_motifs = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motifs = get_motifs_from_seq_and_ss(cseq, css);
    
    auto mt = std::make_shared<motif_data_structure::MotifTree>();
    mt->set_option_value("sterics", false);
    auto m = resources::Manager::instance().bp_step("GG_LL_CC_RR");
    mt->add_motif(m);
    mt->add_motif(ggaa_ttr);
    mt->add_motif(flow_motifs[1], 1, "A7-A22");
    motif_names_.push_back("ggaa_ttr");
    motif_names_.push_back("f1");
    for(int i = 2; i < flow_motifs.size(); i++) {
        motif_names_.push_back("f"+std::to_string(i));
        mt->add_motif(flow_motifs[i]);
    }
    
    mt->add_motif(gaaa_ttr);
    mt->add_motif(chip_motifs[1], -1, "A222-A251");
    motif_names_.push_back("gaaa_ttr");
    motif_names_.push_back("c1");
    for(int i = 2; i < chip_motifs.size(); i++) {
        motif_names_.push_back("c"+std::to_string(i));
        mt->add_motif(chip_motifs[i]);
    }
    
    auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>(mt);
    return mset;
    
}

motif_data_structure::MotifStateEnsembleTreeOP
SimulateTectosApp::get_mset_old_coorigin(
        String const & fseq,
        String const & fss,
        String const & cseq,
        String const & css) {

    auto ggaa_ttr = resources::Manager::instance().motif("GGAA_tetraloop", "", "A14-A15");
    auto gaaa_ttr = resources::Manager::instance().motif("GAAA_tetraloop", "", "A149-A154");

    for(auto const & end : gaaa_ttr->ends()) { gaaa_ttr_end_names_.push_back(end->name()); }
    for(auto const & end : ggaa_ttr->ends()) { ggaa_ttr_end_names_.push_back(end->name()); }

    auto flow_motifs = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motifs = get_motifs_from_seq_and_ss(cseq, css);

    auto mt = std::make_shared<motif_data_structure::MotifTree>();
    mt->set_option_value("sterics", false);
    auto m = resources::Manager::instance().bp_step("GG_LL_CC_RR");
    mt->add_motif(m);
    auto pos = mt->add_motif(ggaa_ttr);
    motif_names_.push_back("ggaa_ttr");
    motif_names_.push_back("f1");
    mt->add_motif(flow_motifs[1], 1, "A7-A22");
    for(int i = 2; i < flow_motifs.size(); i++) {
        motif_names_.push_back("f"+std::to_string(i));
        mt->add_motif(flow_motifs[i]);
    }

    mt->add_motif(gaaa_ttr);
    auto new_chip_motifs = motif::MotifOPs();
    for(int i = chip_motifs.size()-1; i > -1; i--) {
        auto m = resources::Manager::instance().motif(chip_motifs[i]->name(), chip_motifs[i]->end_ids()[1], "");
        new_chip_motifs.push_back(m);
    }

    mt->add_motif(new_chip_motifs[0], pos, "A1-A6");
    int j = new_chip_motifs.size()-1;
    motif_names_.push_back("c"+std::to_string(j));
    j--;
    for(int i = 1; i < new_chip_motifs.size()-1; i++) {
        motif_names_.push_back("c"+std::to_string(j));
        mt->add_motif(new_chip_motifs[i]);
        j--;
    }

    auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>(mt);
    return mset;

}

motif_data_structure::MotifStateEnsembleTreeOP
SimulateTectosApp::get_mset_old_gaaa_coorigin(
        String const & fseq,
        String const & fss,
        String const & cseq,
        String const & css) {

    auto gaaa_ttr = resources::Manager::instance().motif("GAAA_tetraloop", "", "A229-A245");
    auto ggaa_ttr = resources::Manager::instance().motif("GGAA_tetraloop", "", "A7-A22");

    for(auto const & end : gaaa_ttr->ends()) { gaaa_ttr_end_names_.push_back(end->name()); }
    for(auto const & end : ggaa_ttr->ends()) { ggaa_ttr_end_names_.push_back(end->name()); }

    auto flow_motifs = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motifs = get_motifs_from_seq_and_ss(cseq, css);

    auto mt = std::make_shared<motif_data_structure::MotifTree>();
    mt->set_option_value("sterics", false);
    auto m = resources::Manager::instance().bp_step("GG_LL_CC_RR");
    mt->add_motif(m);
    auto pos = mt->add_motif(gaaa_ttr);
    motif_names_.push_back("gaaa_ttr");
    auto new_flow_motifs = motif::MotifOPs();
    for(int i = flow_motifs.size()-1; i > -1; i--) {
        auto m = resources::Manager::instance().motif(flow_motifs[i]->name(), flow_motifs[i]->end_ids()[1], "");
        new_flow_motifs.push_back(m);
    }
    mt->add_motif(new_flow_motifs[0], pos, "A149-A154");
    int j = new_flow_motifs.size()-1;
    motif_names_.push_back("f"+std::to_string(j));
    j--;
    for(int i = 1; i < new_flow_motifs.size()-1; i++) {
        motif_names_.push_back("f"+std::to_string(j));
        mt->add_motif(new_flow_motifs[i]);
        j--;
    }
    mt->add_motif(ggaa_ttr);
    mt->add_motif(chip_motifs[1], pos, "A222-A251");
    motif_names_.push_back("ggaa_ttr");
    motif_names_.push_back("c1");
    for(int i = 2; i < chip_motifs.size(); i++) {
        motif_names_.push_back("c"+std::to_string(i));
        mt->add_motif(chip_motifs[i]);
    }

    auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>(mt);
    return mset;
}

motif::MotifOPs
SimulateTectosApp::get_motifs_from_seq_and_ss(
        String const & seq,
        String const & ss) {
    
    auto parser = secondary_structure::Parser();
    auto ss_motifs = parser.parse_to_motifs(seq, ss);
    auto motifs = motif::MotifOPs();
    
    auto start = 0;
    auto motif = motif::MotifOP(nullptr);
    for(auto const & m : ss_motifs) {
        if(m->mtype() == util::MotifType::TWOWAY && start == 0) {
            start = 1;
            continue;
        }
        
        if(m->mtype() == util::MotifType::HAIRPIN) { break; }
        if(!start) { continue; }
                
        //basepair step
        if(m->mtype() == util::MotifType::HELIX) {
            motif = resources::Manager::instance().bp_step(m->end_ids()[0]);
            motifs.push_back(motif);
        }
        else if(m->mtype() == util::MotifType::TWOWAY) {
            auto end_id = m->end_ids()[0];
            try {
                motif = resources::Manager::instance().motif("", end_id);
            }
            catch(resources::ResourceManagerException const & e) {
                throw SimulateTectosAppException(
                    "cannot find a motif that corresponds to the sequence: " + m->sequence() +
                    " and secondary structure: " + m->dot_bracket() + "for the simulation");
            }
            motifs.push_back(motif);
            
        }
        else {
            throw SimulateTectosAppException(
                    "motif_type: " + type_to_str(m->mtype()) + " is not supported in tecto "
                    "simulations currently only TWOWAY/HELIX are supported");
        }
    }
    
    return motifs;
}

thermo_fluctuation::ThermoFlucSimulationLoggerOP
SimulateTectosApp::_get_logger(
        String const & name) {
    if(name.length() == 0) {
        return std::make_shared<SimulateTectosRecordAllLogger>(
                get_string_option("record_file"),
                motif_names_,
                ggaa_ttr_end_names_,
                gaaa_ttr_end_names_);
    }
    else if (name == "AllLoger") {
        return std::make_shared<SimulateTectosRecordAllLogger>(
                get_string_option("record_file"),
                motif_names_,
                ggaa_ttr_end_names_,
                gaaa_ttr_end_names_);
    }

    else if(name == "Record6D") {
        auto path = base::motif_dirs() + "ref.motif";
        auto ref_motif = motif::file_to_motif(path);
        return std::make_shared<SimulateTectosRecord6D>(
                get_string_option("record_file"),
                get_string_option("record_constraints"),
                ref_motif->basepairs()[0]);
    }
    else if(name == "Record6DHisto" || name == "Record6DHistogram") {
        auto path = base::motif_dirs() + "ref.motif";
        auto ref_motif = motif::file_to_motif(path);
        return std::make_shared<SimulateTectosRecord6DHistogram>(
                get_string_option("record_file"),
                get_string_option("record_constraints"),
                ref_motif->basepairs()[0]);
    }

    else {
        throw SimulateTectosAppException("unknown logger type: " + name);
    }
}

thermo_fluctuation::ThermoFlucScorerOP
SimulateTectosApp::_get_scorer(
        String const & name) {

    if(name.length() == 0) {
        return std::make_shared<thermo_fluctuation::FrameScorerDevel>();
    }

    else if(name == "FrameScorer") {
        return std::make_shared<thermo_fluctuation::FrameScorerDevel>();
    }

    else if(name == "SixDScorer") {
        auto path = base::motif_dirs() + "ref.motif";
        auto ref_motif = motif::file_to_motif(path);
        return std::make_shared<thermo_fluctuation::SixDScorer>(
                get_string_option("constraints"),
                ref_motif->basepairs()[0]);
    }

    else {
        throw SimulateTectosAppException("unknown scorer: " + name);
    }
}



// non-member functions ////////////////////////////////////////////////////////////////////////////


String
remove_Us(String const & seq) {
    auto seq_rna = String("");
    
    //hacky way to convert Ts to Us
    for(auto const & e : seq) {
        if(e == 'T' ) { seq_rna += 'U'; }
        else          { seq_rna += e;   }
    }
    
    return seq_rna;
    
}


// main ////////////////////////////////////////////////////////////////////////////////////////////




