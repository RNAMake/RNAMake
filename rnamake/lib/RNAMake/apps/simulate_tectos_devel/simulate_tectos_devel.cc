//
//  simulate_tectos.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include "base/cl_option.h"
#include "base/settings.h"
#include "structure/residue_type_set_manager.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_tree.h"
#include "simulate_tectos_devel.h"


SimulateTectosApp::SimulateTectosApp() : Application(),
tfs_(ThermoFlucSimulationDevel())
{}


// application setups functions ////////////////////////////////////////////////////////////////////

void
SimulateTectosApp::setup_options() {
    add_option("fseq", "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG", OptionType::STRING);
    add_option("fss",  "((((((....((((((((((((....))))))))))))....))))))", OptionType::STRING);
    add_option("cseq", "CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG",  OptionType::STRING);
    add_option("css",  "(((((((..((((((((((((....))))))))))))...)))))))",  OptionType::STRING);
    add_option("s", 1000000, OptionType::INT);
    add_option("start_pose", false, OptionType::BOOL);
    add_option("start_pdbs", false, OptionType::BOOL);

    
    //extra ensembles
    add_option("extra_me", "", OptionType::STRING);
    
    //record options
    add_option("r", false, OptionType::BOOL);
    add_option("record_file", "test.out", OptionType::STRING);

    //new ggaa loop
    add_option("new_ggaa_model", false, OptionType::BOOL);
    add_option("ggaa_model", "", OptionType::STRING);
    
    add_option("extra_motifs", "", OptionType::STRING);
    
    add_cl_options(tfs_.options(), "simulation");
    
}

void
SimulateTectosApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    
    Application::parse_command_line(argc, argv);
    
    cl_parser_.assign_options(cl_options_, tfs_.options(), "simulation");
    tfs_.update_var_options();
}

void
SimulateTectosApp::run() {
    
    // load extra motifs in resource manager
    if(get_string_option("extra_motifs") != "") {
        auto spl = split_str_by_delimiter(get_string_option("extra_motifs"), ",");
        for(auto const & path : spl) {
            std::cout << path << std::endl;
            auto m = file_to_motif(path);
            RM::instance().add_motif(m);
            
        }
    }
    
    auto fseq = get_string_option("fseq");
    auto fss  = get_string_option("fss");
    auto cseq = get_string_option("cseq");
    auto css  = get_string_option("css");
    auto mset = MotifStateEnsembleTreeOP(nullptr);
    
    if(get_string_option("extra_me") != "") {
        std::cout << "SIMULATE_TECTOS: registered extra motif ensembles from file: ";
        std::cout << get_string_option("extra_me") << std::endl;
        RM::instance().register_extra_motif_ensembles(get_string_option("extra_me"));
    }
    
    if(get_bool_option("new_ggaa_model")) {
        mset = get_mset_new_receptor(remove_Us(fseq), fss, remove_Us(cseq), css);
    }
    else {
        mset = get_mset_old(remove_Us(fseq), fss, remove_Us(cseq), css);
    }
    
    
    if(get_bool_option("start_pose")) {
        auto mt = mset->to_mst()->to_motif_tree();
        mt->to_pdb("start_pose.pdb");
        std::cout << "SIMULATE_TECTOS: outputing starting pose: start_pose.pdb" << std::endl;
    }
    
    if(get_bool_option("start_pdbs")) {
        auto mt = mset->to_mst()->to_motif_tree();
        mt->write_pdbs();
        std::cout << "SIMULATE_TECTOS: outputing each motif as nodes.*.pdb" << std::endl;

    }
    
    auto steric_node_str = String("");
    auto last_node_index = mset->last_node()->index();
    steric_node_str += std::to_string(last_node_index) + "," + std::to_string(last_node_index-1);
    steric_node_str += ":1";
    
    auto end_index = mset->to_mst()->get_node(1)->data()->get_end_index("A1-A6");
    
    tfs_.set_option_value("steps", get_int_option("s"));
    tfs_.set_option_value("steric_nodes", steric_node_str);
    tfs_.set_option_value("record", get_bool_option("r"));
    tfs_.setup(mset, 1, mset->last_node()->index(), end_index, 1);
    auto count = tfs_.run();
    std::cout << count << std::endl;
    
}


MotifStateEnsembleTreeOP
SimulateTectosApp::get_mset_old(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    auto ggaa_ttr = RM::instance().motif("GGAA_tetraloop", "", "A14-A15");
    auto gaaa_ttr = RM::instance().motif("GAAA_tetraloop", "", "A149-A154");
    
    auto flow_motifs = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motifs = get_motifs_from_seq_and_ss(cseq, css);
    
    auto mt = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", false);
    auto m = RM::instance().bp_step("GG_LL_CC_RR");
    mt->add_motif(m);
    mt->add_motif(ggaa_ttr);
    mt->add_motif(flow_motifs[1], 1, "A7-A22");
    for(int i = 2; i < flow_motifs.size(); i++) {
        mt->add_motif(flow_motifs[i]);
    }
    
    mt->add_motif(gaaa_ttr);
    mt->add_motif(chip_motifs[1], -1, "A222-A251");
    for(int i = 2; i < chip_motifs.size(); i++) {
        mt->add_motif(chip_motifs[i]);
    }
    
    //mt->write_pdbs();
    auto mset = std::make_shared<MotifStateEnsembleTree>(mt);
    return mset;
}

MotifStateEnsembleTreeOP
SimulateTectosApp::get_mset_new_receptor(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    if(get_string_option("ggaa_model") == "") {
        String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
        auto lines = get_lines_from_file(base_path+"new_ggaa_tetraloop.motif");
        std::cout << "SIMULATE_TECTOS: using custom ggaa_model: ";
        std::cout << base_path+"new_ggaa_tetraloop.motif" << std::endl;

        auto new_ggaa_tetraloop = std::make_shared<Motif>(
                lines[0],
                ResidueTypeSetManager::getInstance().residue_type_set());
        RM::instance().add_motif(new_ggaa_tetraloop);
    }
    else {
        auto lines = get_lines_from_file(get_string_option("ggaa_model"));
        std::cout << "SIMULATE_TECTOS: using custom ggaa_model: ";
        std::cout << get_string_option("ggaa_model") << std::endl;
        auto new_ggaa_tetraloop = std::make_shared<Motif>(
            lines[0],
            ResidueTypeSetManager::getInstance().residue_type_set());
        RM::instance().add_motif(new_ggaa_tetraloop);

    }
    
    auto ggaa_ttr = RM::instance().motif("new_ggaa_tetraloop", "", "A13-A16");
    auto gaaa_ttr = RM::instance().motif("GAAA_tetraloop", "", "A149-A154");
    
    auto flow_motifs = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motifs = get_motifs_from_seq_and_ss(cseq, css);
    
    auto mt = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", false);
    auto m = RM::instance().bp_step("GG_LL_CC_RR");
    mt->add_motif(m);
    mt->add_motif(ggaa_ttr);
    mt->add_motif(flow_motifs[1], 1, "A7-A22");
    for(int i = 2; i < flow_motifs.size(); i++) {
        mt->add_motif(flow_motifs[i]);
    }
    
    mt->add_motif(gaaa_ttr);
    mt->add_motif(chip_motifs[1], -1, "A222-A251");
    for(int i = 2; i < chip_motifs.size(); i++) {
        mt->add_motif(chip_motifs[i]);
    }
    
    auto mset = std::make_shared<MotifStateEnsembleTree>(mt);
    return mset;
    
}

MotifOPs
SimulateTectosApp::get_motifs_from_seq_and_ss(
    String const & seq,
    String const & ss) {
    
    auto parser = sstruct::SecondaryStructureParser();
    auto ss_motifs = parser.parse_to_motifs(seq, ss);
    auto motifs = MotifOPs();
    
    auto start = 0;
    auto motif = MotifOP(nullptr);
    for(auto const & m : ss_motifs) {
        if(m->mtype() == MotifType::TWOWAY && start == 0) {
            start = 1;
            continue;
        }
        
        if(m->mtype() == MotifType::HAIRPIN) { break; }
        if(!start) { continue; }
                
        //basepair step
        if(m->mtype() == MotifType::HELIX) {
            motif = RM::instance().bp_step(m->end_ids()[0]);
            motifs.push_back(motif);
        }
        else if(m->mtype() == MotifType::TWOWAY) {
            auto end_id = m->end_ids()[0];
            try {
                motif = RM::instance().motif("", end_id);
            }
            catch(ResourceManagerException const & e) {
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




