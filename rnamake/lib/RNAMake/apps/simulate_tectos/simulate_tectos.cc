//
//  simulate_tectos.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/cl_option.h"
#include "base/settings.h"
#include "base/backtrace.hpp"
#include "secondary_structure/secondary_structure_parser.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_tree.h"
#include "simulate_tectos.h"


SimulateTectosApp::SimulateTectosApp() : Application(),
tfs_(ThermoFlucSimulation())
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


// run functions ///////////////////////////////////////////////////////////////////////////////////

void
SimulateTectosApp::run() {
    
    auto fseq = get_string_option("fseq");
    auto fss  = get_string_option("fss");
    auto cseq = get_string_option("cseq");
    auto css  = get_string_option("css");
    
    auto mset = get_mset_old(fseq, fss, cseq, css);
    
    if(get_bool_option("start_pose")) {
        auto mt = mset->to_mst()->to_motif_tree();
        std::cout << "outputing starting pose: start_pose.pdb" << std::endl;
        mt->to_pdb("start_pose.pdb", 1);
    }
    
    
    auto steric_node_str = String("");
    auto last_node_index = mset->last_node()->index();
    steric_node_str += std::to_string(last_node_index) + "," + std::to_string(last_node_index-1);
    steric_node_str += ":1";
    
    tfs_.set_option_value("steps", get_int_option("s"));
    tfs_.set_option_value("steric_nodes", steric_node_str);
    tfs_.setup(mset, 1, mset->last_node()->index(), 1, 1);
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
    auto m = RM::instance().motif("GC=GC");
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
            auto name = bp_name_from_sequence(m->sequence());
            motif = RM::instance().motif(name);
            motifs.push_back(motif);
        }
        else if(m->mtype() == MotifType::TWOWAY) {
            auto end_id = m->end_ids()[0];
            try {
                motif = RM::instance().motif("", end_id);
            }
            catch(ResourceManagerException const & e) {
                throw SimulateTectosAppException(
                    "cannot find a motif that corresponds to the sequence: " + motif->sequence() +
                    " and secondary structure: " + motif->dot_bracket() + "for the simulation");
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
bp_name_from_sequence(
    String const & seq) {
    
    auto name = String("");
    auto name_rna = String("");
    name += seq[0]; name += seq[4];
    name += "=";
    name += seq[1]; name += seq[3];
    
    //hacky way to convert Ts to Us
    for(auto const & e : name) {
        if(e == 'T' ) { name_rna += 'U'; }
        else          { name_rna += e;   }
    }

    return name_rna;
    
}


// main ////////////////////////////////////////////////////////////////////////////////////////////



int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);
    
    //load extra motifs being used
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    RM::instance().add_motif(base_path+"GAAA_tetraloop");
    RM::instance().add_motif(base_path+"GGAA_tetraloop");
    
    auto app = SimulateTectosApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;
}















