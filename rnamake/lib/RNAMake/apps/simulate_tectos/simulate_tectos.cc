//
//  simulate_tectos.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/cl_option.h"
#include "util/settings.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_tree.h"
#include "thermo_fluctuation/thermo_fluc_simulation.h"
#include "secondary_structure/secondary_structure_tree.h"
#include "secondary_structure/secondary_structure_factory.h"
#include "simulate_tectos.h"

Options
parse_command_line(
    int argc,
    const char ** argv) {
    
    CL_Options cl_opts;
    cl_opts.add_option("fseq", "", STRING_TYPE,
                       "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG", false);
    cl_opts.add_option("fss" , "", STRING_TYPE,
                       "((((((....((((((((((((....))))))))))))....))))))", false);
    cl_opts.add_option("cseq", "", STRING_TYPE,
                       "CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG", false);
    cl_opts.add_option("css" , "", STRING_TYPE,
                       "(((((((..((((((((((((....))))))))))))...)))))))", false);
    cl_opts.add_option("s", "steps", FLOAT_TYPE, "1000000", false);
    
    return cl_opts.parse_command_line(argc, argv);
    
}

SimulateTectos::SimulateTectos(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    auto mset = get_mset_old(fseq, fss, cseq, css);
    exit(0);
    ThermoFlucSimulation tfs;
    tfs.setup(mset, 1, mset->last_node()->index(), 1, 1);
    tfs.option("steps", 1000000);
    tfs.option("cutoff", 4.5f);
    //tfs.option("cutoff", 5.0f);

    int count = tfs.run();
    std::cout << count << std::endl;
}



MotifStateEnsembleTreeOP
SimulateTectos::get_mset_old(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
   
    auto flow_motif_names = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motif_names = get_motifs_from_seq_and_ss(cseq, css);
    auto mt = std::make_shared<MotifTree>();
    mt->option("sterics", 0);
    mt->add_motif("GC=GC");
    mt->add_motif("GGAA_tetraloop", "A14-A15");
    mt->add_motif(flow_motif_names[1], -1, "A7-A22");
    /*for(int i = 2; i < flow_motif_names.size(); i++) {
        mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_names[i]));
    }
    mt->add_motif(ResourceManager::getInstance().get_motif("GAAA_tetraloop", "", "A149-A154"));
    mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_names[0]),
                  -1, -1, "A222-A251");
    
    for(int i = 2; i < chip_motif_names.size(); i++) {
        mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_names[i]));
    }
    */
    MotifStateEnsembleTreeOP mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);

    return mset;
}

Strings
SimulateTectos::get_motifs_from_seq_and_ss(
    String const & seq,
    String const & ss) {
    
    auto ssf = sstruct::SecondaryStructureFactory();
    auto p = ssf.pose(seq, ss);
    auto sst = sstruct::tree_from_pose(p);
    auto motif_names = Strings();
    String motif_name, motif_name_rna;

    int start = 0;
    for(auto const & n : *sst) {
        if(n->data()->mtype() == MotifType::TWOWAY) {
            start = 1;
            continue;
        }
        if(n->data()->mtype() == MotifType::HAIRPIN) { break; }
        
        if(!start) { continue; }
        
        if(n->data()->mtype() != MotifType::HELIX) {
            throw std::runtime_error("old method does not have non helical motifs implemented yet!!!!!");
        }
        
        motif_name = bp_name_from_sequence(n->data()->sequence());
        motif_name_rna = "";
        for(auto const & e : motif_name) {
            if(e == 'T' ) { motif_name_rna += 'U'; }
            else          { motif_name_rna += e;   }
        }
        motif_names.push_back(motif_name_rna);

        
    }

    return motif_names;
}

String
bp_name_from_sequence(
    String const & seq) {
    
    String name = "";
    name += seq[0]; name += seq[4];
    name += "=";
    name += seq[1]; name += seq[3];
    return name;
    
}



int main(int argc, const char * argv[]) {
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    ResourceManager::getInstance().add_motif(base_path+"GAAA_tetraloop");
    ResourceManager::getInstance().add_motif(base_path+"GGAA_tetraloop");

    try {
        
        Options opts = parse_command_line(argc, argv);
        SimulateTectos st(opts.option<String>("fseq"),
                          opts.option<String>("fss"),
                          opts.option<String>("cseq"),
                          opts.option<String>("css"));
                                              
    } catch(std::runtime_error e) {
        std::cerr << "caught runtime exception: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
}















