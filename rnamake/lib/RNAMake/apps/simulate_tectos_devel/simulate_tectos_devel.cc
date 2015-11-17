//
//  simulate_tectos.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/cl_option.h"
#include "util/settings.h"
#include "secondary_structure/ss_tree.h"
#include "resources/resource_manager.h"
#include "motif/motif_tree.h"
#include "thermo_fluctuation/thermo_fluc_simulation_devel.h"
#include "simulate_tectos_devel.h"

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
    cl_opts.add_option("s", "steps", INT_TYPE, "1000000", false);
    cl_opts.add_option("c", "cutoff", FLOAT_TYPE, "4.5f", false);
    cl_opts.add_option("t", "temperature", FLOAT_TYPE, "298.15f", false);
    cl_opts.add_option("sr", "steric_radius", FLOAT_TYPE, "2.2f", false);
    cl_opts.add_option("wd", "", FLOAT_TYPE, "1.0f", false);
    cl_opts.add_option("wr", "", FLOAT_TYPE, "1.0f", false);


    cl_opts.add_option("", "static", INT_TYPE, "0", false);
    cl_opts.add_option("r", "record", INT_TYPE, "0", false);
    cl_opts.add_option("rf", "record_file", STRING_TYPE, "test.out", false);
    cl_opts.add_option("rs", "recore_state", INT_TYPE, "0", false);
    cl_opts.add_option("rall", "record_all", INT_TYPE, "0", false);
    cl_opts.add_option("pdbs", "", INT_TYPE, "0", false);
    cl_opts.add_option("ensembles", "", INT_TYPE, "0", false);
    
    return cl_opts.parse_command_line(argc, argv);
    
}

SimulateTectos::SimulateTectos(
    Options & opts) {
    
    auto mset = get_mset_old(opts.option<String>("fseq"),
                             opts.option<String>("fss"),
                             opts.option<String>("cseq"),
                             opts.option<String>("css"));
    ThermoFlucSimulationDevel tfs;
    tfs.setup(mset, 1, mset->last_node()->index(), 1, 1);
    tfs.option("steps",  opts.option<int>("s"));
    tfs.option("cutoff", opts.option<float>("c"));
    tfs.option("temperature", opts.option<float>("t"));
    tfs.option("steric_radius", opts.option<float>("sr"));
    tfs.option("record", opts.option<int>("r"));
    tfs.option("record_file", opts.option<String>("rf"));
    tfs.option("record_state", opts.option<int>("rs"));
    tfs.option("record_all", opts.option<int>("rall"));
    tfs.option("d_weight", opts.option<float>("wd"));
    tfs.option("r_weight", opts.option<float>("wr"));

    if(opts.option<int>("ensembles")) {
        auto lib = MotifStateEnsembleSqliteLibrary("bp_steps");
        lib.load_all();
        for(auto const & mes : lib) {
            std::ofstream out;
            out.open(mes->id() + ".ensemble");
            for (auto const & mem : mes->members()) {
                out << mem->motif_state->name() << " " << mem->energy << std::endl;
            }
            out.close();
        }
        exit(0);
    }
    
    if (opts.option<int>("pdbs")) {
        mset->to_mst()->write_pdbs();
        exit(0);
    }
    
    if(opts.option<int>("static")) {
        std::cout << tfs.static_run() << std::endl;
        exit(0);
    }
    
    int count = tfs.run();
    std::cout << count << std::endl;
}



MotifStateEnsembleTreeOP
SimulateTectos::get_mset_old(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    auto mt = std::make_shared<MotifTree>();
    mt->option("sterics", 0);
    mt->add_motif(ResourceManager::getInstance().get_motif("GC=GC"));
    mt->add_motif(ResourceManager::getInstance().get_motif("GGAA_tetraloop", "", "A14-A15"));
    auto flow_motif_names = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motif_names = get_motifs_from_seq_and_ss(cseq, css);
    mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_names[1]),
                  -1, -1, "A7-A22");
    for(int i = 2; i < flow_motif_names.size(); i++) {
        mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_names[i]));
    }
    mt->add_motif(ResourceManager::getInstance().get_motif("GAAA_tetraloop", "", "A149-A154"));
    mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_names[1]),
                  -1, -1, "A222-A251");
    
    for(int i = 2; i < chip_motif_names.size(); i++) {
        mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_names[i]));
    }
    
    MotifStateEnsembleTreeOP mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);

    return mset;
}

Strings
SimulateTectos::get_motifs_from_seq_and_ss(
    String const & seq,
    String const & ss) {
    
    sstruct::SS_Tree ss_tree(seq, ss);
    sstruct::SS_TreeNodeOP current = nullptr;
    sstruct::SS_TreeNodeOPs required_nodes;
    for(auto const & n : ss_tree) {
        if(n->data()->type() == sstruct::SS_NodeData::SS_Type::SS_BULGE) {
            current = n;
            break;
        }
    }
    
    if(current == nullptr) {
        throw std::runtime_error("cannot find start node in get_motifs_from_seq_and_ss");
    }
    
    current = current->children()[0];
    while(current->data()->type() != sstruct::SS_NodeData::SS_Type::SS_HAIRPIN) {
        required_nodes.push_back(current);
        current = current->children()[0];
    }

    Strings motif_names;
    String motif_name, motif_name_rna;
    String seq1, seq2;
    for(int i = 1; i < required_nodes.size(); i++) {
        if(required_nodes[i]->data()->type()   != sstruct::SS_NodeData::SS_Type::SS_BP ||
           required_nodes[i-1]->data()->type() != sstruct::SS_NodeData::SS_Type::SS_BP) {
            throw std::runtime_error("old method does not have non helical motifs implemented yet!!!!!");
        }
        
        seq1 = required_nodes[i-1]->data()->sequence();
        seq2 = required_nodes[i]->data()->sequence();
        motif_name = "";
        motif_name.push_back(seq1[0]); motif_name.push_back(seq1[2]);
        motif_name.push_back('=');
        motif_name.push_back(seq2[0]); motif_name.push_back(seq2[2]);
        motif_name_rna = "";
        for(auto const & e : motif_name) {
            if(e == 'T' ) { motif_name_rna += 'U'; }
            else          { motif_name_rna += e;   }
        }
        motif_names.push_back(motif_name_rna);
        
    }
    
    return motif_names;
}



int main(int argc, const char * argv[]) {
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos_devel/resources/";
    ResourceManager::getInstance().add_motif(base_path+"GAAA_tetraloop");
    ResourceManager::getInstance().add_motif(base_path+"GGAA_tetraloop");

    try {
        
        Options opts = parse_command_line(argc, argv);
        SimulateTectos st(opts);
                                              
    } catch(std::runtime_error e) {
        std::cerr << "caught runtime exception: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
}