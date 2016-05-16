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
    cl_opts.add_option("rallf", "record_all_file", STRING_TYPE, "test_all.out", false);

    cl_opts.add_option("pdbs", "", INT_TYPE, "0", false);
    cl_opts.add_option("bound_pdb", "", INT_TYPE, "0", false);
    cl_opts.add_option("ensembles", "", INT_TYPE, "0", false);
    cl_opts.add_option("extra_mse", "", STRING_TYPE, "", false);
    cl_opts.add_option("full_seq", "", INT_TYPE, "0", false);
    cl_opts.add_option("reverse", "", INT_TYPE, "0", false);
    cl_opts.add_option("coorigin", "", INT_TYPE, "0", false);

    return cl_opts.parse_command_line(argc, argv);
    
}

SimulateTectos::SimulateTectos(
    Options & opts) {
    
    if(opts.option<String>("extra_mse").length() > 0) {
        ResourceManager::getInstance().register_extra_motif_state_ensembles(
                                            opts.option<String>("extra_mse"));
    }
    
    auto mset = std::make_shared<MotifStateEnsembleTree>();
    ThermoFlucSimulationDevel tfs;
    
    if(opts.option<int>("full_seq")) {
        mset = get_mset_old_full_seq(
                            opts.option<String>("fseq"),
                            opts.option<String>("fss"),
                            opts.option<String>("cseq"),
                            opts.option<String>("css"));
        
        tfs.setup(mset, 4, 25, 1, 1);
        tfs.check_nodes({25, 24});
        tfs.check_nodes_2({4});
    }
    
    else if(opts.option<int>("reverse")) {
        mset = get_mset_old_reverse(
                            opts.option<String>("fseq"),
                            opts.option<String>("fss"),
                            opts.option<String>("cseq"),
                            opts.option<String>("css"));
        tfs.setup(mset, 1, mset->last_node()->index(), 1, 1);
        //Ints nodes = {0};
        //tfs.check_nodes_2(nodes);

    }
    
    else if(opts.option<int>("coorigin")) {
        mset = get_mset_old_coorigin(
                            opts.option<String>("fseq"),
                            opts.option<String>("fss"),
                            opts.option<String>("cseq"),
                            opts.option<String>("css"));
        auto mst = mset->to_mst();
   
        tfs.check_nodes_2({12});
        tfs.setup(mset, 12, mset->last_node()->index(), 1, 1);

    }
    
    else {
        mset = get_mset_old(opts.option<String>("fseq"),
                            opts.option<String>("fss"),
                            opts.option<String>("cseq"),
                            opts.option<String>("css"));
        tfs.setup(mset, 1, mset->last_node()->index(), 1, 1);
        
    }
    

    tfs.option("steps",  opts.option<int>("s"));
    tfs.option("cutoff", opts.option<float>("c"));
    tfs.option("temperature", opts.option<float>("t"));
    tfs.option("steric_radius", opts.option<float>("sr"));
    tfs.option("record", opts.option<int>("r"));
    tfs.option("record_file", opts.option<String>("rf"));
    tfs.option("record_state", opts.option<int>("rs"));
    tfs.option("record_all", opts.option<int>("rall"));
    tfs.option("record_all_file", opts.option<String>("rallf"));
    tfs.option("d_weight", opts.option<float>("wd"));
    tfs.option("r_weight", opts.option<float>("wr"));
    tfs.option("bound_pdb", opts.option<int>("bound_pdb"));

    if(opts.option<int>("ensembles")) {
        auto lib = MotifStateEnsembleSqliteLibrary("bp_steps");
        lib.load_all();
        for(auto const & mes : lib) {
            std::ofstream out;
            out.open(mes->id() + ".ensemble");
            out << "name,energy,d,r" << std::endl;
            for (auto const & mem : mes->members()) {
                out << mem->motif_state->name() << "," << mem->energy << "," << vector_to_str(mem->motif_state->end_states()[1]->d()) << "," << matrix_to_str(mem->motif_state->end_states()[1]->r()) << std::endl;
            }
            out.close();
        }
        exit(0);
    }
    
    if (opts.option<int>("pdbs")) {
        auto mst = mset->to_mst();
        mst->write_pdbs();
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
    auto flow_motif_infos = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motif_infos = get_motifs_from_seq_and_ss(cseq, css);
    mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[1].name),
                  -1, -1, "A7-A22");
    
    for(int i = 2; i < flow_motif_infos.size(); i++) {
        if(flow_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name,
                                                                   flow_motif_infos[i].end_id,
                                                                   flow_motif_infos[i].end_name));
   
        }
    }
    
    int pos = mt->add_motif(ResourceManager::getInstance().get_motif("GAAA_tetraloop", "", "A149-A154"));
    mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[1].name),
                  -1, -1, "A222-A251");
    
    for(int i = 2; i < chip_motif_infos.size(); i++) {
        if(chip_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[i].name,
                                                                   chip_motif_infos[i].end_id,
                                                                   chip_motif_infos[i].end_name));
            
        }
    }
    
    MotifStateEnsembleTreeOP mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);

    return mset;
}


MotifStateEnsembleTreeOP
SimulateTectos::get_mset_old_reverse(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    auto flow_motif_infos = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motif_infos = get_motifs_from_seq_and_ss(cseq, css);
    auto mt = std::make_shared<MotifTree>();
    mt->option("sterics", 0);
    
    mt->add_motif(ResourceManager::getInstance().get_motif("GC=GC"));
    mt->add_motif(ResourceManager::getInstance().get_motif("GAAA_tetraloop", "", "A229-A245"));
    
    mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[1].name),
                  -1, -1, "A222-A251");
    
    for(int i = 2; i < chip_motif_infos.size(); i++) {
        if(chip_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[i].name,
                                                                   chip_motif_infos[i].end_id,
                                                                   chip_motif_infos[i].end_name));
            
        }
    }

    mt->add_motif(ResourceManager::getInstance().get_motif("GGAA_tetraloop", "", "A1-A6"));

    
    mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[1].name),
                  -1, -1, "A7-A22");
    
    for(int i = 2; i < flow_motif_infos.size(); i++) {
        if(flow_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name,
                                                                   flow_motif_infos[i].end_id,
                                                                   flow_motif_infos[i].end_name));
            
        }
    }
    
    MotifStateEnsembleTreeOP mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);
    
    return mset;
}

MotifStateEnsembleTreeOP
SimulateTectos::get_mset_old_full_seq(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    auto mt = std::make_shared<MotifTree>();
    mt->option("sterics", 0);
    auto steps = Strings{"CG=UA", "UA=AU", "AU=GC", "GC=GC"};
    //mt->add_motif(ResourceManager::getInstance().get_motif("GC=GC"));
    for(auto const & s : steps) {
        mt->add_motif(ResourceManager::getInstance().get_motif(s));
    }
    
    mt->add_motif(ResourceManager::getInstance().get_motif("GGAA_tetraloop", "", "A14-A15"));
    auto flow_motif_infos = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motif_infos = get_motifs_from_seq_and_ss(cseq, css);
    mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[1].name),
                  -1, -1, "A7-A22");
    
    for(int i = 2; i < flow_motif_infos.size(); i++) {
        if(flow_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name,
                                                                   flow_motif_infos[i].end_id,
                                                                   flow_motif_infos[i].end_name));
            
        }
    }
    
    int pos = mt->add_motif(ResourceManager::getInstance().get_motif("GAAA_tetraloop", "", "A149-A154"));
    mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[1].name),
                  -1, -1, "A222-A251");
    
    for(int i = 2; i < chip_motif_infos.size(); i++) {
        if(chip_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(chip_motif_infos[i].name,
                                                                   chip_motif_infos[i].end_id,
                                                                   chip_motif_infos[i].end_name));
            
        }
    }
    
    steps = Strings{"CG=CG", "CG=UA", "UA=AU", "AU=GC"};
    int j = 0;
    //mt->add_motif(ResourceManager::getInstance().get_motif("GC=GC"));
    for(auto const & s : steps) {
        if(j == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(s), pos, -1, "A229-A245");
        }
        else {
            mt->add_motif(ResourceManager::getInstance().get_motif(s));
        }
        j++;
    }
    
    MotifStateEnsembleTreeOP mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);
    
    return mset;
}

MotifStateEnsembleTreeOP
SimulateTectos::get_mset_old_coorigin(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    
    auto mt = std::make_shared<MotifTree>();
    mt->option("sterics", 0);
    mt->add_motif(ResourceManager::getInstance().get_motif("GC=GC"));
    int pos = mt->add_motif(ResourceManager::getInstance().get_motif("GGAA_tetraloop", "", "A14-A15"));
    
    auto flow_motif_infos = get_motifs_from_seq_and_ss(fseq, fss);
    auto chip_motif_infos = get_motifs_from_seq_and_ss(cseq, css);
    mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[1].name),
                  -1, -1, "A7-A22");
    
    for(int i = 2; i < flow_motif_infos.size(); i++) {
        if(flow_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(flow_motif_infos[i].name,
                                                                   flow_motif_infos[i].end_id,
                                                                   flow_motif_infos[i].end_name));
            
        }
    }
    
    mt->add_motif(ResourceManager::getInstance().get_motif("GAAA_tetraloop", "", "A149-A154"));
    
    auto fixed_motif_infos = MotifInfos();
    for (int i = chip_motif_infos.size()-1; i > -1; i--) {
        auto spl = split_str_by_delimiter(chip_motif_infos[i].name, "=");
        auto new_info = chip_motif_infos[i];
        new_info.name = spl[1].substr(1) + spl[1].substr(0,1) + "=" +
                        spl[0].substr(1) + spl[0].substr(0,1);
        //std::cout << chip_motif_infos[i].name << " " << new_info.name << std::endl;
        fixed_motif_infos.push_back(new_info);
    }
    
    mt->add_motif(ResourceManager::getInstance().get_motif(fixed_motif_infos[0].name),
                  pos, -1, "A1-A6");
    
    for(int i = 1; i < fixed_motif_infos.size()-1; i++) {
        if(fixed_motif_infos[i].end_id.length() == 0) {
            mt->add_motif(ResourceManager::getInstance().get_motif(fixed_motif_infos[i].name));
        }
        else{
            mt->add_motif(ResourceManager::getInstance().get_motif(fixed_motif_infos[i].name,
                                                                   fixed_motif_infos[i].end_id,
                                                                   fixed_motif_infos[i].end_name));
            
        }
    }

    
    
    MotifStateEnsembleTreeOP mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);
    
    return mset;

}


MotifInfos
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

    auto motif_infos = MotifInfos();
    
    String motif_name, motif_name_rna;
    String seq1, seq2;
    for(int i = 1; i < required_nodes.size(); i++) {
        if(required_nodes[i]->data()->type()   != sstruct::SS_NodeData::SS_Type::SS_BP ||
           required_nodes[i-1]->data()->type() != sstruct::SS_NodeData::SS_Type::SS_BP) {
            String seq1, seq2, ss1, ss2;
            seq1 = required_nodes[i-1]->data()->ss_chains()[0]->sequence() +
            required_nodes[i]->data()->ss_chains()[0]->sequence()   +
            required_nodes[i+1]->data()->ss_chains()[0]->sequence();
            
            seq2 = required_nodes[i+1]->data()->ss_chains()[1]->sequence() +
            required_nodes[i]->data()->ss_chains()[1]->sequence()   +
            required_nodes[i-1]->data()->ss_chains()[1]->sequence();
            
            ss1 = "L";
            for(int j = 0; j < required_nodes[i]->data()->ss_chains()[0]->length(); j++) {
                ss1 += "U";
            }
            ss1 += "L";
            ss2 = "R";
            for(int j = 0; j < required_nodes[i]->data()->ss_chains()[1]->length(); j++) {
                ss2 += "U";
            }
            ss2 += "R";
            String end_id = seq1 + "_" + ss1 + "_" + seq2 + "_" + ss2;
            auto m = ResourceManager::getInstance().get_motif("", end_id);
            motif_infos.push_back(MotifInfo{m->name(), m->end_ids()[0], m->ends()[0]->name()});
            i += 1;
            continue;
            
            //throw std::runtime_error("old method does not have non helical motifs implemented yet!!!!!");
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
        motif_infos.push_back(MotifInfo{motif_name_rna, "", ""});
        
    }
    
    return motif_infos;
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