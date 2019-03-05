//
// Created by Joseph Yesselman on 3/3/19.
//

#include "build_flex_helices/build_flex_helices.h"


#include "base/backtrace.hpp"
#include "util/cartesian_product.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_state_tree.h"
#include "resources/resource_manager.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"

BuildFlexHelicesApp::BuildFlexHelicesApp() {
    disallowed_num_sequences_ = std::vector<Ints>();
    for(int i = 0; i < 4; i++) {
        auto disallowed_sequence = Ints(4);
        for(int j = 0; j < 4; j++) { disallowed_sequence[j] = i; }
        disallowed_num_sequences_.push_back(disallowed_sequence);
    }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
BuildFlexHelicesApp::setup_options() {

}

void
BuildFlexHelicesApp::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);
}

void
BuildFlexHelicesApp::run() {

    auto ideal = RM::instance().motif("HELIX.LE.2");
    ideal->to_pdb("ideal_le.pdb", 1, 1);

    RM::instance().motif("HELIX.IDEAL.2");
    ideal->to_pdb("ideal.pdb", 1, 1);

    get_avg_helix(4);

    /*auto seq = String("CCCCC&GGGGG");
    auto ss  = String("(((((&)))))");

    auto motifs = get_motifs_from_seq_and_ss(seq, ss);
    auto mst = std::make_shared<MotifStateTree>();
    for(auto const & m : motifs) {
        mst->add_state(m);
    }

    mst->write_pdbs();*/

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MotifStateOPs
BuildFlexHelicesApp::get_motifs_from_seq_and_ss(
        String const & seq,
        String const & ss) {

    auto parser = sstruct::SecondaryStructureParser();
    auto ss_motifs = parser.parse_to_motifs(seq, ss);
    auto motifs = MotifStateOPs();

    auto start = 0;
    auto motif = MotifStateOP(nullptr);
    for(auto const & m : ss_motifs) {
        //basepair step
        if(m->mtype() == MotifType::HELIX) {
            motif = RM::instance().bp_step(m->end_ids()[0])->get_state();
            motifs.push_back(motif);
        }
        else {
            throw std::runtime_error("only helices are allowed");
        }
    }

    return motifs;
}

void
BuildFlexHelicesApp::generate_helices(
        int length) {

    auto pairs = std::vector<Strings>();
    pairs.push_back(Strings{"A", "U"});
    pairs.push_back(Strings{"U", "A"});
    pairs.push_back(Strings{"C", "G"});
    pairs.push_back(Strings{"G", "C"});

    auto all_pairs = std::vector<std::vector<Strings>>((int)(length));
    for(int i = 0; i < all_pairs.size(); i++) {
        all_pairs[i] = pairs;
    }

    auto pair_iterator = CartesianProduct<Strings>(all_pairs);
    auto current = std::vector<Strings>();
    auto new_seq = String();
    auto seq1 = String();
    auto seq2 = String();
    auto structure = String();

    for(int i = 0; i < length; i++) { structure += "("; }
    structure += "&";
    for(int i = 0; i < length; i++) { structure += ")"; }

    auto motifs = MotifStateOPs();

    std::ofstream out;
    out.open("length_"+std::to_string(length)+"_helices.dat");

    int i = 0;
    while (!pair_iterator.end()) {
        current = pair_iterator.next();
        seq1 = "";
        seq2 = "";
        for (auto const & p : current) {
            seq1 += p[0];
            seq2 = p[1] + seq2;
        }
        new_seq = seq1 + "&" + seq2;
        motifs = get_motifs_from_seq_and_ss(new_seq, structure);

        auto mst = std::make_shared<MotifStateTree>();
        for(auto const & m : motifs) {
            mst->add_state(m);
        }
        /*mst->to_motif_tree()->to_pdb("helix."+std::to_string(i)+".pdb");
        i++;
        if(i > 100) {
            break;
        }*/
        mst->remove_node_level();

        i++;
        //std::cout << new_seq << " " <<  _get_hit_count(new_start, new_seq, struc) << std::endl;
    }
    std::cout << i << std::endl;

}

MotifOP
BuildFlexHelicesApp::get_avg_helix(
        int length) {

    auto pairs = std::vector<Strings>();
    pairs.push_back(Strings{"A", "U"});
    pairs.push_back(Strings{"U", "A"});
    pairs.push_back(Strings{"C", "G"});
    pairs.push_back(Strings{"G", "C"});

    auto all_pairs = std::vector<std::vector<Strings>>((int)(length));
    for(int i = 0; i < all_pairs.size(); i++) {
        all_pairs[i] = pairs;
    }

    auto pair_iterator = CartesianProduct<Strings>(all_pairs);
    auto current = std::vector<Strings>();
    auto new_seq = String();
    auto seq1 = String();
    auto seq2 = String();
    auto structure = String();

    for(int i = 0; i < length; i++) { structure += "("; }
    structure += "&";
    for(int i = 0; i < length; i++) { structure += ")"; }

    auto motifs = MotifStateOPs();

    std::ofstream out;
    out.open("length_"+std::to_string(length)+"_helices.dat");


    auto num_seq = Ints(length);

    int i = 0;
    while (!pair_iterator.end()) {
        current = pair_iterator.next();
        seq1 = "";
        seq2 = "";
        for (auto const & p : current) {
            seq1 += p[0];
            seq2 = p[1] + seq2;
        }
        for(int j = 0; j < length; j++) {
            num_seq[j] = convert_char_to_res_code(seq1[j]);
        }

        new_seq = seq1 + "&" + seq2;
        /*motifs = get_motifs_from_seq_and_ss(new_seq, structure);

        auto mst = std::make_shared<MotifStateTree>();
        for(auto const & m : motifs) {
            mst->add_state(m);
        }*/


        /*mst->to_motif_tree()->to_pdb("helix."+std::to_string(i)+".pdb");
        i++;
        if(i > 100) {
            break;
        }*/

        i++;
        //std::cout << new_seq << " " <<  _get_hit_count(new_start, new_seq, struc) << std::endl;
    }
    std::cout << i << std::endl;
}

bool
BuildFlexHelicesApp::find_seq_violations(
        Ints const & num_seq) {
    auto pos = 0;

    for (auto const & r : num_seq) {
        for (int j = 0; j < disallowed_res_types_sequences_.size(); j++) {
            if (i < disallowed_res_types_sequences_[j].size()) { continue; }
            auto match = true;
            pos = i - disallowed_res_types_sequences_[j].size();
            for (auto const & e : disallowed_res_types_sequences_[j]) {
                if (res[pos]->res_type() != e) {
                    match = false;
                    break;
                }
                pos += 1;
            }
            if (match) {
                violations[j] += 1;
            }
        }
    }
    return false;
}



int
BuildFlexHelicesApp::convert_char_to_res_code(
        char c) {
    if     (c == 'A') { return 0; }
    else if(c == 'C') { return 1; }
    else if(c == 'G') { return 2; }
    else if(c == 'U') { return 3; }
    else if(c == 'T') { return 3; }
    else if(c == 'N') { return -1; }
    else {
        throw BuildFlexHelicesAppException("incorrect character for secondary string");
    }
    return -1;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(print_backtrace);

    //load extra motifs being used
    //String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/apt_new_interface/resources/";
    //RM::instance().add_motif(base_path+"pRNA_3WJ.pdb", "prna");

    auto app = BuildFlexHelicesApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}