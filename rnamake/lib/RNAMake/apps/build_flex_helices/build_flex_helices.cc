//
// Created by Joseph Yesselman on 3/3/19.
//

#include "build_flex_helices/build_flex_helices.h"


#include "base/backtrace.hpp"
#include "math/quaternion.h"
#include "util/cartesian_product.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_state_tree.h"
#include "resources/resource_manager.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"

BuildFlexHelicesApp::BuildFlexHelicesApp():
        mf_(MotifFactory()){}

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

    auto ideal = RM::instance().motif("HELIX.IDEAL.1");
    ideal->to_pdb("ideal.pdb", 1, 1);

    std::ofstream out;
    out.open("average_helices.dat");

    for(int i = 3; i < 9; i++) {
        std::cout << "GENERATING avg helix with length " << i << std::endl;
        auto m = get_avg_helix(i);
        out << m->to_str() << std::endl;
    }
    out.close();

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

void
BuildFlexHelicesApp::generate_helices(
        int length) {

    /*auto pairs = std::vector<Strings>();
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

        mst->remove_node_level();

        i++;
        //std::cout << new_seq << " " <<  _get_hit_count(new_start, new_seq, struc) << std::endl;
    }
    std::cout << i << std::endl;*/

}

MotifOP
BuildFlexHelicesApp::get_avg_helix(
        int length) {

    auto mst = MotifStateTreeOP(nullptr);
    auto iterator = HelixStructureIterator();
    iterator.setup(length);

    auto q = Quaternion();
    auto q_averager = AverageQuaternionCalculator();
    auto t_average = Point(0, 0, 0);

    int i = 0;
    while(!iterator.end()) {
        mst = iterator.next();

        t_average += mst->last_node()->data()->get_end_state(1)->d();

        q = get_quaternion_from_matrix(mst->last_node()->data()->get_end_state(1)->r());
        q_averager.add_quaternion(q);

        i++;
    }

    t_average = t_average / i;
    auto q_avg = q_averager.get_average();

    auto best_seq = String("");
    auto best_score = 10000.0;

    auto dist = 100.0;
    auto q_dist = 1000.0;
    auto dot = 1000.0;

    iterator = HelixStructureIterator();
    iterator.setup(length);

    while(!iterator.end()) {
        mst = iterator.next();

        dist = t_average.distance(mst->last_node()->data()->get_end_state(1)->d());

        q = get_quaternion_from_matrix(mst->last_node()->data()->get_end_state(1)->r());
        q_dist = 2 * acos(q.dot(q_avg));

        dist += q_dist * 10;

        if (dist < best_score) {
            best_score = dist;
            best_seq = iterator.get_current_sequence();
        }
    }

    auto best_mst = iterator.get_tree_with_sequence(best_seq);
    auto rs = best_mst->to_motif_tree()->get_structure();
    auto m = std::make_shared<Motif>(*rs);
    auto num = 1;
    for(auto & c : m->chains()) {
       for(auto & r : c->residues()) {
           r->num(num);
           r->chain_id("A");
           num++;
       }
    }
    for(auto & bp : m->basepairs()) {
        auto r1 = m->get_residue(bp->res1()->uuid());
        auto r2 = m->get_residue(bp->res2()->uuid());
        bp->res1()->num(r1->num());
        bp->res1()->chain_id("A");
        bp->res2()->num(r2->num());
        bp->res2()->chain_id("A");
    }

    for(auto & bp : m->ends()) {
        auto r1 = m->get_residue(bp->res1()->uuid());
        auto r2 = m->get_residue(bp->res2()->uuid());
        bp->res1()->num(r1->num());
        bp->res1()->chain_id("A");
        bp->res2()->num(r2->num());
        bp->res2()->chain_id("A");
    }

    auto scorer = MotifScorer();

    m->mtype(MotifType::HELIX);
    m->name("HELIX.AVG."+std::to_string(length));
    mf_._setup_secondary_structure(m);
    m->block_end_add(0);
    auto score = scorer.score(m);
    m->score(score);

    //auto s = m->to_str();
    //auto m2 = std::make_shared<Motif>(s, ResidueTypeSetManager::getInstance().residue_type_set());

    //std::cout << i << std::endl;

    return m;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    //load extra motifs being used
    //String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/apt_new_interface/resources/";
    //RM::instance().add_motif(base_path+"pRNA_3WJ.pdb", "prna");

    auto app = BuildFlexHelicesApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}