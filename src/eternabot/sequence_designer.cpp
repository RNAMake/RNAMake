//
//  sequence_designer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <base/log.h>
#include <secondary_structure/sequence_tools.h>
#include <eternabot/sequence_designer.h>

namespace eternabot {


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// setup functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SequenceDesigner::SequenceDesigner():
        scorer_(Scorer()),
        results_(SequenceDesignerResultOPs()),
        rng_(util::RandomNumberGenerator()) {

    parameters_.biased_gc_caps = true;
    v_ = vienna::Vienna();

    possible_bps_ = std::vector<Strings>({{"A", "U"}, {"U", "A"}, {"G", "C"}, {"C", "G"}, {"G", "U"}, {"U", "G"}});
    possible_rt_bps_ = std::vector<secondary_structure::ResTypes>{
            {secondary_structure::ResType::A, secondary_structure::ResType::U},
            {secondary_structure::ResType::U, secondary_structure::ResType::A},
            {secondary_structure::ResType::C, secondary_structure::ResType::G},
            {secondary_structure::ResType::G, secondary_structure::ResType::C},
            {secondary_structure::ResType::U, secondary_structure::ResType::G},
            {secondary_structure::ResType::G, secondary_structure::ResType::U}};

    possible_res_types_ = std::vector<secondary_structure::ResType> {
            secondary_structure::ResType::A,
            secondary_structure::ResType::C,
            secondary_structure::ResType::G,
            secondary_structure::ResType::U
    };


    // generate constraints
    auto disallowed_sequences = Strings{"AAAA", "CCCC", "GGGG", "UUUU"};
    for(auto const & seq : disallowed_sequences) { seq_constraints_.add_disallowed_sequence(seq); }
    seq_constraints_.add_gc_helix_stretch_limit(3);
    current_violations_ = Ints(seq_constraints_.num_constraints());
    next_violations_    = Ints(seq_constraints_.num_constraints());

    setup_options();
    temperature_ = 4.0f;
    mc_ = util::MonteCarlo(temperature_);
}


void
SequenceDesigner::setup() {}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// option functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SequenceDesigner::setup_options() {
    options_.add_option("designs", 1, base::OptionType::INT);
    options_.add_option("steps", 1000, base::OptionType::INT);
    options_.lock_option_adding();
    update_var_options();
    
}
    
void
SequenceDesigner::update_var_options() {
    designs_ = options_.get_int("designs");
    steps_   = options_.get_int("steps");
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


SequenceDesignerResultOPs const &
SequenceDesigner::design(
        secondary_structure::PoseOP const & p) {

    results_ = SequenceDesignerResultOPs();

    pair_map_ = std::vector<std::vector<int>>(p->residues().size()+1);
    for(int i = 0; i < p->residues().size()+1; i++) {
        pair_map_[i] = std::vector<int>(p->residues().size()+1);
    }
    for(auto const & bp : p->basepairs()) {
        pair_map_[bp->res1()->num()][bp->res2()->num()] = 1;
        pair_map_[bp->res2()->num()][bp->res1()->num()] = 1;
    }
    pair_map_entries_ = p->residues().size()*p->residues().size();
    res_type_constraints_ = std::map<int, secondary_structure::ResType>();
    auto motifs = p->helices();

    // boosting!
    for(auto const & m : p->motifs()) {
        if(m->mtype() == util::MotifType::HELIX) { continue; }
        for(auto const & end : m->ends()) {
            closing_bps_[end] = 1;
        }

        if(m->mtype() == util::MotifType::HAIRPIN) {
            int i = -1;
            for(auto const & r : m->residues()) {
                i++;
                if(i == 0 || i == m->residues().size()-1) { continue; }
                if(i == 1 && m->residues().size() > 5) {
                    if(secondary_structure::is_restype_a_ambiguous_code(r->res_type())) {
                        if(secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::A, r->res_type()) &&
                           secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::G, r->res_type())) {
                            res_type_constraints_[r->num()] = secondary_structure::ResType::R;
                        }
                        else if(secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::A, r->res_type())) {
                            res_type_constraints_[r->num()] = secondary_structure::ResType::A;
                        }
                        else if(secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::G, r->res_type())) {
                            res_type_constraints_[r->num()] = secondary_structure::ResType::G;
                        }
                    }
                }
                else if(i == 1) {
                    if(secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::A, r->res_type())) {
                        res_type_constraints_[r->num()] = secondary_structure::ResType::A;
                    }
                }
                else {
                    if (secondary_structure::is_restype_a_ambiguous_code(r->res_type())) {
                        if (secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::A, r->res_type())) {
                            res_type_constraints_[r->num()] = secondary_structure::ResType::A;
                        }
                    }
                }
            }
        }
        else if(m->mtype() == util::MotifType::TWOWAY) {
            // single mismatch boost to G-G
            if(m->residues().size() == 6) {
                for(auto const & r : m->residues()) {
                    auto bps = p->get_basepair(r->uuid());
                    if(bps.size() != 0) { continue; }
                    if(secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::G, r->res_type())) {
                        res_type_constraints_[r->num()] = secondary_structure::ResType::G;
                    }
                }
            }

            else {
                for(auto const & c : m->chains()) {
                    int i = -1;
                    for(auto const & r : c->residues()) {
                        i++;
                        if(i < 2 || i > c->residues().size()-3) { continue;}
                        if (secondary_structure::is_restype_a_ambiguous_code(r->res_type())) {
                            if (secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::A, r->res_type())) {
                                res_type_constraints_[r->num()] = secondary_structure::ResType::A;
                            }
                        }
                    }
                }
            }

        }

        else {
            for(auto const & c : m->chains()) {
                int i = -1;
                for(auto const & r : c->residues()) {
                    i++;
                    if(i < 2 || i > c->residues().size()-3) { continue;}
                    if (secondary_structure::is_restype_a_ambiguous_code(r->res_type())) {
                        if (secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::A, r->res_type())) {
                            res_type_constraints_[r->num()] = secondary_structure::ResType::A;
                        }
                    }
                }
            }
        }

        //std::cout << m->sequence() << std::endl;
    }

    //std::cout << p->sequence() << std::endl;
    for(auto const & r : p->residues()) {
        if(res_type_constraints_.find(r->num()) == res_type_constraints_.end()) {
            res_type_constraints_[r->num()] = r->res_type();
        }
        //std::cout << secondary_structure::convert_res_type_to_str(res_type_constraints_[r->num()]);
    }
    //std::cout << std::endl;
    //std::cout << p->dot_bracket() << std::endl;

    //exit(0);

    _find_designable_bps(p);       // find all basepairs with N-N residues
    for(auto const & r : p->residues()) {
        if ( is_restype_a_ambiguous_code(r->res_type())) {
            designable_res_.push_back(r);
            auto bps = p->get_basepair(r->uuid());
            if(bps.size() == 0) {
                designable_unpaired_res_.push_back(r);
            }

            _assign_new_residue_restype(r);
        }
    }

    scorer_.setup(p);

    for(auto const & h : motifs) {
        _set_initial_helix_sequence(h);
    }

    // get all non bp-step motifs
    for(auto const & m : p->motifs()) {
        if(m->mtype() == util::MotifType::HELIX) { continue; }
        motifs.push_back(m);
    }

    // nothing to design! return just the score
    if(designable_unpaired_res_.size() == 0 && designable_bps_.size() == 0) {
        auto current_score = scorer_.score_secondary_structure(p);
        auto bp_diff_score = _bp_list_diff(p, pair_map_, pair_map_entries_, scorer_.features());
        results_.push_back(std::make_shared<SequenceDesignerResult>(p->sequence(), current_score, bp_diff_score));
        //std::cout << scorer_.print_scores(p) << std::endl;

        return results_;
    }

    //std::cout << _bp_list_diff(p, pair_map_, pair_map_entries_) << std::endl;
    //std::cout << p->sequence() << std::endl;
    //std::cout << _bp_list_diff(p, pair_map_, pair_map_entries_) << std::endl;
    auto score = _optimize_substructure(p, steps_);

    //std::cout << p->sequence() << " " << score << std::endl;
    //std::cout << _bp_list_diff(p, pair_map_, pair_map_entries_) << std::endl;
    score = scorer_.score_secondary_structure(p);
    //std::cout << scorer_.print_scores(p) << std::endl;
    auto bp_diff_score = _bp_list_diff(p, pair_map_, pair_map_entries_, scorer_.features());

    if(score > 100) { score = 100; }


    results_.push_back(std::make_shared<SequenceDesignerResult>(p->sequence(), score, bp_diff_score));
    return results_;

    exit(0);

    return results_;
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
SequenceDesigner::_find_designable_bps(
        secondary_structure::PoseOP p) {
    for(auto const & bp : p->basepairs()) {
        if(is_restype_a_ambiguous_code(bp->res1()->res_type()) || is_restype_a_ambiguous_code(bp->res2()->res_type()) ) {
            designable_bps_.push_back(bp);
        }
    }

}

bool
SequenceDesigner::_new_sequence_violations(
        Ints const & current_violations,
        Ints const & next_violations) {
    for(int i = 0; i < current_violations.size(); i++) {
        if(current_violations[i] < next_violations[i]) { return true; }
    }
    return false;
}

void
SequenceDesigner::_set_initial_helix_sequence(
        secondary_structure::MotifOP h) {
    auto n_bps = h->residues().size()/2;
    auto designable_1 = Ints(n_bps);
    auto designable_2 = Ints(n_bps);
    int i = 0;
    for(auto const & r : h->chains()[0]->residues()) {
        if(std::find(designable_res_.begin(), designable_res_.end(), r) != designable_res_.end()) {
            designable_1[i] = 1;
        }
        else {
            designable_1[i] = 0;
        }
        i++;
    }

    for(auto const & r : h->chains()[1]->residues()) {
        i--;
        if(std::find(designable_res_.begin(), designable_res_.end(), r) != designable_res_.end()) {
            designable_2[i] = 1;
        }
        else {
            designable_2[i] = 0;
        }
    }

    auto res_types_1 = secondary_structure::ResTypes(n_bps);
    auto res_types_2 = secondary_structure::ResTypes(n_bps);
    auto res_type_bp = secondary_structure::ResTypes(2);

    auto gc_count = 0;
    auto c_violations = seq_constraints_.violations(h);
    auto n_violations = Ints(c_violations);
    auto iter = 0;
    auto bp_iter = 0;
    while(true) {
        i = 0;
        for(auto & bp : h->basepairs()) {
            bp_iter = 0;
            auto org_res_type_1 = res_type_constraints_[bp->res1()->num()];
            auto org_res_type_2 = res_type_constraints_[bp->res2()->num()];
            while(true) {
                bp_iter += 1;
                if(bp_iter > 100000) {
                    LOG_ERROR << "cannot find a basepair that satisfies constraints: " +
                    secondary_structure::convert_res_type_to_str(org_res_type_1) + " " +
                    secondary_structure::convert_res_type_to_str(org_res_type_2);
                    exit(0);
                }
                // can design both sides of the helix
                if (designable_1[i] && designable_2[i]) {
                    if ((i == 0 || i == n_bps - 1) && parameters_.biased_gc_caps) {
                        _get_random_res_type_pair_gc_cap(res_type_bp);
                    } else {
                        _get_random_res_type_pair(res_type_bp);
                    }
                }
                    // cannot design either side
                else if (designable_1[i] == 0 && designable_2[i] == 0) {
                    res_type_bp[0] = bp->res1()->res_type();
                    res_type_bp[1] = bp->res2()->res_type();
                }
                    // cannot design first residue
                else if (designable_1[i] == 0) {
                    res_type_bp[0] = bp->res1()->res_type();
                    res_type_bp[1] = secondary_structure::get_complement_res_type(bp->res1()->res_type());
                }
                    // cannot design second residue
                else if (designable_2[i] == 0) {
                    res_type_bp[0] = secondary_structure::get_complement_res_type(bp->res2()->res_type());
                    res_type_bp[1] = bp->res2()->res_type();
                }


                if(!secondary_structure::does_restype_satisfy_constraint(res_type_bp[0], org_res_type_1) ||
                   !secondary_structure::does_restype_satisfy_constraint(res_type_bp[1], org_res_type_2)) {
                    continue;
                }

                res_types_1[i] = res_type_bp[0];
                res_types_2[i] = res_type_bp[1];

                bp->res1()->res_type(res_type_bp[0]);
                bp->res2()->res_type(res_type_bp[1]);

                break;

            }
            i += 1;
        }

        iter += 1;
        if(iter > 1000) {
            LOG_WARNING << "could not find suitable starting sequence for helix";
            break;
        }

        n_violations = seq_constraints_.violations(h);
        if(_new_sequence_violations(c_violations, n_violations)) {
            continue;
        }

        break;


    }

    //std::cout << h->sequence() << std::endl;


}

void
SequenceDesigner::_get_random_res_type_pair_gc_cap(
        secondary_structure::ResTypes & pair) {
    // biased base pair selection for caped ends, 80% chance to be GC/CG over AU/UA
    // will do GCs
    auto rand = rng_.randrange(1000);
    if(rand > 200) {
        // selecting GC
        if(rng_.randrange(1000) > 500) {
            pair[0] = secondary_structure::ResType::G;
            pair[1] = secondary_structure::ResType::C;
        }
        // selecting CG
        else {
            pair[0] = secondary_structure::ResType::C;
            pair[1] = secondary_structure::ResType::G;
        }
    }

    else if(rand > 50){
        // selecting AU
        if(rng_.randrange(1000) > 500) {
            pair[0] = secondary_structure::ResType::A;
            pair[1] = secondary_structure::ResType::U;
        }
        // selecting UA
        else {
            pair[0] = secondary_structure::ResType::U;
            pair[1] = secondary_structure::ResType::A;
        }
    }
    else {
        // selecting AU
        if(rng_.randrange(1000) > 500) {
            pair[0] = secondary_structure::ResType::G;
            pair[1] = secondary_structure::ResType::U;
        }
        // selecting UA
        else {
            pair[0] = secondary_structure::ResType::U;
            pair[1] = secondary_structure::ResType::G;
        }

    }
}

void
SequenceDesigner::_get_random_res_type_pair(
        secondary_structure::ResTypes & pair) {
    if(rng_.randrange(1000) > 200) {
        pair = possible_rt_bps_[rng_.randrange(4)];
    }
    else {
        pair = possible_rt_bps_[rng_.randrange(6)];
    }
}


float
SequenceDesigner::_optimize_substructure(
        secondary_structure::PoseOP p,
        int steps) {


    auto scorer = Scorer();
    scorer.setup(p);
    auto current_designable_bps = secondary_structure::BasepairOPs();
    auto current_designable_unpaired_res = secondary_structure::ResidueOPs();
    for(auto const & bp : p->basepairs()) {
        if(std::find(designable_bps_.begin(), designable_bps_.end(), bp) != designable_bps_.end()) {
            current_designable_bps.push_back(bp);
        }
    }
    for(auto const & r : p->residues()) {
        if(std::find(designable_unpaired_res_.begin(), designable_unpaired_res_.end(), r) != designable_unpaired_res_.end()) {
            current_designable_unpaired_res.push_back(r);
        }
    }

    auto current_move = MonteCarloMoveOP(nullptr);
    auto moves = MonteCarloMoveOPs(2);
    moves[0] = std::make_shared<MutateBPMove>(current_designable_bps, possible_rt_bps_, res_type_constraints_, closing_bps_);
    moves[1] = std::make_shared<MutateUnpairedResMove>(current_designable_unpaired_res, res_type_constraints_);

    auto current_sequence = p->sequence();
    auto next_score = 0.0f;
    //auto best_score = scorer.score_secondary_structure(p);
    auto eternabot_score = 0.0f;
    auto best_score = exp(-_bp_list_diff(p, pair_map_, pair_map_entries_, scorer.features())/10)*scorer.score_secondary_structure(p);
    auto current_score = best_score;

    auto best_sequence = p->sequence();
    auto pos = 0;
    auto tried_again_count = 0;
    auto length = p->sequence().length();

    for(int i = 0; i < steps; i++) {
        if(tried_again_count > 10) {
            tried_again_count = 0;
            i++;
        }
        pos = rng_.randrange(1000);
        if(pos < 500) {
            current_move = moves[0];
        }
        else {
            current_move = moves[1];
        }
        // move didn't do anything, try again
        if(current_move->move(p) == 0) {
            i--;
            tried_again_count++;
            continue;
        }
        tried_again_count = 0;

        //next_score = scorer_.score_secondary_structure(p);
        eternabot_score = scorer.score_secondary_structure(p);
        next_score = exp(-_bp_list_diff(p, pair_map_, pair_map_entries_, scorer.features())/10)*eternabot_score;

        if(mc_.accept(next_score, current_score)) {
            current_score = next_score;
        }
        else {
            current_move->undo(p);
        }

        if(current_score > best_score) {
            auto found = 0;
            for(auto const & previous_sol : previous_solutions_) {
                if(previous_sol == p->sequence()) { found = 1; break; }
            }
            if(found) { continue; }

            //std::cout << best_score << " " << _bp_list_diff(p, pair_map_, pair_map_entries_, scorer.features()) << " " << exp(-_bp_list_diff(p, pair_map_, pair_map_entries_, scorer.features())/5) <<" " << eternabot_score << std::endl;
            best_score = current_score;
            best_sequence = p->sequence();
        }
    }

    p->replace_sequence(best_sequence);

    return best_score;
}

float
SequenceDesigner::_bp_list_diff(
        secondary_structure::PoseOP p,
        std::vector<std::vector<int>> const & pair_map,
        size_t pair_list_size,
        FeaturesOP features) {

    int pi = 0, pj = 0;
    auto score = 0.0f;
    auto & plist = features->dotplot;
    for(int i = 0 ; i < pair_list_size; i++) {
        if(plist[i].p < 0.001) { continue; }
        pi = plist[i].i;
        pj = plist[i].j;
        //score += abs(pair_map[pi][pj] - plist[i].p);
    }

    int i = -1;
    for(auto const & r : p->residues()) {
        i++;
        if(r->dot_bracket()[0] != features->structure[i]) { score += 1;}
    }

    return score;


}

bool
SequenceDesigner::_assign_new_residue_restype(
        secondary_structure::ResidueOP r) {
    auto org_res_type = res_type_constraints_[r->num()];
    auto count = 0;
    while(1) {
        if(std::find(designable_unpaired_res_.begin(), designable_unpaired_res_.end(), r) != designable_unpaired_res_.end()) {
            if(secondary_structure::does_restype_satisfy_constraint(secondary_structure::ResType::A, org_res_type)) {
                r->res_type(secondary_structure::ResType::A);
                return true;
            }
        }

        auto res_type = possible_res_types_[rng_.randrange(possible_res_types_.size())];
        if(secondary_structure::does_restype_satisfy_constraint(res_type, org_res_type)) {
            r->res_type(res_type);
            return true;
        }

        count += 1;
        if(count > 1000) { break; }
    }

    return false;
}



}






















