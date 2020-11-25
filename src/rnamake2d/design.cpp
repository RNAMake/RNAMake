#include <rnamake2d/design.h>


namespace rnamake2d {
    int Design::design_ct = 1;

    void
    Design::initialize_features() {
        best_pose = parser.parse_to_pose(sequence_, target_);
        best = std::make_shared<rnamake2d::Feature2D>(*generator.eternabot::FeatureGenerator::get_features(best_pose));
        // also setup the motifs here
        best->motifs = parse_to_motif2ds(best_pose, target_);
        for(auto& motif : best->motifs) {
            motif->full_sequence( sequence_);
        }
        best->e_pairmap = get_pairmap_from_secstruct( target_ );
        best->elements = get_rna_elemnts_from_secstruct( target_ );
        generator.update_features(best, best_pose);
    }

    void
    Design::accept(double mutant_score) {
        score_ = mutant_score;
        scores.push_back(mutant_score);
        ++num_moves;
        sequence_ = candiate_;
        best = std::make_shared<Feature2D>(*mutant_);
        mutant_ = nullptr;
        best_pose = mutant_pose;
        mutant_pose = nullptr;
    }

    void
    Design::accept_bp(double mutant_score) {
        bp_score_ = mutant_score;
        //scores.push_back(mutant_score);
        ++num_moves;
        sequence_ = candiate_;
        best = mutant_; // ? not sure on this one tbh
        mutant_ = nullptr;
        best_pose = mutant_pose;
        mutant_pose = nullptr;
    }

    void
    Design::initialize_mutant(String& mutated) {
        candiate_ = mutated;
        //sequence_ = mutated;

        mutant_pose = parser.parse_to_pose(candiate_, target_);
        mutant_ = std::make_shared<rnamake2d::Feature2D>(*generator.eternabot::FeatureGenerator::get_features(mutant_pose));
        // also setup the motifs here
        mutant_->motifs = parse_to_motif2ds(best_pose, target_);
        //for(auto& m : mutant->motifs)  m->show();
        for(auto& motif : mutant_->motifs) {
            motif->full_sequence( candiate_ );
        }
        mutant_->e_pairmap = get_pairmap_from_secstruct( target_ );
        mutant_->sequence = mutated;
        mutant_->target = target_;
        mutant_->elements = get_rna_elemnts_from_secstruct(target_);
        generator.update_features(mutant_, mutant_pose);
    }

    Designs
    make_designs(String target, String sequence, int number) {
        /* what does this function do?
         * 1. given the input of a target and sequence => probably NNNN, get a usable design object back
         * 2. do it "N" or number times
        */
        auto designs = Designs();
        while(target.size() - sequence.size()) {
            sequence += 'N';
        }

        for(auto ii = 0 ; ii < number ; ++ii) {
            designs.emplace_back(target,sequence);
        }
        return designs;
    }
}
