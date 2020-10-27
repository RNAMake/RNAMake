#include <rnamake2d/design.h>


namespace rnamake2d {
    int Design::design_ct = 1;


    void
    Design::initialize_features() {
        best_pose = parser.parse_to_pose(sequence,target);
        best = std::make_shared<rnamake2d::Feature2D>(*generator.eternabot::FeatureGenerator::get_features(best_pose));
        // also setup the motifs here
        best->motifs = parse_to_motif2ds(best_pose, target);
        for(auto& motif : best->motifs) {
            motif->full_sequence( sequence );
        }
        best->pairmap = get_pairmap_from_secstruct( target );
        best->elements = get_rna_elemnts_from_secstruct( target );
        generator.update_features(best, best_pose);
    }

    void
    Design::initialize_mutant() {
        mutant_pose = parser.parse_to_pose(candiate,target);
        mutant = std::make_shared<rnamake2d::Feature2D>(*generator.eternabot::FeatureGenerator::get_features(mutant_pose));
        // also setup the motifs here
        mutant->motifs = parse_to_motif2ds(best_pose, target);
        //for(auto& m : mutant->motifs)  m->show();
        for(auto& motif : mutant->motifs) {
            motif->full_sequence( candiate );
        }
        mutant->pairmap = get_pairmap_from_secstruct( target );
        mutant->elements = get_rna_elemnts_from_secstruct(target);
        generator.update_features(mutant, mutant_pose);
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
