#include <rnamake2d/Design.h>


namespace rnamake2d {
    int Design::design_ct = 0 ;


    void
    Design::initialize_features() {
        best_pose = parser.parse_to_pose(sequence,target);
        best = std::make_shared<rnamake2d::Feature2D>(*(generator.eternabot::FeatureGenerator::get_features(best_pose)));
        // also setup the motifs here
        best->motifs = parse_to_motif2ds(best_pose, target);
        for(auto& m : best->motifs)  m->show();
        auto abt = eternabot::ABasicTest();
    }

    void
    Design::initialize_mutant() {
        mutant_pose = parser.parse_to_pose(candiate,target);
        mutant = std::make_shared<rnamake2d::Feature2D>(*(generator.eternabot::FeatureGenerator::get_features(mutant_pose)));
        // also setup the motifs here
        mutant->motifs = parse_to_motif2ds(best_pose, target);
        for(auto& m : mutant->motifs)  m->show();
        auto abt = eternabot::ABasicTest();
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
        designs.emplace_back(Design(target,sequence));
        const auto last = designs.rbegin();

        for(auto ii = 1; ii < number ; ++ii) {
            designs.push_back(*designs.rbegin());
        }
        return designs;
    }
}
