#ifndef RNAMAKE_HELIX_1_H
#define RNAMAKE_HELIX_1_H

#include <rnamake2d/strategy2d.h>


namespace rnamake2d {

    class Helix1 : public Strategy2D {
        float base_ = 92.02556323;
        Floats params_;
    public:
        Helix1() : Strategy2D() {
            name_ = "Helix1";
            params_ = {
                -1.96421263, // au/ua mult
                -0.92136502, // gc/cg mult
                -3.86915512  // ug/gu mult
            };
        }

    public:
        float
        score(Feature2DOP const & feature ) override {
            auto num_au(0);
            auto num_gc(0);
            auto num_ug(0);

            for(const auto& motif : feature->motifs) {
                if(motif->token() != "Helix1") {
                    continue;
                }
                const auto& seq = motif->sequence;
                if(seq == "U&A" || seq == "A&U") {
                    ++num_au;
                } else if( seq == "G&C" || seq == "C&G") {
                    ++num_gc;
                } else if (seq == "U&G" || seq == "G&U") {
                    ++num_ug;
                }
            }

            return base_ + num_au*params_[0] + num_gc*params_[1] + num_ug*params_[2];
        }
    };
}


#endif // RNAMAKE_HELIX_1_H
