#ifndef __BAD_HAIRPIN4_DOUBLET_H__
#define __BAD_HAIRPIN4_DOUBLET_H__

#include <unordered_set>

#include <base/types.h>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadHairpin4Doublet : public Strategy2D {
        std::unordered_set<String> bad_ = {
                "GAAAAC",
                "CAUAAG",
                "UGAGAG",
                "CUUCGG",
                "UGAAAA",
                "AGAAAU",
                "UGAAAG",
                "CAAAAG",
                "AGACGU",
                "CGUAAG",
                "UAAAGA",
                "CUCCUG",
                "GGUGAC"
        };
        float penalty_ = 9.74571441;

    public:
        BadHairpin4Doublet() : Strategy2D() {
            name_ = "BadHairpin4Doublet";
        }

    public:
        float
        score( Feature2DOP const & feature ) override {

            auto bad_ct(0);
            for(const auto& motif : feature->motifs) {
                if(motif->token() != "Hairpin4" || motif->buffer() != 2)  {
                    continue;
                }
                if (bad_.find(motif->sequence) != bad_.end()) {
                    ++bad_ct;
                }
            }

            return 100.f - bad_ct*penalty_;
        }
    };

}

#endif// __BAD_HAIRPIN4_DOUBLET_H__
