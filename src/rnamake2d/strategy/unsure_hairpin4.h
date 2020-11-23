#ifndef RNAMAKE_UNSUREHAIRPIN4_H
#define RNAMAKE_UNSUREHAIRPIN4_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class UnsureHairpin4 : public Strategy2D {
    private:
        std::unordered_set<String> unsure_ = {
                "CGAGAG",
                "CUUCGG",
                "GGAAAU",
                "GGAGAC",
                "GUUCGC",
                "UGAAAA",
                "UGAGAG",
                "UUUCGA",
        };
    public:
        UnsureHairpin4() : Strategy2D() {
            params_ = { 1., 1., 1.};
            name_ = "UnsureHairpin4";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(const auto& hp : feature->motifs) {
                if(hp->mtype() != util::MotifType::HAIRPIN) {
                    continue;
                }
                if(unsure_.find(hp->sequence) != unsure_.end()) {
                    auto pair = String(2,' ');
                    pair[0] = hp->sequence[0];
                    pair[1] = *hp->sequence.rbegin();
                    const float length = hp->sequence.size();
                    if(pair == "AU" || pair == "UA") {
                        result -= params_[0]/length;
                    } else if (pair == "GC" || pair == "CG") {
                        result -= params_[1]/length;
                    } else if (pair == "UG" || pair == "GU") {
                        result -= params_[2]/length;
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_UNSUREHAIRPIN4_H 