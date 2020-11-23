#ifndef RNAMAKE_BADHAIRPIN4_H
#define RNAMAKE_BADHAIRPIN4_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class BadHairpin4 : public Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                "UGAGAG",
                "AGAAAU",
                "AAGGAU",
                "CAAAAG",
                "CGAAAG",
                "UCCAGA",
                "GUAAUU",
                "GGAAAU",
                "UUUCGA",
                "UUUCGG",
                "AUUCGU",
                "UGAAAG",
                "UAUAAA",
                "CUUCGG",
        };
    public:
        BadHairpin4() : Strategy2D() {
            params_ = {1., 1., 1.};
            name_ = "BadHairpin4";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(const auto& hp : feature->motifs) {
                if(hp->token() == "hairpin4" && bad_.find(hp->sequence) != bad_.end()) {
                    auto pair = String(2,' ');
                    pair[0] = hp->sequence[0];
                    pair[1] = *hp->sequence.rbegin();
                    if(pair == "AU" || pair == "UA") {
                        result -= params_[0];
                    } else if (pair == "GC" || pair == "CG") {
                        result -= params_[1];
                    } else if (pair == "UG" || pair == "GU") {
                        result -= params_[2];
                    }
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_BADHAIRPIN4_H 