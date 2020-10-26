#ifndef RNAMAKE_GOODHAIRPIN4_H
#define RNAMAKE_GOODHAIRPIN4_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class GoodHairpin4 : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                                    "CAAAAG",
                                    "CGAAAG",
                                    "CUAAUG",
                                    "GAAAAC",
                                    "GGAAAC",
                                    "UCCAGA",
                                    };
    public:
        GoodHairpin4() : Strategy2D() {
            params_ = {91.99024139, -4.94683135, 30.25168973};
            name_ = "GoodHairpin4";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(params_[0]);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() == util::MotifType::HAIRPIN &&
                    motif->sequence.size() == 6   /* check that it is a hairpin4 */
                        ) {
                    if( good_.find(motif->sequence) != good_.end()) {
                        const String end = String{*motif->sequence.begin()} + String{*motif->sequence.rbegin()};

                        if(end == "AU" || end == "UA") {
                            result -= params_[1];
                        } else if (end == "GC" || end == "CG") {
                            result -= params_[2];
                        }
                    }
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_GOODHAIRPIN4_H 