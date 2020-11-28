#ifndef RNAMAKE_GOOD_HAIRPIN4_H
#define RNAMAKE_GOOD_HAIRPIN4_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class GoodHairpin4Singlet : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                                "CUAGCG",
                                "GGAGAC",
                                "GGAAAC",
                                "CGAAAG",
                            };
    public:
        GoodHairpin4Singlet() : Strategy2D() {
            name_ = "GoodTetraloopSinglet";
            params_ = {91.87045255f, -2.06462586f};
        }

    public:
        float
        score( Feature2DOP const & feature ) override  {
            auto result(params_[0]);

            for(const auto& hairpin : feature->hairpins()) {
                if(hairpin->buffer() == 1 && good_.find(hairpin->sequence) != good_.end()) {
                    result += params_[1];
                }
            }


            return result;
        }
    };
}

#endif // RNAMAKE_GOOD_HAIRPIN4_H
