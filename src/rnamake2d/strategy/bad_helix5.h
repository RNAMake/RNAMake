#ifndef RNAMAKE_BADHELIX5_H
#define RNAMAKE_BADHELIX5_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class BadHelix5 : public Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                "GUCCU&AGGAC",
                "CUAUG&CAUAG"
        };
    public:
        BadHelix5() : Strategy2D() {
            params_ = {1.};
            name_ = "BadHelix5";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(const auto& helix : feature->motifs) {
                if(helix->mtype() != util::MotifType::HELIX || helix->sequence.size() != 5) {
                    continue;
                }

                if (bad_.find(helix->sequence) != bad_.end()) {
                    result -= params_[0];
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_BADHELIX5_H 