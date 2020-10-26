#ifndef RNAMAKE_BAD_0_1_BULGE
#define RNAMAKE_BAD_0_1_BULGE

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class Bad_0_1_Bulge : public rnamake2d::Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                "GAC&GC",
                "CAU&AG",
                "GCG&CC",
                "GA&UAC",
                "GAG&CC",
                "GG&CGC",
                "CGC&GG",
                "AG&CGU",
                "GU&GUA",
        };
    public:
        Bad_0_1_Bulge() : Strategy2D() {
            name_ = "BadBulge_0_1";
            params_ = {6.61528926f};
        }

    public:
        float
        score(Feature2DOP const &feature) override {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
               if(motif->mtype() == util::MotifType::NWAY &&
                    bad_.find(motif->sequence) != bad_.end()) {
                    result -= params_[0];
               }
            }

            return result;
        }
    };
}



#endif // RNAMAKE_BAD_0_1_BULGE
