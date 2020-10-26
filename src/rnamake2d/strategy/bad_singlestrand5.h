#ifndef RNAMAKE_SINGLESTRAND5
#define RNAMAKE_SINGLESTRAND5

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadSingleStrand5 : public Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                "ACAGA",
                "CGAGA",
                "GCAGA",
                "GGAGA"
        };
    public:
        BadSingleStrand5() : Strategy2D() {
            name_ = "BadSingleStrand5";
            params_ = {9.74393939f};
        }

    public:
        float
        score( Feature2DOP const & feature ) override {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() == util::MotifType::SSTRAND
                && bad_.find(motif->sequence) != bad_.end()) {
                    result -= params_[0];
                }
            }

            return result;
        }
    };

}

#endif // RNAMAKE_SINGLESTRAND5
