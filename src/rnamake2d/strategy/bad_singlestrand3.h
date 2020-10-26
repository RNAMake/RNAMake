#ifndef RNAMAKE_BAD_SINGLESTRAND3
#define RNAMAKE_BAD_SINGLESTRAND3

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadSingleStrand3 : public rnamake2d::Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                    "CCU",
                    "ACC",
                    "AAU",
                    "GCA",
                    "GCU",
                    "GUU",
                    "UAC",
                    "GUU",
                    "AAU",
                    "CUG",
        };

    public:
        BadSingleStrand3() : Strategy2D() {
            name_ = "BadSingleStrand3";
            params_ = {16.4};
        }

    public:
        float
        score(Feature2DOP const & feature) const {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
               if(motif->mtype() == util::MotifType::SSTRAND) {
                   if(bad_.find(motif->sequence) != bad_.end()) {
                       result -= params_[0];
                   }
               }
            }

            return result;
        }
    };
}

#endif // RNAMAKE_BAD_SINGLESTRAND3
