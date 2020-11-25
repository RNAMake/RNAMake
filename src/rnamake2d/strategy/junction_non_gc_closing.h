#ifndef RNAMAKE_NON_GC_CLOSING_H
#define RNAMAKE_NON_GC_CLOSING_H

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class JunctionNonGCClosing : public Strategy2D {
    public:
        JunctionNonGCClosing() : Strategy2D() {
            name_ = "JunctionNonGCClosing";
            params_ = {17.29684124};
        }

    public:
        float
        score(Feature2DOP const& feature) override {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::NWAY && motif->mtype() != util::MotifType::TWOWAY) {
                    continue;
                }
                std::cout<<motif<<std::endl;
                //const auto& jnc = std::dynamic_pointer_cast<rnamake2d::Junction>(motif);
                const auto jnc = (rnamake2d::Junction*)motif.get();
                result -= params_[0]*(jnc->au + jnc->gu + jnc->unknown);
            }

            return result;
        }
    };
}

#endif // RNAMAKE_NON_GC_CLOSING_H