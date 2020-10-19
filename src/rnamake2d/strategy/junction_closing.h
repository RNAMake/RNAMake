#ifndef RNAMAKE_JUNCTION_CLOSING_H
#define RNAMAKE_JUNCTION_CLOSING_H

#include <base/types.h>

#include <rnamake2d/strategy2d.h>
namespace rnamake2d {

    class JunctionClosing : public Strategy2D {
         float base_ = 93.60409583;
         std::vector<float> params_;

    public:
            JunctionClosing() : Strategy2D() {
                name_ = "JunctionClosing";
                params_ = {
                        -5.62182494, // pct au
                        -1.67773463, // pct gc
                        -5.66339269, // pct gu
                        -6.45617608  // pct other
                };
            }

    public:
        float
        score( Feature2DOP const & feature ) override {
            float num_ua(0.f);
            float num_gc(0.f);
            float num_ug(0.f);
            float unknown(ook0.f);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::NWAY && motif->mtype() != util::MotifType::TWOWAY) {
                    continue;
                }
                num_ua += float(std::dynamic_pointer_cast<rnamake2d::Junction>(motif)->au);
                num_gc += float(std::dynamic_pointer_cast<rnamake2d::Junction>(motif)->gc);
                num_ug += float(std::dynamic_pointer_cast<rnamake2d::Junction>(motif)->gu);
                unknown += float(std::dynamic_pointer_cast<rnamake2d::Junction>(motif)->unknown);
            }

            const float total = num_ua + num_ug + num_gc + unknown;
            num_gc /= total;
            num_ug /= total;
            num_ua /= total;
            unknown /= total;

            return base_ + num_ua*params_[0] + num_gc*params_[1] + num_ug*params_[2] + unknown*params_[3];
        }
    };

}


#endif // RNAMAKE_JUNCTION_CLOSING_H