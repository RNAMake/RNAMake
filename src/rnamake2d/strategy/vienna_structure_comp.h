#ifndef RNAMAKE_VIENNASTRUCTURECOMP_H
#define RNAMAKE_VIENNASTRUCTURECOMP_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>
#include <vienna/vienna.h>

namespace rnamake2d {
    class ViennaStructureComp : public Strategy2D {
    private:
        vienna::Vienna vienna_;
    public:
        ViennaStructureComp() : Strategy2D() {
            params_ = {85.30683151,  6.89077945};
            name_ = "ViennaStructureComp";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(params_[0]);
            vienna_.fold(feature->sequence);

            const auto predicted = vienna_.get_structure();
            auto agree(0);
            for(auto ii = 0; ii < feature->sequence.size(); ++ii) {
                const auto targ = feature->target[ii];
                const auto pred = predicted[ii];

                if(targ == pred && targ == '.') {
                    ++agree;
                } else if (targ != '.' && pred != '.') {
                    ++agree;
                }

            }

            return result + params_[2] *(float(agree) / float(feature->sequence.size()));
        }
    };

}
#endif // RNAMAKE_VIENNASTRUCTURECOMP_H 