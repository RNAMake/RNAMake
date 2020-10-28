#ifndef RNAMAKE_VIENNAMFENORMALIZED_H
#define RNAMAKE_VIENNAMFENORMALIZED_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>
#include <vienna/vienna.h>

namespace rnamake2d {
    class ViennaMFENormalized : public Strategy2D {
    private:
        vienna::Vienna vienna_;

    public:
        ViennaMFENormalized() : Strategy2D() {
            params_ = {-28.83431937};
            name_ = "ViennaMFENormalized";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            vienna_.fold(feature->sequence);
            const auto mfe_normalized = vienna_.free_energy() / float(feature->sequence.size());

            return 100.f - mfe_normalized*params_[0];

        }
    };

}
#endif // RNAMAKE_VIENNAMFENORMALIZED_H 