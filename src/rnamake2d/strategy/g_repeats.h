#ifndef RNAMAKE_G_REPEATS_H
#define RNAMAKE_G_REPEATS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class GRepeats : public Strategy2D  {

    public:
        GRepeats() : Strategy2D() {
            name_ = "GRepeats";
            params_ =  { 5.17491348f };
        }

    public:
        float
        score( Feature2DOP const & feature ) override {
            auto result(100.f);
            // NEEDS the eternabot

            return result;
        }
    };
}

#endif // RNAMAKE_G_REPEATS_H