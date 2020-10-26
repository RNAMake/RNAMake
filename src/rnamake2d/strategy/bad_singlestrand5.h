#ifndef RNAMAKE_SINGLESTRAND5
#define RNAMAKE_SINGLESTRAND5

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadSingleStrand5 : Strategy2D {
    private:
    public:
        BadSingleStrand5() : Strategy2D() {

        }

    public:
        float
        score( Feature2DOP const & feature ) override {

        }
    };

}

#endif // RNAMAKE_SINGLESTRAND5
