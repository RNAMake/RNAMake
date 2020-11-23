#ifndef RNAMAKE_MERRYSKIESONLYASINTHELOOPS_H
#define RNAMAKE_MERRYSKIESONLYASINTHELOOPS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class MerrySkiesOnlyAsInTheLoops : public Strategy2D {
    private:

    public:
        MerrySkiesOnlyAsInTheLoops() : Strategy2D() {
            params_ = {};
            name_ = "merryskies_only_as_in_the_loops";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            float result(100.f);
            //const auto seq = feature->sequence;
            //const auto target = feature->target;

            //for(auto ii = 0 ; ii<seq.size(); ++ii){
            //    if( target[ii] == '.' && seq[ii] != 'A') {
            //        result -= 1.f;
            //    }
            //}

            return result;
        }
    };

}
#endif // RNAMAKE_MERRYSKIESONLYASINTHELOOPS_H 