#ifndef __BAD_SINGLE_STRAND3_H__
#define __BAD_SINGLE_STRAND3_H__

#include <unordered_set>

#include <base/types.h>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadSingleStrand3 : public Strategy2D {
    private:
        float penalty_;
    private:
        std::unordered_set< String > bad_;
    public:
        BadSingleStrand3() : Strategy2D() {
            bad_ = {
                    "CCU",
                    "ACC",
                    "AAU",
                    "GCA",
                    "GCU",
                    "GUU",
                    "UAC",
                    "GUU",
                    "AAU",
                    "CUG"
            };
            name_ = "BadSingleStrand3";
            penalty_ = 16.4;
        }

    public:
        float
        score(Feature2DOP const& feature) override {
            auto bad_ct(0);
            for(const auto& motif : feature->motifs) {
                if(motif->token() != "SingleStrand3") {
                    continue;
                }
                if(bad_.find(motif->sequence) != bad_.end()) {
                    ++bad_ct;
                }
            }

            return 100.f - bad_ct*penalty_;
        }
    };

}

#endif // __BAD_SINGLE_STRAND3_H__