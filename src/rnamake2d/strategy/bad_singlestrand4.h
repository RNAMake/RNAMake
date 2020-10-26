#ifndef RNAMAKE_BADSINGLESTRAND4_H
#define RNAMAKE_BADSINGLESTRAND4_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class BadSingleStrand4 : public Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                                "GAUA",
                                "UUUA",
                                "UUGG",
                                "AUUC",
                                "GAGG",
                                "GUGU",
                                "UCUU",
                                "UUUC",
                                "CCAA",
                                "AUUA",
                                "GGGA",
                            };
    public:
        BadSingleStrand4() : Strategy2D() {
            params_ = {14.13636364};
            name_ = "BadSingleStrand4";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() == util::MotifType::SSTRAND) {
                   if(bad_.find(motif->sequence) != bad_.end()) {
                       result -= params_[1];
                   }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_BADSINGLESTRAND4_H 