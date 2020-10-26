#ifndef RNAMAKE_BAD_JUNCTION3_3_3_3_H
#define RNAMAKE_BAD_JUNCTION3_3_3_3_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadJunction3_3_3_3 : public Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                "GAACG&CAACG&CAAAU",
                "CAAAG&CAAAG&CAAAG",
                "UACAC&GACAC&GACAG",
                "GAAAG&UAAAG&UAAAU",
        };
    public:
        BadJunction3_3_3_3() : Strategy2D() {
            name_ = "BadJunction3_3_3_3";
            params_ = { 7.06296296 };
        }

    public:
        float
        score(Feature2DOP const & feature ) override {
            auto result(100.f);
            for(const auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::NWAY) {
                    continue;
                }
                if(bad_.find(motif->sequence) != bad_.end())  {
                    result -= params_[0];
                }
            }

            return result;
        }
    };
}

#endif // RNAMAKE_BAD_JUNCTION3_3_3_3_H