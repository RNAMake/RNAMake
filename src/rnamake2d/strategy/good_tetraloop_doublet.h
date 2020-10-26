#ifndef RNAMAKE_TETRALOOP_DOUBLET_H
#define RNAMAKE_TETRALOOP_DOUBLET_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class GoodTetraloopDoublet : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                "CAGGCC",
                "CUAACG",
                "CUAGUG",
                "CGAGAG",
                "GGAAAC",
                "CGAAAG",
        };
    public:
        GoodTetraloopDoublet() : Strategy2D() {
            params_ = {91.86279615, -0.58356472};
            name_ = "GoodTetraloopDoublet";
        }

    public:
        float
        score( Feature2DOP const & feature ) override {
            auto result(params_[0]);

            for(const auto& motif : feature->motifs) {
               if(motif->mtype() != util::MotifType::HAIRPIN) {
                   continue;
               }

               if(good_.find(motif->sequence) != good_.end()) {
                    result += params_[1];
               }
            }

            return result;
        }
    };
}

#endif // RNAMAKE_TETRALOOP_DOUBLET_H
