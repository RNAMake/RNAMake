#ifndef RNAMAKE_HELIX_4_H
#define RNAMAKE_HELIX_4_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class Helix4 : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                                           "CAGC&GCUG",
                                           "CUAC&GUAG",
                                           "GCAG&CUGC",
                                           "GCGC&GCGC",
                                           "GCUC&GAGC",
                                           "GCUG&CAGC",
                                           "GGAC&GUCC",
                                           "GGCC&GGCC",
                                           "UCAC&GUGA",
                                           "UCUG&CAGA",
                                           "UGAC&GUCA",
                                        };
    public:
        Helix4() : Strategy2D() {
            params_ = { 92.35644786,  0.33235129,  1.04569803 };
            name_ =  "Helix4";
        }

    public:
        float
        score(Feature2DOP const & feature ) override {
            auto result(params_[0]);
            for(const auto& motif : feature->motifs) {
                if(motif->sequence.size() != 9)  {
                    continue;
                }

                if(good_.find(motif->sequence) != good_.end()) {
                    result += params_[1];
                } else {
                    result -= params_[2];
                }

            }

            return result;
        }
    };

}

#endif// RNAMAKE_HELIX_4_H
