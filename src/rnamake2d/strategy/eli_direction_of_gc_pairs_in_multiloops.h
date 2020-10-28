#ifndef RNAMAKE_ELIDIRECTIONOFGCPAIRSINMULTILOOPS_H
#define RNAMAKE_ELIDIRECTIONOFGCPAIRSINMULTILOOPS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliDirectionOfGCPairsInMultiloops : public Strategy2D {
    private:

    public:
        EliDirectionOfGCPairsInMultiloops() : Strategy2D() {
            params_ = {-0.5961869866498919, -4.650528887868891};
            name_ = "eli_direction_of_gc_pairs_in_multiloops";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto penalty(0.f);

            auto multiloop_found(false);

            for(const auto& element : feature->elements) {
                if(element.type_ == RNAELEMENT::LOOP) {
                    continue;
                }
                const auto closing_pairs = element.get_loop_closing_pairs(
                        feature->sequence, feature->e_pairmap );

                if(closing_pairs.size() < 3)  {
                    continue;
                }
                multiloop_found = true;

                for(auto index = 0; index < closing_pairs.size(); ++index) {
                    const auto& pair = closing_pairs[index];
                    if(index == 0) {
                       if(pair == "CG") {
                           continue;
                       } else if (pair == "GC") {
                           penalty += params_[0];
                       } else {
                           penalty += params_[1];
                       }
                    } else {
                        if(pair == "CG") {
                            continue;
                        } else if (pair == "GC") {
                            penalty += params_[0];
                        } else {
                            penalty += params_[1];
                        }
                    }
                }
            }

            if(!multiloop_found) {
                return UNSCORABLE;
            }
            return 100 + penalty;
        }
    };

}
#endif // RNAMAKE_ELIDIRECTIONOFGCPAIRSINMULTILOOPS_H 