#ifndef RNAMAKE_ELIDIRECTIONOFGCPAIRSINMULTILOOPSNECKAREA_H
#define RNAMAKE_ELIDIRECTIONOFGCPAIRSINMULTILOOPSNECKAREA_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliDirectionOfGCPairsInMultiloopsNeckarea : public Strategy2D {
    private:

    public:
        EliDirectionOfGCPairsInMultiloopsNeckarea() : Strategy2D() {
            params_ = {-5.772660822831785,
                       -4.7011313433002595};
            name_ = "eli_direction_of_gc_pairs_in_multiloops_neckarea";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto penalty(0.f);
            auto mutliloop_found(false);

            for(auto e_index = 0; e_index < feature->elements.size(); ++e_index) {
                const auto& element = feature->elements[e_index];
                if(element.type_ == RNAELEMENT::LOOP) {
                    continue;
                }
                const auto closing_pairs = element.get_loop_closing_pairs(
                        feature->sequence, feature->pairmap);
                if(closing_pairs.size() < 3) {
                    continue;
                }
                mutliloop_found = true;
                for(auto index = 0; index < closing_pairs.size(); ++index) {
                    const auto& pair = closing_pairs[index];
                    if(e_index == 1 && index == 0) {
                        continue;
                    }
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

            if (!mutliloop_found) {
                return UNSCORABLE;
            }
            return 100.f + penalty;
        }
    };

}
#endif // RNAMAKE_ELIDIRECTIONOFGCPAIRSINMULTILOOPSNECKAREA_H 