#ifndef RNAMAKE_ELIWRONGDIRECTIONOFGCPAIRSINMULTILOOPS_H
#define RNAMAKE_ELIWRONGDIRECTIONOFGCPAIRSINMULTILOOPS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliWrongDirectionOfGCPairsInMultiloops : public Strategy2D {
    private:

    public:
        EliWrongDirectionOfGCPairsInMultiloops() : Strategy2D() {
            params_ = {-2.1793408203125};
            name_ = "eli_wrong_direction_of_gc_pairs_in_multiloops";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto penalty(0.f);

            auto multiloop_found(false);

            for(const auto& element : feature->elements) {
                if(element.type_ != RNAELEMENT::LOOP) {
                    continue;
                }
                const auto closing_pairs = element.get_loop_closing_pairs(
                       feature->sequence, feature->e_pairmap);
                if(closing_pairs.size() < 3) {
                    continue;
                }
                multiloop_found = true;
                for(auto index = 0; index < closing_pairs.size(); ++index) {
                    const auto& pair = closing_pairs[index];
                    if(index == 0) {
                        if(pair == "GC") {
                            continue;
                        } else if (pair == "CG") {
                            ++penalty;
                        } else {
                            penalty += 2;
                        }
                    } else {
                        if(pair == "CG") {
                            continue;
                        } else if (pair == "GC") {
                            ++penalty;
                        } else {
                            penalty += 2;
                        }
                    }
                }
            }

            if(!multiloop_found) {
                return UNSCORABLE;
            }
            return 100.f + params_[0]*penalty;
        }
    };

}
#endif // RNAMAKE_ELIWRONGDIRECTIONOFGCPAIRSINMULTILOOPS_H 