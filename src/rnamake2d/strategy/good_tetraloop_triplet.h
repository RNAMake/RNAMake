#ifndef RNAMAKE_GOODTETRALOOPTRIPLET_H
#define RNAMAKE_GOODTETRALOOPTRIPLET_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class GoodTetraloopTriplet : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {"CGGCAG",
                                            "CAUACG",
                                            "GGAAUC",
                                            "GAAACC",
                                            "GGAUAC",
                                            "GGAAAC",
                                            "CGAAGG",
                                            "CGAUAG",
                                            "GGUAAC",
                                            "CGACAG",
                                            "CGAAAG",
                                            "CGAGAG",
                                            "CGGAAG",
                                            "AGAACU",
                                            "GAGAAC",
                                            "CGUGAG",
                                            "CUUAUG",
                                            "CGCCAG",
                                            "CGACGG",
                                            "CAUAAG",
                                            };
    public:
        GoodTetraloopTriplet() : Strategy2D() {
            params_ = {91.90750453, -0.34001696};
            name_ = "GoodTetraloopTriplet";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(params_[0]);

            for(const auto& motif : feature->motifs) {
               if(motif->mtype() != util::MotifType::HAIRPIN ||
                    motif->buffer_ != 3
               ) {
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
#endif // RNAMAKE_GOODTETRALOOPTRIPLET_H 