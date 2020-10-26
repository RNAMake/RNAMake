#ifndef RNAMAKE_SINGLESTRAND5_H
#define RNAMAKE_SINGLESTRAND5_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class SingleStrand5 : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                            "AAAAA",
                            "AAACG",
                            "AAAUA",
                            "AACAA",
                            "AAGAA",
                            "AAUAA",
                            "GAAAA",
                            "GGAAC",
                            "GGAAG",
                        };

    public:
        SingleStrand5() : Strategy2D() {
            params_ = {7.41142174};
            name_ = "SingleStrand5";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);
            for(const auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::SSTRAND ||
                    motif->sequence.size() != 5
                ) {
                    continue;
                }

                if(good_.find(motif->sequence) == good_.end()) {
                    result -= params_[0];
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_SINGLESTRAND5_H 