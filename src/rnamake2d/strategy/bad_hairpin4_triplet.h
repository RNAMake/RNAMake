#ifndef RNAMAKE_BAD_HAIRPIN4_TRIPLET
#define RNAMAKE_BAD_HAIRPIN4_TRIPLET

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadHaripin4Triplet : public Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                                "UGAGAG",
                                "AGAAAU",
                                "UGAAAA",
                                "GAACUC",
                                "CGCAAG",
                                "GUAAUC",
                                "GUAAGC",
                                "AUAUUU",
                                "GGCAAC",
                                "AAUCCA",
                                "CAAAAG",
                                "GAAAAC",
                                "UGAAAG",
                                "AUAUAU",
                                "GAUCAC",
                                "AAACCU",
                                "CUUCGG",
                                "AUAAUU",
                                "CUCUUU",
                                "AAUAAU",
                                "CACCAG",
                                "UUUCGA",
                                "CUCCGG",
                                "GGUGAC",};

    public:
        BadHaripin4Triplet() : Strategy2D() {
            name_ = "BadHairpin4Triplet";
            params_ =  { 6.67566766f };
        }

    public:
        float
        score(Feature2DOP const & feature ) override {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() == util::MotifType::HAIRPIN
                    && motif->buffer_ == 3
                    && bad_.find(motif->sequence) != bad_.end()
                    ) {
                   result -= params_[0];
                }
            }

            return result;
        }
    };

}

#endif // RNAMAKE_BAD_HAIRPIN4_TRIPLET