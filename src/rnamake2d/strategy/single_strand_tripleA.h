#ifndef RNAMAKE_SINGLE_STRAND_TRIPLEA
#define RNAMAKE_SINGLE_STRAND_TRIPLEA

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class SingleStrandTripleA : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                        "AAA",
                        "AAC",
                        "AUA",
                        "GAA"
        };

    public:
        SingleStrandTripleA() : Strategy2D() {
            name_ = "SingleStrandTripleA";
            params_ = {90.5849651f,  0.1915561f };
        }

    public:
        float
        score( Feature2DOP const & feature  ) override {
           auto result = params_[0];

           for(const auto& motif : feature->motifs) {
               if(motif->mtype() == util::MotifType::SSTRAND) {
                   result += params_[1];
               }
           }

           return result;
        }
    };
}

#endif //  RNAMAKE_SINGLE_STRAND_TRIPLEA