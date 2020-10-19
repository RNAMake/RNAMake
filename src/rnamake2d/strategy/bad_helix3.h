#ifndef __BAD_HELIX3_H__
#define __BAD_HELIX3_H__

#include <base/types.h>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadHelix3 : public Strategy2D  {
        float penalty_ = 8.26555556;
        std::unordered_set<String> bad_ = {
                "ACU&AGU",
                "AGA&UCU",
                "AGU&AUU",
                "CAA&UUG",
                "GCG&UGC",
                "GGC&GUC",
                "GGG&GGG",
                "GUA&UAC",
        };

    public:
        BadHelix3() : Strategy2D()  {
            name_ = "BadHelix3";
        }

    public:
        ~BadHelix3() = default;

    public:
        float
        score(Feature2DOP const & features ) override  {
            auto bad_helix_ct(0);

            for(const auto& motif : features->motifs) {
                if(motif->token() != "Helix3") {
                    continue;
                }
                if(bad_.find(motif->sequence) != bad_.end())  {
                    ++bad_helix_ct;
                }
            }

            return 100.f - bad_helix_ct*penalty_;
        }
    };

}


#endif // __BAD_HELIX3_H__
