#ifndef RNAMAKE_GOOD_JUNCTION3_3_3_3_H
#define RNAMAKE_GOOD_JUNCTION3_3_3_3_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
   class GoodJunction3_3_3_3 : public Strategy2D  {
   private:
       std::unordered_set<String> good_ = {
               "UAAAG&CAAAC&GAAAG",
               "UAAAG&CAAAG&CAAAG",
               "UAAAG&CAAAG&CAAAA",
               "UAACG&CAACG&CAAAG",
               "GAAAC&GAAAC&GAAAC",
               "GAAAC&GAAAC&GAAAU",
       };

   public:
       GoodJunction3_3_3_3() : Strategy2D() {
            params_ = {91.82503933,  1.09075014};
            name_ = "GoodJunction3_3_3_3";
       }

   public:
       float
       score( Feature2DOP const & feature) override {
            auto result(params_[0]);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::NWAY)  {
                    continue;
                }
                if(good_.find(motif->sequence) != good_.end()) {
                    result -= params_[1];
                }
            }

            return  result;
       }
   };
}

#endif // RNAMAKE_GOOD_JUNCTION3_3_3_3_H