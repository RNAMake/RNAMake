#ifndef RNAMAKE_DEVIADDEVIADSTRATEGY_H
#define RNAMAKE_DEVIADDEVIADSTRATEGY_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class DeviadDeviadStrategy : public Strategy2D {
    private:

    public:
        DeviadDeviadStrategy() : Strategy2D() {
            params_ = {-0.5038620441323505,
                       1.1420981943541504,
                       18.593301693968414,
                       -1.044371808090831,
                       1.019593202686596,
                       0.7568781945870047,
                       2.6152963254215056,
                       -4.9203291617294695,
                       -0.13899354254374402,
                       -32.061818182312564,
                       2.570969770531801,
                       -1.3614060249420725,
                       9.876280826765};

            name_ = "deviad_deviad_strategy";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto score(80.f);
            const auto total_pairs = feature->gc + feature->gu + feature->ua;
            const auto gc_rate = feature->gc / total_pairs;
            const auto gu_rate = feature->gu / total_pairs;

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                if(feature->elements[ii].type_ == RNAELEMENT::LOOP) {
                   if(feature->elements[ii].score_ < params_[0]) {
                        score += params_[1];
                   } else if (feature->elements[ii].score_ > params_[2]) {
                       score += params_[3];
                   } else if (feature->elements[ii].score_ > params_[6]) {
                       auto closing_pairs = feature->elements[ii].get_loop_closing_pairs(
                                    feature->sequence, feature->pairmap
                               );
                       auto is_there_gu(false);

                       for(auto jj = 0; jj < closing_pairs.size(); ++jj) {
                           const auto pair = closing_pairs[jj];
                           if(pair == "GU" || pair == "UG") {
                               is_there_gu = true;
                               break;
                           }
                       }

                       if(is_there_gu) {
                           score += params_[7];
                       }
                   }
                } else {
                    const auto stacklen = feature->elements[ii].indices_.size() /2;
                    auto gc_count(0);
                    auto gu_count(0);

                    if(stacklen <= 5 ) {
                        for(auto jj = 0; jj < feature->elements[ii].indices_.size(); jj += 2 ) {
                            const String pair = String{feature->sequence[feature->elements[ii].indices_[jj]]} +
                                    String{feature->sequence[feature->elements[ii].indices_[jj + 1]]};

                            if(pair == "GC" || pair == "CG")  {
                                ++gc_count;
                            }
                            if (pair == "GU" || pair == "UG") {
                                ++gu_count;
                            }
                        }

                        if(gc_count > 1) {
                            score += params_[4];
                        }
                        if(gu_count > 0) {
                            score += params_[5];
                        }

                    }
                }
            }

            for(auto ii = 0; ii < feature->sequence.size(); ++ii) {
                if(feature->sequence[ii] == 'C' && feature->pairmap[ii] < 0) {
                    score += params_[8];
                }
            }

            score += feature->fe / params_[9];

            if(gc_rate > 0.5) {
                score += params_[10];
            }

            if(gu_rate > 0.15) {
                score += params_[11];
            }

            score += feature->meltpoint / params_[12];

            return score;
        }
    };

}
#endif // RNAMAKE_DEVIADDEVIADSTRATEGY_H 