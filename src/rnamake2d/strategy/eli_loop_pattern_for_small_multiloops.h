#ifndef RNAMAKE_ELILOOPPATTERNFORSMALLMULTILOOPS_H
#define RNAMAKE_ELILOOPPATTERNFORSMALLMULTILOOPS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliLoopPatternForSmallMultiloops : public Strategy2D {
    private:
        std::unordered_set<int> skippable_indices_;
        Strings good_patterns = {
                "CGGCGCCG",
                "GCCGCGGC",
                "CGCGCGCG",
        };
    public:
        EliLoopPatternForSmallMultiloops() : Strategy2D() {
            params_ = {100., -1., -2.};
            name_ = "eli_loop_pattern_for_small_multiloops";
        }

    private:
        static
        bool
        had_non_gc_(String const& str) {
            if(str.find('A') != String::npos) {
                return true;
            } else if (str.find('U') != String::npos) {
                return true;
            }
            return false;
        }

    private:
        static
        float
        compare_pairs_(String const& string1, String const& string2, std::vector<float> const& params) {
            auto errors(0.f);

            for(auto ii = 0; ii <= 6; ii += 2) {
                if(string1.substr(ii,2) != string2.substr(ii,2)) {
                    if(had_non_gc_(string2.substr(ii, 2))) {
                        errors += params[2];
                    } else {
                        errors += params[1];
                    }
                }
            }
            return errors;
        }

    private:
        String
        get_gc_multiloop_members_(int index, Feature2DOP const& design) {
            auto returnstring = String("");
            const auto& sequence = design->sequence;
            const auto& pairmap = design->e_pairmap;
            const auto paired_index = pairmap[index];
            returnstring += sequence[index];
            if(paired_index != -1 ) {
                returnstring += sequence[paired_index];
                const auto index2 = paired_index - 1;
                if(skippable_indices_.find(index2) != skippable_indices_.end()) {
                    return "";
                }
                returnstring += sequence[index2];
                const auto paired_index2 = pairmap[index2];
                if(paired_index2 != -1) {
                    returnstring += sequence[paired_index2];
                    const auto index3 = paired_index2 - 1;
                    if(skippable_indices_.find(index3) != skippable_indices_.end() or index3 == index) {
                        return "";
                    }
                    returnstring += sequence[index3];
                    const auto paired_index3 = pairmap[index3];
                    if(paired_index3 != -1 ) {
                        returnstring += sequence[paired_index3];
                        const auto index4 = paired_index3 - 1;
                        if(skippable_indices_.find(index4) != skippable_indices_.end() or index4 == index) {
                            return "";
                        }
                        returnstring += sequence[index4];
                        const auto paired_index4 = pairmap[index4];
                        if(paired_index4 != -1 ) {
                            returnstring += sequence[paired_index4];
                            return returnstring;
                        }
                    }
                }
            }
            return "";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            const auto& pairmap = feature->e_pairmap;

            auto score = params_[0];
            auto gc_multiloop_found(false);
            if(!skippable_indices_.empty())  {
                skippable_indices_.clear();
            }

            for(auto ii = 0; ii < pairmap.size(); ++ii) {
                if(skippable_indices_.find(ii) != skippable_indices_.end()) {
                    continue;
                }
                const auto gc_multiloop_string = get_gc_multiloop_members_(ii, feature);

                if(gc_multiloop_string != "") {
                    if(!gc_multiloop_found) {
                        gc_multiloop_found = true;
                    }

                    skippable_indices_.insert(ii);
                    skippable_indices_.insert(pairmap[ii]);

                    if (std::find(good_patterns.begin(), good_patterns.end() ,gc_multiloop_string) != good_patterns.end()) {
                            continue;
                        }
                    } else {
                        auto errors = compare_pairs_( good_patterns[0], gc_multiloop_string, params_);
                        errors = std::max(errors, compare_pairs_( good_patterns[1], gc_multiloop_string, params_));
                        errors = std::max(errors, compare_pairs_( good_patterns[2], gc_multiloop_string, params_));
                        score += errors;
                    }
                }
            if(!gc_multiloop_found) {
                return UNSCORABLE;
            }
            return score;
        }
    };

}
#endif // RNAMAKE_ELILOOPPATTERNFORSMALLMULTILOOPS_H 