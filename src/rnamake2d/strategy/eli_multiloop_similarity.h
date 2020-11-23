#ifndef RNAMAKE_ELIMULTILOOPSIMILARITY_H
#define RNAMAKE_ELIMULTILOOPSIMILARITY_H

#include <unordered_set>
#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliMultiloopSimilarity : public Strategy2D {
    private:

    public:
        EliMultiloopSimilarity() : Strategy2D() {
            params_ = {0.441015625, -1.112109375};
            name_ = "eli_multiloop_similarity";
        }

    private:
        bool
        isExist(std::vector<Ints> const & array , Ints const& element) {
            for(auto ii = 0; ii < array.size(); ++ii) {
                if(comp_vectors(array[ii], element)) {
                    return true;
                }
            }
            return false;
        }

    private:
        void
        getChildInfo(RNAElement* elem, Strings& out, Strings& energy) {
            out.push_back(std::to_string(elem->type_) + std::to_string(elem->indices_.size()));
            energy.push_back(std::to_string(elem->score_));
            for(auto ii = 0; ii < elem->children_.size(); ++ii) {
                getChildInfo(elem->children_[ii], out, energy);
            }
        }


    public:
        float
        score(Feature2DOP const & feature) override {
            auto score(100.f);
            const auto& elements = feature->elements;
            const auto root = *elements.rbegin();
            auto infoGrp = std::vector<Strings>{};
            auto energyGrp = std::vector<Strings>{};
            for(auto ii = 0; ii < root.children_.size(); ++ii) {
                auto parent = root.children_[ii]->children_[0];
                const auto loop_groups = parent->get_loop_groups();
                for(auto jj = 0; jj < parent->children_.size(); ++jj) {
                    auto info = Strings{};
                    auto energy = Strings{};
                    getChildInfo(parent->children_[jj] , info, energy);
                    infoGrp.push_back(info);
                    energyGrp.push_back(energy);
                }

                if(loop_groups.size() == 0) {
                    continue;
                }
                for(auto jj = 0; jj < loop_groups.size() - 1; ++jj) {
                    if ( loop_groups[jj].size() != loop_groups[jj+1].size()) {
                        return score;
                    }
                }
            }

            auto identicalGrp = std::vector<Ints>{};
            for(auto ii = 0; ii < infoGrp.size(); ++ii) {
                for(auto jj = 0; jj < infoGrp.size(); ++jj) {
                    if(ii != jj && infoGrp[ii] == infoGrp[jj]) {
                        auto tmp = Ints{};
                        tmp.push_back(std::min(ii,jj));
                        tmp.push_back(std::max(ii,jj));

                        if(isExist(identicalGrp, tmp) == false ) {
                            identicalGrp.push_back(tmp);
                        }
                    }
                }
            }
            auto penalty(0.f);
            for(auto ii = 0; ii < identicalGrp.size(); ++ii) {
                auto e1(0.f);
                for( auto jj = 0;  jj < energyGrp[identicalGrp[ii][0]].size(); ++jj) {
                    e1 = e1 + std::stof(energyGrp[identicalGrp[ii][0]][jj]);
                }
                auto e2(0.f);
                for(auto jj = 0; jj < energyGrp[identicalGrp[ii][1]].size(); ++jj)  {
                    e2 = e2 + std::stof(energyGrp[identicalGrp[ii][1]][jj]);
                }
                if( std::abs(e1 - e2) > TOLERANCE) {
                    auto e = std::abs(e1 - e2);
                    penalty = int(penalty) + int((e - float(int(e*1000) % int(params_[0]*1000))))/1000.f / params_[0];
                }
            }
            score = score + penalty + params_[1];
            return score;
         }
    };

}
#endif // RNAMAKE_ELIMULTILOOPSIMILARITY_H 