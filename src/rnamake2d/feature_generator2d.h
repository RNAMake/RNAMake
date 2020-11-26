#ifndef __FEATURE_GENERATOR2D_H__
#define __FEATURE_GENERATOR2D_H__

#include <eternabot/feature_generator.h>
#include <rnamake2d/ss_motif.h>
#include <rnamake2d/rna_element.h>

namespace rnamake2d {

    struct Feature2D : eternabot::Features {

        Motif2DOPs motifs;
        String target;
        String sequence;
        RNAElems elements;
        Ints e_pairmap;
    private:
        JunctionOPs junctions_;
        HairpinOPs hairpins_;
        SingleStrandOPs singlestrands_;
        HelixOPs helices_;
    public:
        explicit Feature2D(eternabot::Features &features) {
            this->gu = features.gu;
            this->gc = features.gc;
            this->ua = features.ua;
            this->meltpoint = features.meltpoint;
            this->fe = features.fe;
            this->a_count = features.a_count;
            this->c_count = features.c_count;
            this->g_count = features.g_count;
            this->u_count = features.u_count;
            this->pairmap = std::move(features.pairmap);
            this->length = features.length;
        }

        Feature2D() = default;

        Feature2D &
        operator=(eternabot::Features &features) {
            this->gu = features.gu;
            this->gc = features.gc;
            this->ua = features.ua;
            this->meltpoint = features.meltpoint;
            this->fe = features.fe;
            this->a_count = features.a_count;
            this->c_count = features.c_count;
            this->g_count = features.g_count;
            this->u_count = features.u_count;
            this->pairmap = std::move(features.pairmap);
            this->length = features.length;
            return *this;
        }

    public:
        JunctionOPs const &
        junctions() const {
            return junctions_;
        }

    public:
        void
        reset_motifs() {
            if(!motifs.empty() ){
                motifs.clear();
            }

            if(!helices_.empty()) {
                helices_.clear();
            }

            if(!junctions_.empty()) {
                junctions_.clear();
            }

            if(!hairpins_.empty()) {
                hairpins_.clear();
            }

            if(!singlestrands_.empty()) {
                singlestrands_.clear();
            }
        }

        void
        junctions(JunctionOPs& junctions) {
            for(auto& junc : junctions) {
                junctions_.push_back(junc);
                motifs.push_back(junc);
            }
        }


    public:
        HairpinOPs const &
        hairpins() const {
            return hairpins_;
        }

        void
        hairpins(HairpinOPs& hairpins) {
            for(auto& hp : hairpins) {
                hairpins_.push_back(hp);
                motifs.push_back(hp);
            }
        }

    public:
        SingleStrandOPs const &
        singlestrands() const {
            return singlestrands_;
        }

        void
        singlestrands(SingleStrandOPs& singlestrands) {
            for(auto& ss : singlestrands) {
                singlestrands_.push_back(ss);
                motifs.push_back(ss);
            }
        }

    public:
        HelixOPs const &
        helices() const {
            return helices_;
        }

        void
        helices(HelixOPs& helices) {
            for(auto& hel : helices) {
                helices_.push_back(hel);
                motifs.push_back(hel);
            }
        }
    };
} // namespace rnamake2d

using Feature2DOP = std::shared_ptr<rnamake2d::Feature2D>;
using Feature2DOPs = std::vector<Feature2DOP>;

namespace rnamake2d {

    class FeatureGenerator2d : public eternabot::FeatureGenerator {
    public:
        FeatureGenerator2d() : eternabot::FeatureGenerator() {

        }

    public:
        void
        update_features(
                Feature2DOP &,
                secondary_structure::PoseOP const &);
    };
}


#endif// __FEATURE_GENERATOR2D_H__