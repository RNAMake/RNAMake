#ifndef __FEATURE_GENERATOR2D_H__
#define __FEATURE_GENERATOR2D_H__

#include <eternabot/feature_generator.h>
#include <rnamake2d/SSMotif.h>

namespace rnamake2d {

    struct Feature2D : eternabot::Features {
        Motif2DOPs motifs;

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
