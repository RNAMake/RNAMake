#include <rnamake2d/feature_generator2d.h>



namespace rnamake2d {
    void
    FeatureGenerator2d::update_features(
             Feature2DOP   & feature,
            secondary_structure::PoseOP const & pose ) {
            // original eternabot stuff
            this->eternabot::FeatureGenerator::update_features(&(*feature), pose);
    }
} // namespace rnamake2d
