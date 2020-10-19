#ifndef __STRATEGY2D_H__
#define __STRATEGY2D_H__

#include <memory>

#include <base/types.h>
#include <eternabot/strategy.h>

#include <rnamake2d/feature_generator2d.h>

namespace rnamake2d {

    class Strategy2D : public eternabot::Strategy {
        protected:
            String name_;
        private:
            float
            score(const eternabot::FeaturesOP & features) final {
                return 0.f;
            }

        public:
            virtual
            float
            score(const Feature2DOP& ) = 0;

        public:
            String
            name() const {
                return name_;
            }
    };


}

using Strategy2DOP = std::shared_ptr<rnamake2d::Strategy2D>;

namespace rnamake2d {
    Strategy2DOP
    get_strategy(String const&);
}


#endif// __STRATEGY2D_H__