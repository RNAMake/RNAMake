

#include "../common.hpp"

#include "util/steric_lookup.hpp"
#include "util/random_number_generator.h"

class PointGenerator {
public:
    PointGenerator():
    test_p_(math::Vector3()),
    rng_(util::RandomNumberGenerator())
    {}

    ~PointGenerator() {}

public:
  inline math::Vector3 const & rand_point(int scale) {
    test_p_.set_x(rng_.rand() * rng_.randrange(scale));
    test_p_.set_y(rng_.rand() * rng_.randrange(scale));
    test_p_.set_z(rng_.rand() * rng_.randrange(scale));

    if(rng_.randrange(1000) < 500) { test_p_.x(-test_p_.get_x()); }
    if(rng_.randrange(1000) < 500) { test_p_.y(-test_p_.get_y()); }
    if(rng_.randrange(1000) < 500) { test_p_.z(-test_p_.get_z()); }

    return test_p_;
  }


private:
    math::Vector3 test_p_;
    util::RandomNumberGenerator rng_;
};


TEST_CASE( "Test Steric Lookup for quick Sterics " ) {
    
    auto p_generator = PointGenerator();

    SUBCASE("Test adding points to lookup") {
        auto sl = util::StericLookupNew();
        auto p = math::Vector3();
        sl.add_point(p);
        sl.to_pdb("grid.pdb");
        
        auto clash = 0, sl_clash = 0;
        double dist = 0;
        auto test_p = math::Vector3();
        int count = 0;
        for(int i = 0; i < 1000; i++) {
            clash = 0;
            test_p = p_generator.rand_point(10);
            dist = p.distance(test_p);
            
            if(dist < 2.65) { clash = 1; }
            sl_clash = sl.clash(test_p);
            
            if(sl_clash == 0 && clash == 1) {
                count += 1;
            }
            
        }

        CHECK(count < 100);
    }
    
    SUBCASE("Test adding a set of points to lookup") {
        auto points = math::Vector3s();
        for(int i = 0; i < 1000; i++) {
            points.push_back(p_generator.rand_point(100));
        }
        
        auto sl = util::StericLookup();
        sl.add_points(points);
        
        auto dist = 0.0f;
        auto test_p = math::Vector3();
        auto clash = 0, miss_count = 0, sl_clash = 0;
        
        for(int i = 0; i < 10000; i++) {
            clash = 0;
            test_p = p_generator.rand_point(100);
            
            for(auto const & p : points) {
                dist = test_p.distance(p);
                if(dist < 2.65) {
                    clash = 1;
                    break;
                }
            }
            
            sl_clash = sl.clash(test_p);
            
            if(sl_clash != clash) { miss_count += 1; }
        }

        CHECK(miss_count < 200);
    }
    
}
