//
// Created by Joseph Yesselman on 3/15/19.
//

#ifndef RNAMAKE_NEW_HASH_MOTIFS_H
#define RNAMAKE_NEW_HASH_MOTIFS_H

#include <unordered_map>

#include <stdio.h>
#include "base/application.hpp"
#include "math/hashing.h"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_state_graph.hpp"


class MotifPathHash {
public:
    MotifPathHash(
            math::BoundingBox const & bounding_box,
            math::Real6 const & bin_widths):
            binner_(math::SixDCoordinateBinner(bounding_box, bin_widths)) {
        additions_ = std::vector<math::Real6>();
        dummy_ = math::Real6();
        _setup_additions();
    }

    MotifPathHash(
            std::ifstream & in) :
            binner_(math::BoundingBox(), math::Real6{0.1, 0.1, 0.1, 0.1, 0.1, 0.1}) {
        additions_ = std::vector<math::Real6>();
        dummy_ = math::Real6();
        _setup_additions();

        auto datas = math::Real3();
        for (int i = 0; i < 3; i++) {
            in.read(reinterpret_cast<char *>(&datas[i]), sizeof(datas[i]));
        }
        auto lower = math::Point(datas[0], datas[1], datas[2]);
        for (int i = 0; i < 3; i++) {
            in.read(reinterpret_cast<char *>(&datas[i]), sizeof(datas[i]));
        }
        auto upper = math::Point(datas[0], datas[1], datas[2]);
        auto bin_widths = math::Real6();
        for (int i = 0; i < 6; i++) {
            in.read(reinterpret_cast<char *>(&bin_widths[i]), sizeof(bin_widths[i]));
        }
        binner_ = math::SixDCoordinateBinner(math::BoundingBox(lower, upper), bin_widths);
        u_int64_t num;
        u_int64_t key;
        in.read(reinterpret_cast<char *>(&num), sizeof(num));
        for (int i = 0; i < num; i++) {
            unsigned size;
            in.read(reinterpret_cast<char *>(&key), sizeof(key));
            in.read(reinterpret_cast<char *>(&size), sizeof(size));
            std::vector<int> values(size);
            in.read(reinterpret_cast<char *>(&values[0]), size*sizeof(int));
            stored_[key] = values;
        }
    }

private:
    void
    _setup_additions() {
        float grid_size = 0.5;
        float cutoff = 3;
        float radius = 10;

        auto add = Floats();
        for (int i = 1; i < radius; i++) {
            add.push_back(float(-i * grid_size));
        }
        add.push_back(0);
        for (int i = 1; i < radius; i++) {
            add.push_back(float(i * grid_size));
        }

        auto deg_add = Floats();
        for (int i = 1; i < 2; i++) {
            deg_add.push_back(float(-i * 10));
        }
        deg_add.push_back(0);
        for (int i = 1; i < 2; i++) {
            deg_add.push_back(float(i * 10));
        }

        float dist = 0;
        math::Point p, rad_p;
        math::Point origin(0, 0, 0);
        math::Real6 values;
        for (auto const & x : add) {
            for (auto const & y : add) {
                for (auto const & z : add) {
                    p = math::Point(x, y, z);
                    dist = p.distance(origin);
                    if (dist > cutoff) { continue; }

                    for(auto const & a : deg_add) {
                        for(auto const & b : deg_add) {
                            for(auto const & g : deg_add) {
                                rad_p = math::Point(a, b, g);
                                dist = p.distance(origin);
                                if(dist > 15) { continue; }

                                values[0] = p[0]; values[1] = p[1]; values[2] = p[2];
                                values[3] = a;    values[4] = b;    values[5] = g;
                                additions_.push_back(values);
                            }
                        }
                    }
                }
            }
        }
    }


public:
    std::vector<Ints>
    search(
            math::Real6 const & values) {
        auto all_vals = std::vector<Ints>();
        auto bin_index = u_int64_t();
        for(auto const & add : additions_) {
            for (int i = 0; i < 6; i++) {
                dummy_[i] = values[i] + add[i];
            }
            bin_index = binner_.bin_index(dummy_);
            if(stored_.find(bin_index) != stored_.end()) {
                all_vals.push_back(stored_[bin_index]);
            }
        }
        return all_vals;
    }


public:
    void
    add(
            math::Real6 const & values,
            Ints const & motifs) {
        auto bin_index = binner_.bin_index(values);
        if (stored_.find(bin_index) == stored_.end()) {
            stored_[bin_index] = std::vector<int>();
        }
        for(auto const & i : motifs) {
            stored_[bin_index].push_back(i);
        }
    }

    void
    add_value(
            math::Real6 const & values) {
        /*for(auto const & add : additions_) {
            for(int i = 0; i < 6; i++) {
                dummy_[i] = values[i] + add[i];
            }
            MotifPathHash::add(dummy_);
        }*/


    }


    void
    output_binary(
            std::ofstream & out) {
        auto lower = binner_.get_bounding_box().lower();
        auto upper = binner_.get_bounding_box().upper();
        out.write((const char *) &lower.x(), sizeof(lower.x()));
        out.write((const char *) &lower.y(), sizeof(lower.y()));
        out.write((const char *) &lower.z(), sizeof(lower.z()));
        out.write((const char *) &upper.x(), sizeof(upper.x()));
        out.write((const char *) &upper.y(), sizeof(upper.y()));
        out.write((const char *) &upper.z(), sizeof(upper.z()));
        for (int i = 0; i < 6; i++) {
            out.write((const char *) &binner_.get_bin_widths()[i], sizeof(binner_.get_bin_widths()[i]));
        }
        u_int64_t num = 0;
        for (auto const & kv : stored_) {
            num += 1;
        }
        unsigned size;
        out.write((const char *) &num, sizeof(num));
        for (auto const & kv : stored_) {
            out.write((const char *) &kv.first, sizeof(kv.first));
            size = kv.second.size();
            out.write((const char *) &size, sizeof(unsigned));
            out.write((const char *) &kv.second[0], kv.second.size()*sizeof(int));
        }

    }

private:
    std::vector<math::Real6> additions_;
    math::Real6 dummy_;
    math::SixDCoordinateBinner binner_;
    std::unordered_map<u_int64_t, std::vector<int>> stored_;

};


class HashMotifs : base::Application {
public:
    HashMotifs();

    ~HashMotifs() {}

public: // application functions

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

public:
    void
    _generate_pairwise_map();

    void
    _old_buildup();

};


#endif //RNAMAKE_NEW_HASH_MOTIFS_H
