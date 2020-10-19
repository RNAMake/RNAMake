//
//  feature_generator.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/1/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "eternabot/feature_generator.h"


namespace eternabot {

    FeaturesOP
    FeatureGenerator::get_features(
            secondary_structure::PoseOP const &p) {

        auto features = std::make_shared<Features>();
        for (auto const &bp : p->basepairs()) {
            features->pairmap[bp->res1()->num()] = bp->res2()->num();
            features->pairmap[bp->res2()->num()] = bp->res1()->num();
        }

        features->length = (int) p->residues().size();
        update_features(features, p);
        return features;

    }


    void
    FeatureGenerator::update_features(
            FeaturesOP &features,
            secondary_structure::PoseOP const &p) {

        features->gc = 0;
        features->ua = 0;
        features->gu = 0;
        features->a_count = 0;
        features->c_count = 0;
        features->g_count = 0;
        features->u_count = 0;

        for (auto const &bp : p->basepairs()) {
            if (secondary_structure::is_gc_pair(bp)) { features->gc += 1; }
            else if (secondary_structure::is_au_pair(bp)) { features->ua += 1; }
            else if (secondary_structure::is_gu_pair(bp)) { features->gu += 1; }
        }

        for (auto const &r : p->residues()) {
            if (r->res_type() == secondary_structure::ResType::A) { features->a_count += 1; }
            else if (r->res_type() == secondary_structure::ResType::C) { features->c_count += 1; }
            else if (r->res_type() == secondary_structure::ResType::G) { features->g_count += 1; }
            else if (r->res_type() == secondary_structure::ResType::U) { features->u_count += 1; }
        }

        features->helices = p->helices();
        features->dotplot = v_.bp_probabilities(p->sequence());
        features->fe = v_.free_energy();
        features->multi_loops.resize(0);
        features->structure = v_.get_structure();

        for (auto const &m : p->motifs()) {
            if (m->mtype() == util::MotifType::NWAY) {
                features->multi_loops.push_back(m);
            }
        }


    }

    void
    FeatureGenerator::update_features(
            Features *features,
            secondary_structure::PoseOP const &p) {

        features->gc = 0;
        features->ua = 0;
        features->gu = 0;
        features->a_count = 0;
        features->c_count = 0;
        features->g_count = 0;
        features->u_count = 0;

        for (auto const &bp : p->basepairs()) {
            if (secondary_structure::is_gc_pair(bp)) { features->gc += 1; }
            else if (secondary_structure::is_au_pair(bp)) { features->ua += 1; }
            else if (secondary_structure::is_gu_pair(bp)) { features->gu += 1; }
        }

        for (auto const &r : p->residues()) {
            if (r->res_type() == secondary_structure::ResType::A) { features->a_count += 1; }
            else if (r->res_type() == secondary_structure::ResType::C) { features->c_count += 1; }
            else if (r->res_type() == secondary_structure::ResType::G) { features->g_count += 1; }
            else if (r->res_type() == secondary_structure::ResType::U) { features->u_count += 1; }
        }

        features->helices = p->helices();
        features->dotplot = v_.bp_probabilities(p->sequence());
        features->fe = v_.free_energy();
        features->multi_loops.resize(0);
        features->structure = v_.get_structure();

        for (auto const &m : p->motifs()) {
            if (m->mtype() == util::MotifType::NWAY) {
                features->multi_loops.push_back(m);
            }
        }
    }
}