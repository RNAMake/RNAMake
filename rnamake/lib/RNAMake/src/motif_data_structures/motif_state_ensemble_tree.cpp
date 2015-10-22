//
//  motif_state_ensemble_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "data_structure/graph/graph_node.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "resources/resource_manager.h"
#include "util/cartesian_product.h"
#include "math/xyz_matrix.h"


int
MotifStateEnsembleTree::add_ensemble(
    MotifStateEnsembleOP const & ensemble,
    int parent_index,
    int parent_end_index) {
    
    auto parent = tree_.last_node();
    if(parent_index != -1) {
        parent = tree_.get_node(parent_index);
    }
    
    if(parent == nullptr) {
        return tree_.add_data(std::make_shared<MotifStateEnsemble>(ensemble->copy()),
                              ensemble->num_end_states());
    }
    
    auto avail_pos = tree_.get_available_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->block_end_add()) { continue; }
        
        return tree_.add_data(std::make_shared<MotifStateEnsemble>(ensemble->copy()),
                              ensemble->num_end_states(),
                              parent->index(),
                              p);
    }
    
    return -1;
    
}

MotifStateTreeOP
MotifStateEnsembleTree::to_mst() {
    
    auto mst = std::make_shared<MotifStateTree>();
    mst->option("sterics", 0);
    
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for(auto const & n : tree_) {
        i++;
        auto state = n->data()->most_populated();
        if(i == 0) {
            mst->add_state(state);
            continue;
        }
        
        parent_index = n->parent_index();
        parent_end_index = n->parent_end_index();
        j = mst->add_state(state, parent_index, parent_end_index);
        if(j == -1) {
            std::runtime_error("can not build motif state tree from mset");
        }
    }
    
    return mst;
}

void
MotifStateEnsembleTree::setup_from_mt(
    MotifTreeOP const & mt) {
    
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for(auto const & n : *mt) {
        i++;
        MotifStateEnsembleOP mse;
        try {
            mse = ResourceManager::getInstance().get_motif_state_ensemble(n->data()->end_ids()[0]);
        }
        //cannot find ensemble build one from motif
        catch(ResourceManagerException const & e) {
            auto m = ResourceManager::getInstance().get_motif(n->data()->name(),
                                                              n->data()->end_ids()[0]);
            mse = std::make_shared<MotifStateEnsemble>(m->get_state());
        }
        
        if(i == 0) {
            add_ensemble(mse);
        }
        else {
            parent_index = n->parent()->index();
            parent_end_index = n->parent_end_index();
            if(parent_end_index == -1) {
                std::runtime_error("cannot setup_from_mt in MotifStateEnsembleTree");
            }
            j = add_ensemble(mse, parent_index, parent_end_index);
            if(j == -1) {
                std::runtime_error("failed to add ensemble in setup_from_mt");
            }
        }
        
        
    }
    
    
}


void
MotifStateEnsembleTreeEnumerator::record() {
    
    std::vector<Ints> ranges(mtst_->size());
    for(int i = 0; i < mtst_->size(); i++) {
        int max = (int)mtst_->get_node(i)->data()->members().size();
        Ints range(max);
        for(int j = 0; j < max; j++) {
            range[j] = j;
        }
        ranges[i] = range;
    }
    
    auto mst = mtst_->to_mst();
    CartesianProduct<int> iterator(ranges);
    Ints c, last_combo;
    Matrix r;
    Point d;
    Vector euler;
    int j = 0;
    float mag;
    float r_diff, r_diff_flip, r_best;
    
    Matrix I = Matrix::identity();
    Matrix I_flip = I.get_flip_orientation();

    std::ofstream out("test.out");
    while(!iterator.end()) {
        c = iterator.next();
        if(last_combo.size() == 0) { last_combo = c; }
        
        for(int i = 0; i < c.size(); i++) {
            if(c[i] == last_combo[i]) { continue; }
            else{
                mst->replace_state(i, mtst_->get_node(i)->data()->members()[c[i]]->motif_state);
            }
         }
        
        d = mst->last_node()->data()->cur_state->end_states()[1]->d();
        r = mst->last_node()->data()->cur_state->end_states()[1]->r();
        mag = d.magnitude();
        
        r_diff = I.difference(r);
        r_diff_flip = I_flip.difference(r);
        
        r_best = r_diff;
        if(r_diff > r_diff_flip) {
            r_best = r_diff_flip;
        }
        
        calc_euler(r, euler);
        out << vector_to_str(d) << "," << matrix_to_str(r) << "," << mag << "," << r_diff << "," << r_best << "," << euler[0] << "," << euler[1] << "," << euler[2] << std::endl;

        j++;
        if(j % 1000 == 0) { std::cout << j << std::endl; }
        
        last_combo = c;
        
    }
    
    out.close();


    
    
}
























