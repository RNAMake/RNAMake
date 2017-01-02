//
//  motif_state_ensemble_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "data_structure/graph/graph_node.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "motif/motif_ensemble.h"
#include "resources/resource_manager.h"
#include "util/cartesian_product.h"
#include "math/xyz_matrix.h"

MotifStateEnsembleTree::MotifStateEnsembleTree():
    connections_(MotifConnections()),
    tree_( TreeStatic<MotifStateEnsembleOP>()){}

MotifStateEnsembleTree::MotifStateEnsembleTree(
    MotifTreeOP const & mt):
    MotifStateEnsembleTree() {
    
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for(auto const & n : *mt) {
        i++;
        
        auto mse = MotifStateEnsembleOP();
        auto found_supplied = RM::instance().has_supplied_motif_ensemble(
                                    n->data()->name(), n->data()->end_name(0));
                
        if(n->data()->mtype() == MotifType::HELIX) {
            if(n->data()->residues().size() > 4) {
                throw MotifStateTreeEnsembleException(
                    "helix: " + n->data()->name() + " has more then 2 basepairs "
                    "must be broken up into basepair steps before being converted "
                    "into motif state ensemble");
            }
            
            
            try {
                mse = RM::instance().motif_state_ensemble(n->data()->end_ids()[0]);
            }
            catch(ResourceManagerException const & e) {
                throw MotifStateTreeEnsembleException(
                    "cannot find motif state ensemble for a basepair with id: " +
                    n->data()->end_ids()[0] + " is this a WC basepair?");
            }
        }
        
        // extra motif ensemble supplied by user
        else if(found_supplied) {
            std::cout << "MOTIF STATE ENSEMBLE TREE: found supplied ensemble for name=";
            std::cout << n->data()->name() << " endname=" << n->data()->end_name(0) << std::endl;
            mse = RM::instance().get_supplied_motif_ensemble(
                                    n->data()->name(), n->data()->end_name(0))->get_state();
        }
        
        else {
            mse = std::make_shared<MotifStateEnsemble>(n->data()->get_state());
        }
        
    
        if(i == 0) {
            add_ensemble(mse);
        }
        else {
            parent_index = n->parent()->index();
            parent_end_index = n->parent_end_index();
            if(parent_end_index == -1) {
                throw MotifStateTreeEnsembleException(
                    "cannot setup_from_mt in MotifStateEnsembleTree");
            }
            j = add_ensemble(mse, parent_index, parent_end_index);
            if(j == -1) {
                throw MotifStateTreeEnsembleException("failed to add ensemble in setup_from_mt");
            }
        }
    }
        
    for(auto const & c : mt->connections()) {
        connections_.add_connection(c->i(), c->j(), c->name_i(), c->name_j());
        
    }
}


MotifStateEnsembleTree::MotifStateEnsembleTree(
    MotifStateTreeOP const & mst):
    MotifStateEnsembleTree(mst->to_motif_tree()) {}

//add functions ////////////////////////////////////////////////////////////////////////////////////


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
        return tree_.add_data(std::make_shared<MotifStateEnsemble>(*ensemble),
                              ensemble->num_end_states());
    }
    
    auto avail_pos = tree_.get_available_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->block_end_add()) { continue; }
        
        return tree_.add_data(std::make_shared<MotifStateEnsemble>(*ensemble),
                              ensemble->num_end_states(),
                              parent->index(),
                              p);
    }
    
    return -1;
    
}

MotifStateTreeOP
MotifStateEnsembleTree::to_mst() {
    
    auto mst = std::make_shared<MotifStateTree>();
    mst->set_option_value("sterics", false);    
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for(auto const & n : tree_) {
        i++;
        auto state = n->data()->most_populated();
        state->new_uuids();
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
    mst->set_option_value("sterics", true);
    
    return mst;
}


//enumerator  functions ////////////////////////////////////////////////////////////////////////////

void
MotifStateEnsembleTreeEnumerator::record(
    String fname) {
    

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

    std::ofstream out(fname + ".out");
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
        out << vector_to_str(d) << "," << matrix_to_str(r) << "," << mag << "," << r_diff << "," << r_best << "," << euler[0] << "," << euler[1] << "," << euler[2] << " " <<
            mtst_->get_node(0)->data()->members()[c[0]]->energy << std::endl;

        j++;
        if(j % 1000 == 0) { std::cout << j << std::endl; }
        
        last_combo = c;
        
    }
    
    out.close();


    
    
}
























