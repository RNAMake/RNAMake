//
//  motif_graph.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resources/resource_manager.h"
#include "motif_data_structures/motif_graph.h"


void
MotifGraph::setup_options() {
    options_ = Options();
    options_.add_option(Option("sterics", 1));
    options_.add_option(Option("clash_radius", 2.9f));
}

void
MotifGraph::update_var_options() {
    sterics_              = options_.option<int>("sterics");
    clash_radius_         = options_.option<float>("clash_radius");
}

int
MotifGraph::add_motif(
    String const & m_name,
    int parent_index,
    String const & p_end_name) {
    
    auto parent = graph_.last_node();
    if(parent_index != -1) {
        parent = graph_.get_node(parent_index);
    }
    auto parent_end_index = parent->data()->end_index(p_end_name);
    return add_motif(m_name, parent_index, parent_end_index);
}

int
MotifGraph::add_motif(
    String const & m_name,
    String const & m_end_name,
    int parent_index,
    int parent_end_index) {
    
    auto m = MotifOP();
    try {
        m = ResourceManager::getInstance().get_motif(m_name, "", m_end_name);
    }
    catch(ResourceManagerException const & e) {
        throw MotifGraphException("failed to retrieve motif by name in add_motif: "
                                 + String(e.what()));
    }
    
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifGraph::add_motif(
    String const & m_name,
    int parent_index,
    int parent_end_index) {

    auto m = MotifOP();
    try {
        m = ResourceManager::getInstance().get_motif(m_name);
    }
    catch(ResourceManagerException const & e) {
        throw MotifGraphException("failed to retrieve motif by name in add_motif: "
                                   + String(e.what()));
    }
    
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifGraph::add_motif(
    MotifOP const & m,
    int parent_index,
    int parent_end_index) {
    
    auto parent = graph_.last_node();
    
    
    //catch out of bounds node index
    try {
        if(parent_index != -1) {
            parent = graph_.get_node(parent_index);
        }
    }
    catch(GraphException e) {
        throw MotifGraphException("could not add motif: " + m->name() + " with parent: "
                                  + std::to_string(parent_index) + "there is no node with" +
                                  "that index");
    }

    if(parent == nullptr) {
        auto m_copy = std::make_shared<Motif>(*m);
        m_copy->get_beads(m_copy->ends());
        m_copy->new_res_uuids();
        int pos = graph_.add_data(m_copy, -1, -1, -1, (int)m_copy->ends().size());
        merger_.add_motif(m_copy);
        return pos;
    }
    
    Ints avail_pos;
    try {
        avail_pos = graph_.get_available_pos(parent, parent_end_index);
    }
    catch(GraphException e) {
        throw MotifGraphException("could not add motif: " + m->name() + " with parent: "
                                  + std::to_string(parent_index));
    }
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->block_end_add()) { continue; }
        auto m_added = get_aligned_motif(parent->data()->ends()[p], m->ends()[0], m);
        if(sterics_ && _steric_clash(m_added)) { continue; }
        
        m_added->new_res_uuids();
        int pos = graph_.add_data(m_added, parent->index(), p, 0, (int)m_added->ends().size());
        merger_.add_motif(m_added, m_added->ends()[0], parent->data(), parent->data()->ends()[p]);
        return pos;
        
    }
    
    
    return 0;
}

void
MotifGraph::add_motif_tree(
    MotifTreeOP const & mt,
    int parent_index,
    String const & parent_end_name) {
    
    auto parent = graph_.last_node();
    if(parent_index != -1) {
        parent = graph_.get_node(parent_index);
    }
    auto parent_end_index = parent->data()->end_index(parent_end_name);
    int i = 0;
    for(auto const & n : *mt) {
        if(i == 0) { add_motif(n->data(), parent_index, parent_end_index); }
        else       { add_motif(n->data()); }
        i++;
    }
    
}


BasepairOP const &
MotifGraph::get_end(int pos) {
    auto n = graph_.get_node(pos);
    auto avail_pos = n->available_children_pos();
    assert(avail_pos.size() == 1 &&
           "called get_end with only node pos there are more then one ends");
    return n->data()->ends()[avail_pos[0]];
}

BasepairOP const &
MotifGraph::get_end(
        int pos,
        String const & end_name) {
    
    auto n = graph_.get_node(pos);
    auto end_index = n->data()->end_index(end_name);
    auto avail_pos = n->available_children_pos();
    
    assert(std::find(avail_pos.begin(), avail_pos.end(), end_index) != avail_pos.end() &&
           "end_name specified is not availiable");
    
    return n->data()->ends()[end_index];
}

int
MotifGraph::_steric_clash(MotifOP const & m) {
    float dist = 0;
    for(auto const & n : graph_) {
        for(auto const & c1 : n->data()->beads()) {
            if(c1.btype() == BeadType::PHOS) { continue; }
            for(auto const & c2 : m->beads()) {
                if(c2.btype() == BeadType::PHOS) { continue; }
                dist = c1.center().distance(c2.center());
                if(dist < clash_radius_) { return 1; }
            }
        }
    }
    return 0;
}

void
MotifGraph::write_pdbs(String const & fname) {
    std::stringstream ss;
    for( auto const & n : graph_.nodes()) {
        ss << fname << "." << n->index() << ".pdb";
        n->data()->to_pdb(ss.str());
        ss.str("");
    }
}

void
MotifGraph::remove_motif(int pos) {
    auto n = graph_.get_node(pos);
    merger_.remove_motif(n->data());
    graph_.remove_node(pos);
}

void
MotifGraph::replace_ideal_helices() {
    for(auto const & n : graph_) {
        if(n->data()->mtype() != MotifType::HELIX) { continue; }
        if(n->data()->residues().size() == 4) { continue; }
        
        auto parent = GraphNodeOP<MotifOP>(nullptr);
        auto parent_end_index = 0;
        auto other  = GraphNodeOP<MotifOP>(nullptr);
        auto other_end_index = 0;
        
        if(n->connections()[0] != nullptr) {
            parent = n->connections()[0]->partner(n->index());
            parent_end_index = n->connections()[0]->end_index(parent->index());
        }
        if(n->connections()[1] != nullptr) {
            other = n->connections()[1]->partner(n->index());
            other_end_index = n->connections()[1]->end_index(other->index());
        }
        
        auto name_spl = split_str_by_delimiter(n->data()->name(), ".");
        int count = 1;
        if(name_spl.size() == 3) {
            count = std::stoi(name_spl[2]);
        }
        
        remove_motif(n->index());
        auto h = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
        int pos = 0;
        if(parent == nullptr) { pos = _add_motif_to_graph(h); }
        else                  { pos = _add_motif_to_graph(h, parent->index(), parent_end_index); }
        
        
        int old_pos = pos;
        for(int j = 0; j < count; j++) {
            pos = _add_motif_to_graph(h, old_pos, 1);
            old_pos = pos;
        }
        
        if(other != nullptr) {
            graph_.connect(pos, other->index(), 1, other_end_index);
            auto node = graph_.get_node(pos);
            merger_.connect_motifs(node->data(), other->data(),
                                   node->data()->ends()[1],
                                   other->data()->ends()[other_end_index]);
        }
        
    }
    
}

int
MotifGraph::_add_motif_to_graph(
    MotifOP const & m,
    int parent_index,
    int parent_end_index) {
    
    if(parent_index == -1) {
        auto m_copy = std::make_shared<Motif>(*m);
        m_copy->new_res_uuids();
        m_copy->get_beads(m_copy->ends());
        
        int pos = graph_.add_data(m_copy, -1, -1, -1, (int)m_copy->ends().size(), 1);
        merger_.add_motif(m_copy);
        return pos;
    }
    
    else {
        auto parent = graph_.get_node(parent_index);
        auto m_added = get_aligned_motif(parent->data()->ends()[parent_end_index],
                                         m->ends()[0], m);
        m_added->new_res_uuids();
        int pos = graph_.add_data(m_added, parent->index(), parent_end_index, 0,
                                  (int)m_added->ends().size());
        merger_.add_motif(m_added, m_added->ends()[0], parent->data(),
                          parent->data()->ends()[parent_end_index]);
        return pos;
    }
    
}

void
MotifGraph::replace_helical_sequence(sstruct::PoseOP const & ss) {
    for(auto & n : graph_.nodes()) {
        if(n->data()->mtype() != MotifType::HELIX) { continue; }
        
        auto ss_m = ss->motif(n->data()->id());
        if(ss_m == nullptr) {
            throw MotifGraphException("could not find ss motif, cannot update helical sequence");
        }
        
        auto spl = split_str_by_delimiter(ss_m->end_ids()[0], "_");
        auto new_name = String();
        new_name += spl[0][0]; new_name += spl[2][1]; new_name += "=";
        new_name += spl[0][1]; new_name += spl[2][0];
        auto m = ResourceManager::getInstance().get_motif(new_name);
        auto org_res = n->data()->residues();
        auto new_res = m->residues();
        for(int i = 0; i < org_res.size(); i++) {
            new_res[i]->uuid(org_res[i]->uuid());
        }
        
        auto parent = GraphNodeOP<MotifOP>(nullptr);
        auto parent_end_index = 0;
        auto other  = GraphNodeOP<MotifOP>(nullptr);
        
        if(n->connections()[0] != nullptr) {
            parent = n->connections()[0]->partner(n->index());
            parent_end_index = n->connections()[0]->end_index(parent->index());
        }
        if(n->connections()[1] != nullptr) {
            other = n->connections()[1]->partner(n->index());
        }
        
        if(parent != nullptr) {
            auto m_added = get_aligned_motif(parent->data()->ends()[parent_end_index],
                                             m->ends()[0],
                                             m);
            n->data() = m_added;
        }
        else {
            n->data() = m;
        }
        
        
        if(other == nullptr) { continue; }
        if(other->data()->mtype() == MotifType::HELIX) { continue; }
        
        auto m_added = get_aligned_motif(n->data()->ends()[1],
                                         other->data()->ends()[0],
                                         other->data());
        other->data() = m_added;
        merger_.update_motif(other->data());
    }
    
}




































































