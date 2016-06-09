//
//  motif_graph.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

//RNAMake Headers
#include "resources/resource_manager.h"
#include "structure/residue_type_set_manager.h"
#include "motif_data_structures/motif_graph.h"

MotifGraph::MotifGraph(String const & s):
graph_(GraphStatic<MotifOP>()),
merger_(MotifMerger()),
clash_radius_(2.5),
options_(Options("MotifGraphOptions")),
sterics_(0),
aligned_(std::map<int, int>())
{
    auto spl = split_str_by_delimiter(s, "&");
    auto node_spl = split_str_by_delimiter(spl[0], "|");
    Strings sspl;
    int i = 0;
    auto pos = 0;
    for(auto const & n_spl : node_spl) {
        sspl = split_str_by_delimiter(n_spl, ",");
        auto m = RM::instance().motif(sspl[0], sspl[2], sspl[1]);

        if(i == 0 ) {
            pos = add_motif(m);
        }
        else {
            pos = add_motif(m, std::stoi(sspl[3]), sspl[4]);
        }
        assert(pos != -1 && "cannot build motif_graph from topology");
        i++;
    }
    
    setup_options();
    
    /*auto spl_indices = split_str_by_delimiter(spl[1], ",");
    auto indices = Ints();
    for(auto const & s : spl_indices) {
        indices.push_back(std::stoi(s));
    }
    graph_.set_node_indexes(indices);*/
    
    if(spl.size() == 2) { return; }
    
    auto connection_spl = split_str_by_delimiter(spl[2], "|");
    for(auto const & c_spl : connection_spl) {
        sspl = split_str_by_delimiter(c_spl, " ");
        add_connection(std::stoi(sspl[0]), std::stoi(sspl[1]), sspl[2], sspl[3]);
    }
}


MotifGraph::MotifGraph(
    String const & s,
    MotifGraphStringType const & s_type):
graph_(GraphStatic<MotifOP>()),
merger_(MotifMerger()),
clash_radius_(2.5),
options_(Options("MotifGraphOptions")),
sterics_(0),
aligned_(std::map<int, int>()) {
    setup_options();

    if      (s_type == MotifGraphStringType::OLD) { MotifGraph(s); }
    else if (s_type == MotifGraphStringType::MG)  { _setup_from_str(s); }
    else if (s_type == MotifGraphStringType::TOP) { _setup_from_top_str(s); }
}

void
MotifGraph::_setup_from_top_str(String const & s) {
    options_.set_value("sterics", false);
    auto spl = split_str_by_delimiter(s, "&");
    auto node_spl = split_str_by_delimiter(spl[0], "|");
    auto sspl = Strings();
    int i = 0;
    int max_index = 0;
    for(auto const & n_spl : node_spl) {
        sspl = split_str_by_delimiter(n_spl, ",");
        auto m = MotifOP(nullptr);
        try {
            m = RM::instance().motif(sspl[0], "", sspl[1]);
        } catch(...) { }
        if(m == nullptr) {
            m = RM::instance().motif(sspl[0]);
  
        }
        m->get_beads(m->ends());
        m->new_res_uuids();
        graph_.add_data(m, -1, -1, -1, (int)m->ends().size(), 1,
                        std::stoi(sspl[2]));
        aligned_[std::stoi(sspl[2])] = std::stoi(sspl[3]);
        
        if(std::stoi(sspl[2]) > max_index) {
            max_index = std::stoi(sspl[2]);
        }
    }
    
    graph_.index(max_index+1);
    
    auto con_spl = split_str_by_delimiter(spl[1], "|");
    for(auto const & c_str : con_spl) {
        auto c_spl = split_str_by_delimiter(c_str, ",");
        graph_.connect(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                       std::stoi(c_spl[2]), std::stoi(c_spl[3]));
    }
    
    int start = -1;
    for(auto const kv : aligned_) {
        if(kv.second == 0) { start = kv.first; }
    }
    
    auto n = GraphNodeOP<MotifOP>();
    for(auto it = graph_.transverse(graph_.get_node(start));
        it != graph_.end();
        ++it) {
        
        n = (*it);
        if(n->index() == start) {
            merger_.add_motif(n->data());
            continue;
        }
        if(n->connections()[0] == nullptr) { continue; }
        auto c = n->connections()[0];
        auto parent = c->partner(n->index());
        auto parent_end_index = c->end_index(parent->index());
        auto m_added = get_aligned_motif(parent->data()->ends()[parent_end_index],
                                         n->data()->ends()[0],
                                         n->data());
        n->data() = m_added;
        merger_.add_motif(n->data(), n->data()->ends()[0],
                          parent->data(), parent->data()->ends()[parent_end_index]);
    }
    
    
}

void
MotifGraph::_setup_from_str(String const & s) {
    options_.set_value("sterics", false);
    auto spl = split_str_by_delimiter(s, "FAF");
    auto node_spl = split_str_by_delimiter(spl[0], "KAK");
    int max_index = 0;
    for(auto const & n_str : node_spl) {
        if(n_str.length() < 10) { break; }
        auto n_spl = split_str_by_delimiter(n_str, "^");
        auto m = std::make_shared<Motif>(n_spl[0],
                                         ResidueTypeSetManager::getInstance().residue_type_set());
        graph_.add_data(m, -1, -1, -1, (int)m->ends().size(), 1, std::stoi(n_spl[1]));
        aligned_[std::stoi(n_spl[1])] = std::stoi(n_spl[2]);
        if(std::stoi(n_spl[2]) > max_index) {
            max_index = std::stoi(n_spl[2]);
        }

    }

    graph_.index(max_index+1);
    
    auto con_spl = split_str_by_delimiter(spl[1], "|");
    for(auto const & c_str : con_spl) {
        if(c_str.length() < 4) { break; }
        auto c_spl = split_str_by_delimiter(c_str, ",");
        graph_.connect(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                       std::stoi(c_spl[2]), std::stoi(c_spl[3]));
        
    }
    
    int start = -1;
    for(auto const kv : aligned_) {
        if(kv.second == 0) { start = kv.first; }
    }
    
    auto n = GraphNodeOP<MotifOP>();
    for(auto it = graph_.transverse(graph_.get_node(start));
        it != graph_.end();
        ++it) {
        
        n = (*it);
        
        if(n->index() == start) {
            merger_.add_motif(n->data());
            continue;
        }
        if(n->connections()[0] == nullptr) { std::cout << "no parent: " << n->index() << std::endl;
            continue;}
        
        auto c = n->connections()[0];
        auto parent = c->partner(n->index());
        auto parent_end_index = c->end_index(parent->index());
        merger_.add_motif(n->data(), n->data()->ends()[0],
                          parent->data(), parent->data()->ends()[parent_end_index]);
    }

}

void
MotifGraph::setup_options() {
    options_.add_option("sterics", true, OptionType::BOOL);
    options_.add_option("clash_radius", 2.9f, OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
MotifGraph::update_var_options() {
    sterics_              = options_.get_bool("sterics");
    clash_radius_         = options_.get_int("clash_radius");
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
    MotifOP const & m,
    int parent_index,
    String const & p_end_name) {
    
    auto parent = graph_.last_node();
    if(parent_index != -1) {
        parent = graph_.get_node(parent_index);
    }
    auto parent_end_index = parent->data()->end_index(p_end_name);
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifGraph::add_motif(
    String const & m_name,
    String const & m_end_name,
    int parent_index,
    int parent_end_index) {
    
    auto m = MotifOP();
    try {
        m = RM::instance().motif(m_name, "", m_end_name);
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
        m = RM::instance().motif(m_name);
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
        aligned_[pos] = 0;
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
        if(pos != -1) {
            merger_.add_motif(m_added, m_added->ends()[0], parent->data(), parent->data()->ends()[p]);
        }
        aligned_[pos] = 1;
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
    int pos = 0;
    for(auto const & n : *mt) {
        auto m = RM::instance().motif(n->data()->name(), n->data()->end_ids()[0]);
        if(i == 0) { pos = add_motif(n->data(), parent_index, parent_end_index); }
        else       { pos = add_motif(n->data()); }
        i++;
    }
    
}

void
MotifGraph::add_motif_tree(
    MotifTreeOP const & mt,
    int parent_index) {
    
    int parent_end_index = -1;
    if(parent_index != -1) {
        auto parent = graph_.get_node(parent_index);
        auto avail_pos = parent->available_children_pos();
        assert(avail_pos.size() > 0 && "cannot add_motif_tree to this node no free ends");
        parent_end_index = avail_pos[0];
    }
    
    int i = 0, j = 0;
    auto index_hash = std::map<int, int>();
    for(auto const & n : *mt) {
        if(i == 0) {
            j = add_motif(n->data(), parent_index, parent_end_index);
        }
        else       {
            int pi = index_hash[n->parent_index()];
            j = add_motif(n->data(), pi, n->parent_end_index());
        }
        index_hash[n->index()] = j;
        i++;
    }
}


void
MotifGraph::add_connection(
    int i,
    int j,
    String const & i_bp_name) {
    
    auto node_i = graph_.get_node(i);
    auto node_j = graph_.get_node(j);
    
    auto ei = node_i->data()->end_index(i_bp_name);
    assert(node_i->available_pos(ei) &&
           "cannot add_connection");
    
    auto node_j_indexes = node_j->available_children_pos();
    assert(node_j_indexes.size() != 0 &&
           "cannot add_connection no available pos");
    
    graph_.connect(i, j, ei, node_j_indexes[0]);
    
    merger_.connect_motifs(node_i->data(), node_j->data(),
                           node_i->data()->ends()[ei],
                           node_j->data()->ends()[node_j_indexes[0]]);

}

void
MotifGraph::add_connection(
    int i,
    int j,
    String const & i_bp_name,
    String const & j_bp_name) {
    
    auto node_i = graph_.get_node(i);
    auto node_j = graph_.get_node(j);
    
    auto name_i = String("");
    auto name_j = String("");
    auto ei = -1;
    auto ej = -1;

    
    if (i_bp_name != "") {
        ei = node_i->data()->end_index(i_bp_name);
        assert(node_i->available_pos(ei) && "cannot add_connection");
        name_i = i_bp_name;
    }
    else {
        auto node_i_indexes = node_i->available_children_pos();
        assert(node_i_indexes.size() != 0 && "cannot add_connection no available spots");
        ei = node_i_indexes[0];
        name_i = node_i->data()->ends()[ei]->name();
    }
    
    if (j_bp_name != "") {
        ej = node_j->data()->end_index(j_bp_name);
        assert(node_j->available_pos(ej) && "cannot add_connection");
        name_j = j_bp_name;
    }
    else {
        auto node_j_indexes = node_j->available_children_pos();
        assert(node_j_indexes.size() != 0 && "cannot add_connection no available spots");
        ej = node_j_indexes[0];
        name_j = node_j->data()->ends()[ej]->name();
    }
        
    graph_.connect(i, j, ei, ej);
    
    merger_.connect_motifs(node_i->data(), node_j->data(),
                           node_i->data()->ends()[ei],
                           node_j->data()->ends()[ej]);
    
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

String
MotifGraph::topology_to_str() {
    String s = "", con_str = "";
    String key1 = "", key2 = "";
    std::map<String, int> seen_connections;
    for(auto const & n : graph_.nodes()) {
        s += n->data()->name() + "," + n->data()->ends()[0]->name() + ",";
        s += std::to_string(n->index()) + "," + std::to_string(aligned_[n->index()]) + "|";
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }
            key1 = std::to_string(c->node_1()->index()) + " " + std::to_string(c->node_2()->index());
            key2 = std::to_string(c->node_2()->index()) + " " + std::to_string(c->node_1()->index());
            if(seen_connections.find(key1) != seen_connections.end() ||
               seen_connections.find(key2) != seen_connections.end()) {
                continue;
            }
            seen_connections[key1] = 1;
            con_str += std::to_string(c->node_1()->index()) + "," + std::to_string(c->node_2()->index()) + ",";
            con_str += std::to_string(c->end_index_1()) + "," + std::to_string(c->end_index_2()) + "|";
        }
    }
    s += "&";
    s += con_str;
    return s;
}

String
MotifGraph::to_str() {
    String s = "", con_str = "";
    String key1 = "", key2 = "";
    std::map<String, int> seen_connections;
    for(auto const & n : graph_.nodes()) {
        s += n->data()->to_str() + "^" + std::to_string(n->index()) + "^";
        s += std::to_string(aligned_[n->index()]) + " KAK ";
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }
            key1 = std::to_string(c->node_1()->index()) + " " + std::to_string(c->node_2()->index());
            key2 = std::to_string(c->node_2()->index()) + " " + std::to_string(c->node_1()->index());
            if(seen_connections.find(key1) != seen_connections.end() ||
               seen_connections.find(key2) != seen_connections.end()) {
                continue;
            }
            seen_connections[key1] = 1;
            con_str += std::to_string(c->node_1()->index()) + "," + std::to_string(c->node_2()->index()) + ",";
            con_str += std::to_string(c->end_index_1()) + "," + std::to_string(c->end_index_2()) + "|";
        }
    }
    s += " FAF ";
    s += con_str;
    return s;

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
    aligned_.erase(pos);
}

void
MotifGraph::remove_level(int level) {
    int pos = 0;
    while(pos < graph_.nodes().size()) {
        auto n = graph_.nodes()[pos];
        if(n->level() >= level) {
            remove_motif(n->index());
            continue;
        }
        pos++;
    }
}

void
MotifGraph::replace_ideal_helices() {
    int found = 1;
    while(found) {
        found =0;
        for(auto const & n : graph_) {
            if(n->data()->mtype() != MotifType::HELIX) { continue; }
            if(n->data()->residues().size() == 4) { continue; }
            
            found = 1;
            
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
            auto h = RM::instance().motif("HELIX.IDEAL");
            int pos = 0;
            if(parent == nullptr) {
                pos = _add_motif_to_graph(h);
                aligned_[pos] = 0 ;
            }
            else                  {
                pos = _add_motif_to_graph(h, parent->index(), parent_end_index);
                aligned_[pos] = 1 ;
            }
            
            int old_pos = pos;
            for(int j = 0; j < count; j++) {
                pos = _add_motif_to_graph(h, old_pos, 1);
                aligned_[pos] = 1;
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
        if(n->data()->name() == new_name) {
            continue;
        }
        auto m = RM::instance().motif(new_name);
        m->id(n->data()->id());
        auto org_res = n->data()->residues();
        auto new_res = m->residues();
        auto new_bps = m->basepairs();
        for(int i = 0; i < org_res.size(); i++) {
            new_res[i]->uuid(org_res[i]->uuid());
        }
        
        for(int i = 0; i < n->data()->basepairs().size(); i++) {
            new_bps[i]->uuid(n->data()->basepairs()[i]->uuid());
        }
        
        n->data() = m;

    }
    
    _align_motifs_all_motifs();
    
}

GraphNodeOPs<MotifOP>
MotifGraph::unaligned_nodes() {
    auto nodes = GraphNodeOPs<MotifOP>();
    for(auto const & kv : aligned_) {
        if(kv.second == 0) { nodes.push_back(get_node(kv.first)); }
    }
    return nodes;
}

void
MotifGraph::_align_motifs_all_motifs() {
    int start = -1;
    for(auto const kv : aligned_) {
        if(kv.second == 0) { start = kv.first; }
    }
    
    auto n = GraphNodeOP<MotifOP>();
    for(auto it = graph_.transverse(graph_.get_node(start));
        it != graph_.end();
        ++it) {
        
        n = (*it);
        if(n->index() == start) {
            merger_.add_motif(n->data());
            continue;
        }
        if(n->connections()[0] == nullptr) { continue; }
        auto c = n->connections()[0];
        auto parent = c->partner(n->index());
        auto parent_end_index = c->end_index(parent->index());
        auto m_added = get_aligned_motif(parent->data()->ends()[parent_end_index],
                                         n->data()->ends()[0],
                                         n->data());
        n->data() = m_added;
        merger_.update_motif(n->data());
    }
}
































































