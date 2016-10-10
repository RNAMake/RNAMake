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


MotifGraph::MotifGraph():
    graph_(GraphStatic<MotifOP>()),
    merger_(std::make_shared<MotifMerger>()),
    clash_radius_(2.5),
    sterics_(1),
    options_(Options()),
    aligned_(std::map<int, int>()) {  setup_options(); }


MotifGraph::MotifGraph(
    String const & s,
    MotifGraphStringType const & s_type) : MotifGraph() {

    if      (s_type == MotifGraphStringType::MG)  { _setup_from_str(s); }
    else if (s_type == MotifGraphStringType::TOP) { _setup_from_top_str(s); }
    else {
        throw MotifGraphException("cannot destringify motif graph wrong string type");
    }
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
        auto m = RM::instance().motif(sspl[0], "", sspl[1]);

        m->get_beads(m->ends());
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
            merger_->add_motif(n->data());
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
        merger_->add_motif(n->data(), n->data()->ends()[0],
                          parent->data(), parent->data()->ends()[parent_end_index]);
    }
    
    options_.set_value("sterics", true);

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
        
        /*try {
            auto m2 = RM::instance().motif(m->name());
        } catch(ResourceManagerException const & e) {
            RM::instance().add_motif(m);
        }*/

    

        if(m->ends().size() > 1) {
            m->get_beads(m->ends()[0]);
        }
        graph_.add_data(m, -1, -1, -1, (int)m->ends().size(), 1, std::stoi(n_spl[1]));
        aligned_[std::stoi(n_spl[1])] = std::stoi(n_spl[2]);
        if(std::stoi(n_spl[1]) > max_index) {
            max_index = std::stoi(n_spl[1]);
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
    
    auto seen_connections = std::map<String, int>();
    auto n = GraphNodeOP<MotifOP>();
    for(auto it = graph_.transverse(graph_.get_node(start));
        it != graph_.end();
        ++it) {
        
        n = (*it);
        
        if(n->index() == start) {
            merger_->add_motif(n->data());
            continue;
        }
        if(n->connections()[0] == nullptr) { std::cout << "no parent: " << n->index() << std::endl;
            continue;}
        
        
        auto c = n->connections()[0];
        auto parent = c->partner(n->index());
        auto parent_end_index = c->end_index(parent->index());
        seen_connections[std::to_string(n->index()) + " " + std::to_string(parent->index())] = 1;
        merger_->add_motif(n->data(), n->data()->ends()[0],
                          parent->data(), parent->data()->ends()[parent_end_index]);
        
    }
    
    // catch connections not used in alignement for for chain connections
    for(auto const & n : graph_) {
        for(auto const & c : n->connections()) {
            if(c != nullptr) { continue; }
            auto partner = c->partner(n->index());
            auto key1 = std::to_string(n->index()) + " " + std::to_string(partner->index());
            auto key2 = std::to_string(partner->index()) + " " + std::to_string(n->index());
            if(seen_connections.find(key1) != seen_connections.end()) { continue; }
            if(seen_connections.find(key2) != seen_connections.end()) { continue; }
            auto end1 = n->data()->ends()[c->end_index(n->index())];
            auto end2 = partner->data()->ends()[c->end_index(partner->index())];
            merger_->connect_motifs(n->data(), partner->data(), end1, end2);
            seen_connections[key1] = 1;
            
        }
    }
    
    options_.set_value("sterics", true);

    
}


MotifGraph::MotifGraph(
    MotifGraph const & mg):
    options_(Options()),
    graph_(GraphStatic<MotifOP>(mg.graph_)) {
        
    auto motifs = MotifOPs();
    // dear god this is horrible but cant figure out a better way to do a copy
    for(auto const & n : mg.graph_.nodes()) {
        graph_.get_node(n->index())->data() = std::make_shared<Motif>(*n->data());
        motifs.push_back(graph_.get_node(n->index())->data());
    }
        
    options_ = Options(mg.options_);
    merger_ = std::make_shared<MotifMerger>(*mg.merger_, motifs);
    aligned_ = mg.aligned_;
}


//add function helpers /////////////////////////////////////////////////////////////////////////////


GraphNodeOP<MotifOP>
MotifGraph::_get_parent(
    String const & m_name,
    int parent_index) {
    
    auto parent = graph_.last_node();
    
    //catch non existant parent
    try {
        if(parent_index != -1) { parent = graph_.get_node(parent_index); }
    }
    catch(GraphException const & e) {
        throw MotifGraphException(
            "could not add motif: " + m_name + " with parent index: " +
            std::to_string(parent_index) + "there is no node with that index");
    }
    
    return parent;
}

Ints
MotifGraph::_get_available_parent_end_pos(
    GraphNodeOP<MotifOP> const & parent,
    int parent_end_index) {
    
    auto avail_pos = Ints();
    
    if(parent_end_index != -1) {
        int avail = parent->available_pos(parent_end_index);
        if(!avail) {
            throw MotifGraphException(
                "cannot add motif to tree as the end with index: " +
                std::to_string(parent_end_index) +" you are trying to "
                "add it to is already filled or does not exist");
        }
        
        if(parent_end_index == parent->data()->block_end_add()) {
            throw MotifGraphException(
                "cannot add motif: to tree as the parent_end_index"
                " supplied is blocked see class Motif");
        }
    
        avail_pos.push_back(parent_end_index);
    }
    
    else {
        auto avail_pos_temp = parent->available_children_pos();
        for(auto const & p : avail_pos_temp) {
            if(p == parent->data()->block_end_add()) { continue; }
            avail_pos.push_back(p);
        }
    }
    
    return avail_pos;
    
}

int
MotifGraph::_add_motif_to_graph(
    MotifOP & m,
    GraphNodeOP<MotifOP> const & parent,
    int parent_end_index) {
    
    m->new_res_uuids();
    int pos = -1;
    
    if(parent == nullptr) {
        pos = graph_.add_data(m, -1, -1, -1, (int)m->ends().size(), 1);
        if(pos != -1) {
            merger_->add_motif(m);
            aligned_[pos] = 0;
        }
    }
    
    else {
        pos = graph_.add_data(m, parent->index(), parent_end_index, 0,
                              (int)m->ends().size());
        merger_->add_motif(m, m->ends()[0], parent->data(),
                           parent->data()->ends()[parent_end_index]);

        aligned_[pos] = 1;

    }
    
    
    return pos;
    
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

int
MotifGraph::_get_connection_end(
    GraphNodeOP<MotifOP> const & node,
    String const & bp_name) {
    
    int node_end_index = -1;
    
    if(bp_name != "") {
        auto ei = node->data()->get_end_index(bp_name);
        
        if(!node->available_pos(ei)) {
            throw MotifGraphException(
                "cannot add connection with " + std::to_string(node->index()) + " and "
                "end name " + bp_name + " as the end is blocked");
        }
        
        if(ei == node->data()->block_end_add()) {
            throw MotifGraphException(
                "cannot add connection with " + std::to_string(node->index()) + " and "
                "end name " + bp_name + " as the end is blocked");
        }

        node_end_index = ei;
        
    }
    
    else {
        auto node_indexes = node->available_children_pos();
        
        if(node_indexes.size() > 1) {
            throw MotifGraphException(
                "cannot connect nodes " + std::to_string(node->index()) + " its unclear "
                " which ends to attach as there is more then one possibility");
        }
        
        if(node_indexes.size() == 0) {
            throw MotifGraphException(
                "cannot connect nodes " + std::to_string(node->index())  + " there are "
                "no ends free ends to attach too");
        }
        
        auto node_index_name = node->data()->end_name(node_indexes[0]);
        node_end_index = node_indexes[0];
        
    }
    
    return node_end_index;
    
}

void
MotifGraph::_align_motifs_all_motifs() {
    int start = -1;
    for(auto const kv : aligned_) {
        if(kv.second == 0) {
            start = kv.first;
            break;
        }
    }
    
    if(start == -1) {
        throw MotifGraphException(
            "could not find start position in graph to perform alignment!");
    }
    
    auto n = GraphNodeOP<MotifOP>();
    for(auto it = graph_.transverse(graph_.get_node(start));
        it != graph_.end();
        ++it) {
        
        n = (*it);
        if(n->index() == start) {
            merger_->update_motif(n->data());
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
        merger_->update_motif(n->data());
    }
}



//add functions ////////////////////////////////////////////////////////////////////////////////////

int
MotifGraph::add_motif(
    MotifOP const & m,
    int parent_index,
    String const & p_end_name) {
    
    auto parent = _get_parent(m->name(), parent_index);
    auto parent_end_index = 0;
    
    if(parent == nullptr) {
        throw MotifGraphException(
            "cannot add motif: " + m->name() + " as there are no open parents to add too");
    }
    
    try {
        parent_end_index = parent->data()->get_end_index(p_end_name);
    }
    catch (RNAStructureException const & e) {
        throw MotifGraphException(
            "cannot find parent_end_name: " + p_end_name + " in "
            "parent motif: " + parent->data()->name());
    }
    
    if(parent_end_index == m->block_end_add()) {
        throw MotifGraphException(
            "cannot add motif: to graph as the parent_end_name"
            " supplied is blocked see class Motif");
    }
    
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifGraph::add_motif(
    MotifOP const & m,
    int parent_index,
    int parent_end_index) {
    
    for(auto const & n : graph_.nodes()) {
        if(n->data()->id() == m->id()) {
            throw MotifGraphException(
                "cannot add motif there is already a motif in this graph with its unique "
                " indentifier");
        }
    }
    
    auto parent = _get_parent(m->name(), parent_index);

    if(parent == nullptr) {
        auto m_copy = std::make_shared<Motif>(*m);
        m_copy->get_beads(m_copy->ends()[0]);
        return _add_motif_to_graph(m_copy, nullptr, -1);
    }
    
    auto avail_pos = _get_available_parent_end_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        auto m_added = get_aligned_motif(parent->data()->ends()[p], m->ends()[0], m);
        if(sterics_ && _steric_clash(m_added)) { continue; }
        
        auto pos = _add_motif_to_graph(m_added, parent, p);
        if(pos != -1) { return pos; }
        
    }
    
    return -1;
}

void
MotifGraph::add_motif_tree(
    MotifTreeOP const & mt,
    int parent_index,
    String const & parent_end_name) {
    
    auto parent = _get_parent("motif_tree", parent_index);
    auto parent_end_index = parent->data()->get_end_index(parent_end_name);
    _add_motif_tree(mt, parent_index, parent_end_index);
}

void
MotifGraph::add_motif_tree(
    MotifTreeOP const & mt,
    int parent_index,
    int parent_end_index) {
    
    _add_motif_tree(mt, parent_index, parent_end_index);
}

void
MotifGraph::_add_motif_tree(
    MotifTreeOP const & mt,
    int parent_index,
    int parent_end_index) {
    
    auto parent = _get_parent("motif_tree", parent_index);
    
    auto avail_pos = Ints();
    try { avail_pos = graph_.get_available_pos(parent, parent_end_index); }
    catch(GraphException const & e) {
        throw MotifGraphException("could not add motif_tree with parent: "
                                  + std::to_string(parent_index));
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
    String const & i_bp_name,
    String const & j_bp_name) {
    
    auto node_i = GraphNodeOP<MotifOP>(nullptr);
    auto node_j = GraphNodeOP<MotifOP>(nullptr);
    
    try {  node_i = graph_.get_node(i); }
    catch(TreeException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(i) +" does not exist");
    }
    
    try {  node_j = graph_.get_node(j); }
    catch(TreeException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(j) +" does not exist");
    }
    
    auto node_i_ei = _get_connection_end(node_i, i_bp_name);
    auto node_j_ei = _get_connection_end(node_j, j_bp_name);
    
    auto node_i_end_name = node_i->data()->end_name(node_i_ei);
    auto node_j_end_name = node_j->data()->end_name(node_j_ei);
    
    
    graph_.connect(i, j, node_i_ei, node_j_ei);
    
    merger_->connect_motifs(node_i->data(), node_j->data(),
                            node_i->data()->ends()[node_i_ei],
                            node_j->data()->ends()[node_j_ei]);

}


//remove functions /////////////////////////////////////////////////////////////////////////////////


void
MotifGraph::remove_motif(int pos) {
    auto n = graph_.get_node(pos);
    merger_->remove_motif(n->data());
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


//designing functions //////////////////////////////////////////////////////////////////////////////

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
            int pos = 0;
            auto h = RM::instance().motif("HELIX.IDEAL");
            if(parent == nullptr) {
                h->get_beads(h->ends()[0]);
                pos = _add_motif_to_graph(h, nullptr, -1);
            }
            else                  {
                auto m_added = get_aligned_motif(parent->data()->ends()[parent_end_index],
                                                 h->ends()[0], h);
                pos = _add_motif_to_graph(m_added, parent, parent_end_index);
            }
            
            
            int old_pos = pos;
            for(int j = 0; j < count; j++) {
                auto h = RM::instance().motif("HELIX.IDEAL");
                auto parent = get_node(old_pos);
                auto m_added = get_aligned_motif(parent->data()->ends()[1],
                                                 h->ends()[0], h);
                pos = _add_motif_to_graph(m_added, parent, 1);
                old_pos = pos;
            }
            
      
            if(other != nullptr) {
                graph_.connect(pos, other->index(), 1, other_end_index);
                auto node = graph_.get_node(pos);
                merger_->connect_motifs(node->data(), other->data(),
                                        node->data()->ends()[1],
                                        other->data()->ends()[other_end_index]);
            }
            
        }
    }
    
    //re align all motifs
    _align_motifs_all_motifs();
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
        //std::cout << ss_m->end_ids()[0] << " " << ss_m->sequence() << " " << new_name << std::endl;
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


//outputting functions /////////////////////////////////////////////////////////////////////////////


void
MotifGraph::write_pdbs(String const & fname) {
    std::stringstream ss;
    for( auto const & n : graph_.nodes()) {
        ss << fname << "." << n->index() << ".pdb";
        n->data()->to_pdb(ss.str());
        ss.str("");
    }
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


//getters functions ////////////////////////////////////////////////////////////////////////////////


BasepairOP const &
MotifGraph::get_available_end(int pos) {
    auto n = graph_.get_node(pos);
    auto avail_pos = n->available_children_pos();
   
    if(avail_pos.size() == 1) {
        throw MotifGraphException(
            "attempting to get end with pos " + std::to_string(pos) +
            " but this node has no free ends to build from");
    }
    
    return n->data()->ends()[avail_pos[1]];
}

BasepairOP const &
MotifGraph::get_available_end(
    int pos,
    String const & end_name) {
    
    auto n = graph_.get_node(pos);
    auto end_index = n->data()->get_end_index(end_name);
    if(!n->available_pos(end_index)) {
        throw MotifGraphException(
            "attempting to get end with pos " + std::to_string(pos) + "and end_name " + end_name +
            " but the end with this end name is not available to build from");
    }
    
    return n->data()->ends()[end_index];
}

BasepairOP 
MotifGraph::get_available_end(
    String const & m_name,
    String const & end_name) {
    
    auto node = GraphNodeOP<MotifOP>(nullptr);
    for(auto const & n : graph_.nodes()) {
        if(n->data()->name() == m_name) {
            if(node != nullptr) {
                throw MotifGraphException(
                    "attempting to get end with m_name: " + m_name + " and end name " + end_name +
                    " but there are MULTIPLE nodes with this name, dont know which to pick");
            }
            node = n;
        }
    }
    
    if(node == nullptr) {
        throw MotifGraphException(
            "attempting to get end with m_name: " + m_name + " and end name " + end_name +
            " but there are NO nodes with this name");

    }
    
    auto end = BasepairOP(nullptr);
    for(auto const & e : node->data()->ends()) {
        if(e->name() == end_name) {
            if(end != nullptr) {
                throw MotifGraphException(
                    "attempting to get end with m_name: " + m_name + " and end name " + end_name +
                    " but there are MULTIPLE ends with this name in this motif");
            }
            end = e;
        }
    }
    
    if(end == nullptr) {
        throw MotifGraphException(
            "attempting to get end with m_name: " + m_name + " and end name " + end_name +
            " but there are NO ends with this name");
        
    }
    
    auto end_index = node->data()->get_end_index(end_name);
    if(!node->available_pos(end_index)) {
        throw MotifGraphException(
            "attempting to get end with m_name " + m_name + "and end_name " + end_name +
            " but the end with this end name is not available to build from");
    }
    
    
    return end;
    
}


GraphNodeOPs<MotifOP> const
MotifGraph::unaligned_nodes() const {
    auto nodes = GraphNodeOPs<MotifOP>();
    for(auto const & kv : aligned_) {
        if(kv.second == 0) { nodes.push_back(get_node(kv.first)); }
    }
    return nodes;
}


MotifGraph::_MotifGraphBuildPointOPs
MotifGraph::get_build_points() {
    auto build_points = _MotifGraphBuildPointOPs();
    for(auto const & n : graph_.nodes()) {
        int i = -1;
        for(auto c : n->connections()) {
            i++;
            if(c == nullptr && i != 0) {
                build_points.push_back(std::make_shared<_MotifGraphBuildPoint>(n, i));
            }
        }
    }
    return build_points;
}


//option functions ////////////////////////////////////////////////////////////////////////////////


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




























































