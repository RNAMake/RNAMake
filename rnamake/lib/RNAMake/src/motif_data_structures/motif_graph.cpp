//
//  motif_graph.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include <list>
#include <queue>          // std::queue

//RNAMake Headers
#include "resources/resource_manager.h"
#include "structure/residue_type_set_manager.h"
#include "motif_data_structures/motif_graph.h"


MotifGraph::MotifGraph():
    graph_(data_structure::graph::GraphStatic<MotifOP>()),
    merger_(nullptr),
    clash_radius_(2.5),
    sterics_(1),
    options_(base::Options()),
    update_merger_(1),
    update_align_list_(1),
    align_list_(data_structure::graph::GraphNodeOPs<MotifOP>()),
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
    auto spl = base::split_str_by_delimiter(s, "&");
    auto node_spl = base::split_str_by_delimiter(spl[0], "|");
    auto sspl = Strings();
    //int i = 0;
    int max_index = 0;
    for(auto const & n_spl : node_spl) {
        sspl = base::split_str_by_delimiter(n_spl, ",");
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
    
    auto con_spl = base::split_str_by_delimiter(spl[1], "|");
    for(auto const & c_str : con_spl) {
        auto c_spl = base::split_str_by_delimiter(c_str, ",");
        graph_.connect(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                       std::stoi(c_spl[2]), std::stoi(c_spl[3]));
    }
 
    options_.set_value("sterics", true);
    _align_motifs_all_motifs();
}

void
MotifGraph::_setup_from_str(String const & s) {
    options_.set_value("sterics", false);
    auto spl = base::split_str_by_delimiter(s, "FAF");
    auto node_spl = base::split_str_by_delimiter(spl[0], "KAK");
    int max_index = 0;
    for(auto const & n_str : node_spl) {
        if(n_str.length() < 10) { break; }
        auto n_spl = base::split_str_by_delimiter(n_str, "^");
        auto m = std::make_shared<Motif>(n_spl[0],
                                         structure::ResidueTypeSetManager::getInstance().residue_type_set());
        
        try {
            auto m2 = RM::instance().motif(m->name());
        } catch(ResourceManagerException const & e) {
            RM::instance().register_motif(m);
        }

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
    
    auto con_spl = base::split_str_by_delimiter(spl[1], "|");
    for(auto const & c_str : con_spl) {
        if(c_str.length() < 4) { break; }
        auto c_spl = base::split_str_by_delimiter(c_str, ",");
        graph_.connect(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                       std::stoi(c_spl[2]), std::stoi(c_spl[3]));
        
    }
    
    options_.set_value("sterics", true);
    
}


MotifGraph::MotifGraph(
    MotifGraph const & mg):
    options_(base::Options()),
    graph_(data_structure::graph::GraphStatic<MotifOP>(mg.graph_)) {
        
    // dear god this is horrible but cant figure out a better way to do a copy
    for(auto const & n : mg.graph_.nodes()) {
        graph_.get_node(n->index())->data() = std::make_shared<Motif>(*n->data());
    }
        
    options_ = base::Options(mg.options_);
    aligned_ = mg.aligned_;
    update_merger_ = 1;
    update_align_list_ = 1;
}


//setup helpers ////////////////////////////////////////////////////////////////////////////////////

void
MotifGraph::update_indexes(std::map<int, int> const & index_hash) {
    auto invert_hash = std::map<int, int>();
    for(auto const & kv : index_hash) {
        invert_hash[kv.second] = kv.first;
    }
    
    auto largest = 0;
    auto new_aligned = std::map<int, int>();
    auto nodes = data_structure::graph::GraphNodeOPs<MotifOP>();
    for(auto const & kv : invert_hash) {
        auto n = get_node(kv.first);
        nodes.push_back(n);
    }
    
    for(auto & n : nodes) {
        auto new_i = invert_hash[n->index()];
        auto aligned = aligned_[n->index()];
        
        if(new_i > largest) { largest = new_i; }
        n->index(new_i);
        new_aligned[new_i] = aligned;
    }
    
    aligned_ = new_aligned;
    update_merger_ = 1;
    update_align_list_ = 1;
    graph_.index(largest+1);
}


//add function helpers /////////////////////////////////////////////////////////////////////////////


data_structure::graph::GraphNodeOP<MotifOP>
MotifGraph::_get_parent(
    String const & m_name,
    int parent_index) {
    
    auto parent = graph_.last_node();
    
    //catch non existant parent
    try {
        if(parent_index != -1) { parent = graph_.get_node(parent_index); }
    }
    catch(data_structure::graph::GraphException const & e) {
        throw MotifGraphException(
            "could not add motif: " + m_name + " with parent index: " +
            std::to_string(parent_index) + "there is no node with that index");
    }
    
    return parent;
}

Ints
MotifGraph::_get_available_parent_end_pos(
    data_structure::graph::GraphNodeOP<MotifOP> const & parent,
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
    data_structure::graph::GraphNodeOP<MotifOP> const & parent,
    int parent_end_index) {
    
    m->new_res_uuids();
    int pos = -1;
    
    if(parent == nullptr) {
        pos = graph_.add_data(m, -1, -1, -1, (int)m->ends().size(), 1);
        if(pos != -1) {
            aligned_[pos] = 0;
            update_align_list_ = 1;
            update_merger_ = 1;
        }
    }
    
    else {
        pos = graph_.add_data(m, parent->index(), parent_end_index, 0,
                              (int)m->ends().size());
        
        if(pos != -1) {
            aligned_[pos] = 1;
            update_align_list_ = 1;
            update_merger_ = 1;
        }
        

    }
    
    
    return pos;
    
}

int
MotifGraph::_steric_clash(MotifOP const & m) {
    float dist = 0;
    for(auto const & n : graph_) {
        for(auto const & c1 : n->data()->beads()) {
            if(c1.btype() == structure::BeadType::PHOS) { continue; }
            for(auto const & c2 : m->beads()) {
                if(c2.btype() == structure::BeadType::PHOS) { continue; }
                dist = c1.center().distance(c2.center());
                if(dist < clash_radius_) { return 1; }
            }
        }
    }
    return 0;
}

int
MotifGraph::_get_connection_end(
    data_structure::graph::GraphNodeOP<MotifOP> const & node,
    String const & bp_name) {
    
    int node_end_index = -1;
    
    if(bp_name != "") {
        auto ei = node->data()->get_end_index(bp_name);
        
        if(!node->available_pos(ei)) {
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
    _update_align_list();
    auto non_aligned_nodes = unaligned_nodes();
    for(auto & n : align_list_) {
        if(std::find(non_aligned_nodes.begin(), non_aligned_nodes.end(), n) != non_aligned_nodes.end() ) { continue; }
        
        auto parent = n->connections()[0]->partner(n->index());
        auto pei = n->connections()[0]->end_index(parent->index());
        
        auto m_added = get_aligned_motif(parent->data()->ends()[pei],
                                         n->data()->ends()[0],
                                         n->data());
        n->data() = m_added;
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
    catch (structure::RNAStructureException const & e) {
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
    int parent_end_index,
    int orphan) {
    
    for(auto const & n : graph_.nodes()) {
        if(n->data()->id() == m->id()) {
            throw MotifGraphException(
                "cannot add motif there is already a motif in this graph with its unique "
                " indentifier");
        }
    }
    
    auto parent = _get_parent(m->name(), parent_index);

    if(parent == nullptr || orphan) {
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
    catch(data_structure::graph::GraphException const & e) {
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
    
    auto node_i = data_structure::graph::GraphNodeOP<MotifOP>(nullptr);
    auto node_j = data_structure::graph::GraphNodeOP<MotifOP>(nullptr);
    
    try {  node_i = graph_.get_node(i); }
    catch(data_structure::graph::GraphException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(i) +" does not exist");
    }
    
    try {  node_j = graph_.get_node(j); }
    catch(data_structure::graph::GraphException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(j) +" does not exist");
    }
    
    auto node_i_ei = _get_connection_end(node_i, i_bp_name);
    auto node_j_ei = _get_connection_end(node_j, j_bp_name);
    
    auto node_i_end_name = node_i->data()->end_name(node_i_ei);
    auto node_j_end_name = node_j->data()->end_name(node_j_ei);
    
    
    graph_.connect(i, j, node_i_ei, node_j_ei);
    update_merger_ = 1;
}


//remove functions /////////////////////////////////////////////////////////////////////////////////


void
MotifGraph::remove_motif(int pos) {
    auto n = graph_.get_node(pos);
    graph_.remove_node(pos);
    aligned_.erase(pos);
    update_merger_ = 1;
    update_align_list_ = 1;
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
            if(n->data()->mtype() != util::MotifType::HELIX) { continue; }
            if(n->data()->residues().size() == 4) { continue; }
            
            found = 1;
            
            auto parent = data_structure::graph::GraphNodeOP<MotifOP>(nullptr);
            auto parent_end_index = 0;
            auto other  = data_structure::graph::GraphNodeOP<MotifOP>(nullptr);
            auto other_end_index = 0;
            
            if(n->connections()[0] != nullptr) {
                parent = n->connections()[0]->partner(n->index());
                parent_end_index = n->connections()[0]->end_index(parent->index());
            }
            
            if(n->connections()[1] != nullptr) {
                other = n->connections()[1]->partner(n->index());
                other_end_index = n->connections()[1]->end_index(other->index());
            }

            //TODO look at resiude size instead of this mess!!
            auto name_spl = base::split_str_by_delimiter(n->data()->name(), ".");
            int count = 1;
            if(name_spl.size() == 3 && name_spl[1] != "AVG") {
                count = std::stoi(name_spl[2]);
            }
            else if(name_spl.size() == 3) {
                count = std::stoi(name_spl[2])-2;
            }
            else if(name_spl.size() == 4) {
                count = std::stoi(name_spl[2])-2;
            }
            
            auto old_n_aligned = aligned_[n->index()];
            auto old_n = n;
            
            remove_motif(n->index());
            int pos = 0;
            auto h = RM::instance().motif("HELIX.IDEAL");
            if(parent == nullptr) {
                h->get_beads(h->ends()[0]);
                pos = _add_motif_to_graph(h, nullptr, -1);
                
                if(old_n_aligned == 0) {
                    auto new_n = get_node(pos);
                    auto m_added = get_aligned_motif(old_n->data()->ends()[0],
                                                     new_n->data()->ends()[0],
                                                     new_n->data());
                    n->data() = m_added;
                }

            }
            
            else {
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
            }
            
            break;
            
        }
    }
    
    //re align all motifs
    update_align_list_ = 1;
    update_merger_ = 1;
    _align_motifs_all_motifs();
}

void
MotifGraph::replace_helical_sequence(secondary_structure::PoseOP const & ss) {
    for(auto & n : graph_.nodes()) {
        if(n->data()->mtype() != util::MotifType::HELIX) { continue; }
        
        auto ss_m = ss->motif(n->data()->id());
        if(ss_m == nullptr) {
            throw MotifGraphException("could not find ss motif, cannot update helical sequence");
        }
        
        //std::cout << ss_m->end_ids()[0] << " " << ss_m->sequence() << " " << new_name << std::endl;
        if(n->data()->end_ids()[0] == ss_m->end_ids()[0] && n->data()->name() != "HELIX.IDEAL") {
            continue;
        }
        auto m = RM::instance().bp_step(ss_m->end_ids()[0]);
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
        
        if(aligned_[n->index()] == 0) {
            auto m_added = get_aligned_motif(n->data()->ends()[0], m->ends()[0], m);
            n->data() = m_added;
        }
        else {
            n->data() = m;
        }
    }
    
    update_merger_ = 1;
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


//misc functions ///////////////////////////////////////////////////////////////////////////////////

void
MotifGraph::_update_align_list() {
    if(!update_align_list_) { return; }
    
    auto non_aligned_nodes = unaligned_nodes();
    auto open = std::queue<data_structure::graph::GraphNodeOP<MotifOP>>();
    auto used_nodes = std::map<data_structure::graph::GraphNodeOP<MotifOP>, int>();

    align_list_ = data_structure::graph::GraphNodeOPs<MotifOP>();
    
    for(auto const & start : non_aligned_nodes) {
        open.push(start);
        auto seen_nodes = std::map<data_structure::graph::GraphNodeOP<MotifOP>, int>();

        while (!open.empty()) {
            auto n = open.front();
            open.pop();
            
            seen_nodes[n] = 1;
            if(n->index() == start->index()) {
                align_list_.push_back(n);
            }
            else {
                // should this be block end?
                if(n->connections()[0] == nullptr) { continue; }
                auto c = n->connections()[0];
                auto parent = c->partner(n->index());
                if(used_nodes.find(parent) == used_nodes.end()) { continue; }

                // don't add unaligned nodes twice
                auto found = 0;
                for(auto const & na_node : non_aligned_nodes) {
                    if(na_node == n) { found = 1; break; }
                }
                if(found) { continue; }
                align_list_.push_back(n);
            }
            
            used_nodes[n] = 1;
            int i = -1;
            for(auto const & c : n->connections()) {
                i++;
                if(i == n->data()->block_end_add() || c == nullptr) { continue; }
                auto partner_n = c->partner(n->index());
                if(seen_nodes.find(partner_n) != seen_nodes.end() ||
                   used_nodes.find(partner_n) != used_nodes.end()) { continue; }
                if(c->end_index(partner_n->index()) == partner_n->data()->block_end_add()) {
                    open.push(partner_n);
                }
                else if(partner_n->data()->ends().size() == 1) {
                    open.push(partner_n);
                }
            }
        }
    }
    
    update_align_list_ = 0;
}

void
MotifGraph::_update_merger() {
    if(!update_merger_) { return; }
    
    //make sure align list is up to date
    _update_align_list();
    merger_ = std::make_shared<MotifMerger>();
    auto non_aligned_nodes = unaligned_nodes();
    auto seen_connections = std::map<data_structure::graph::GraphConnectionOP<MotifOP>, int>();
    
    for(auto const & n : align_list_) {
        if(std::find(non_aligned_nodes.begin(), non_aligned_nodes.end(), n) != non_aligned_nodes.end() ) {
            merger_->add_motif(n->data());
            continue;
        }
        
        auto c = n->connections()[0];
        auto parent = c->partner(n->index());
        auto pei = c->end_index(parent->index());
        seen_connections[c] = 1;
        merger_->add_motif(n->data(), n->data()->ends()[0],
                           parent->data(), parent->data()->ends()[pei]);
        
    }
    
    for(auto const & n : align_list_) {
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }
            if(seen_connections.find(c) != seen_connections.end()) { continue; }
            auto partner = c->partner(n->index());
            auto end1 = n->data()->ends()[c->end_index(n->index())];
            auto end2 = partner->data()->ends()[c->end_index(partner->index())];
            merger_->connect_motifs(n->data(), partner->data(), end1, end2);
            seen_connections[c] = 1;
        }
    }
    
    update_merger_ = 0;
}


//getters functions ////////////////////////////////////////////////////////////////////////////////


structure::BasepairOP const &
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

structure::BasepairOP const &
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

structure::BasepairOP
MotifGraph::get_available_end(
    String const & m_name,
    String const & end_name) {
    
    auto node = data_structure::graph::GraphNodeOP<MotifOP>(nullptr);
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
    
    auto end = structure::BasepairOP(nullptr);
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

data_structure::graph::GraphNodeOPs<MotifOP> const
MotifGraph::unaligned_nodes() const {
    auto nodes = data_structure::graph::GraphNodeOPs<MotifOP>();
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
    options_.add_option("sterics", true, base::OptionType::BOOL);
    options_.add_option("clash_radius", 2.9f, base::OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
MotifGraph::update_var_options() {
    sterics_              = options_.get_bool("sterics");
    clash_radius_         = options_.get_int("clash_radius");
}




























































