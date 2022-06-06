//
//  motif_state_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structure/motif_state_tree.h"
#include "resources/resource_manager.h"

namespace motif_data_structure {

MotifStateTree::MotifStateTree() :
        tree_(data_structure::tree::TreeStatic<MSNodeDataOP>()),
        aligner_(motif::MotifStateAligner()),
        queue_(std::queue<MotifStateTreeNodeOP>()),
        connections_(MotifConnections()),
        options_(base::Options()) {
    setup_options();
}

MotifStateTree::MotifStateTree(
        MotifTreeOP const & mt) :
        MotifStateTree() {

    set_option_value("sterics", mt->get_bool_option("sterics"));

    int i = -1;
    for (auto const & n : *mt) {
        i++;

        auto ms = resources::Manager::instance().motif_state(n->data()->name(),
                                                             n->data()->end_ids()[0],
                                                             n->data()->ends()[0]->get_name_str());

        ms->uuid(n->data()->id());
        if (i == 0) {
            add_state(ms);
        } else {
            int j = add_state(ms, n->parent_index(), n->parent_end_index());
            //std::cout << ms->name() << " " << n->parent_index() << " " << n->parent_end_index() << std::endl;
            if (j == -1) {
                throw MotifStateTreeException(
                        "could not convert motif tree to motif state tree");
            }
        }

    }

    for (auto const & c : mt->connections()) {
        connections_.add_connection(c->i(), c->j(), c->name_i(), c->name_j());
    }

}

MotifStateTree::MotifStateTree(
        MotifStateTree const
& mst):

aligner_ (motif::MotifStateAligner()),
queue_(std::queue<MotifStateTreeNodeOP>()),
options_(base::Options(mst.options_)),
connections_(MotifConnections(mst.connections_)),
tree_(data_structure::tree::TreeStatic<MSNodeDataOP>(mst.tree_)) {

    for (auto const & n : mst) {
        tree_.get_node(n->index())->data() = std::make_shared<MSNodeData>(*n->data());
    }

    update_var_options();
}

MotifStateTree::MotifStateTree(
        String const & s) :
        MotifStateTree() {

    set_option_value("sterics", false);

    auto spl = base::split_str_by_delimiter(s, "|");
    auto node_spl = base::split_str_by_delimiter(spl[0], " ");
    int i = -1;
    int pos = 0;
    for (auto const & e : node_spl) {
        i++;
        auto n_spl = base::split_str_by_delimiter(e, ",");
        auto ms = motif::MotifStateOP(nullptr);

        try {
            ms = resources::Manager::instance().motif_state(n_spl[0], n_spl[2], n_spl[1]);
        }
        catch (resources::ResourceManagerException const & e) {
            throw MotifStateTreeException(
                    String("could not get state did you forget to add it to the resource manager: ") +
                    e.what());
        }

        if (i == 0) {
            pos = add_state(ms);
        } else {
            pos = add_state(ms, std::stoi(n_spl[3]), std::stoi(n_spl[4]));
        }

        if (pos == -1) {
            throw MotifStateTreeException(
                    "failed to add " + ms->name() + " pos " + std::to_string(i) + " in the tree "
                            "during rebuild from string");
        }

    }

    if (spl.size() == 1) { return; }

    auto connection_spl = base::split_str_by_delimiter(spl[1], " ");
    for (auto const & c_str : connection_spl) {
        auto c_spl = base::split_str_by_delimiter(c_str, ",");
        connections_.add_connection(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                                    c_spl[2], c_spl[3]);
    }

    set_option_value("sterics", true);


}


//add function helpers /////////////////////////////////////////////////////////////////////////////

MotifStateTreeNodeOP
MotifStateTree::_get_parent(
        int parent_index) {

    auto parent = tree_.last_node();

    //catch non existant parent
    try {
        if (parent_index != -1) { parent = tree_.get_node(parent_index); }
    }
    catch (data_structure::tree::TreeException const & e) {
        throw MotifStateTreeException(
                "could not add motif state with parent index: " + std::to_string(parent_index) +
                "there is no node with that index");
    }

    return parent;
}

Indexes
MotifStateTree::_get_available_parent_end_pos(
        MotifStateTreeNodeOP const & parent,
        int parent_end_index) {

    auto avail_pos = Indexes();

    if (parent == nullptr) { return avail_pos; }

    if (parent_end_index != -1) {
        int avail = parent->available_pos(parent_end_index);
        if (avail == 0) {
            throw MotifStateTreeException(
                    "cannot add state to tree as the end with index: " +
                    std::to_string(parent_end_index) + " you are trying to "
                            "add it to is already filled or does not exist");
        }

        if (parent_end_index == parent->data()->block_end_add()) {
            throw MotifStateTreeException(
                    "cannot add motif: to tree as the parent_end_index"
                            " supplied is blocked see class Motif");
        }

        auto name = parent->data()->end_name(parent_end_index);
        if (connections_.in_connection(parent->index(), name)) {
            throw MotifStateTreeException(
                    "cannot add state to tree as the end "
                            "you are trying to add it to is in a connection");
        }

        avail_pos.push_back(parent_end_index);
    } else {
        auto avail_pos_temp = parent->available_children_pos();
        for (auto const & p : avail_pos_temp) {
            if (p == parent->data()->block_end_add()) { continue; }
            auto parent_end_name = parent->data()->end_name(p);
            if (connections_.in_connection(parent->index(), parent_end_name)) {
                continue;
            }
            avail_pos.push_back(p);
        }
    }

    return avail_pos;


}

int
MotifStateTree::_get_parent_index_from_name(
        MotifStateTreeNodeOP const & parent,
        String const & parent_end_name) {

    auto parent_end_index = -1;

    try {
        parent_end_index = parent->data()->get_end_index(parent_end_name);
    }
    catch (motif::MotifStateException) {
        throw MotifStateTreeException(
                "cannot find parent_end_name: " + parent_end_name + " while trying to "
                        "add a state to tree");
    }

    int avail = parent->available_pos(parent_end_index);
    if (!avail) {
        throw MotifStateTreeException(
                "cannot add state to tree as the end with name: " +
                parent_end_name + " you are trying to "
                        "add it to is already filled or does not exist");
    }

    if (parent_end_index == parent->data()->block_end_add()) {
        throw MotifStateTreeException(
                "cannot add state: to tree as the parent_name_index"
                        " supplied is blocked see class Motif");
    }

    if (connections_.in_connection(parent->index(), parent_end_name)) {
        throw MotifStateTreeException(
                "cannot add state to tree as the end with name " + parent_end_name +
                "you are trying to add it to is in a connection");
    }

    return parent_end_index;


}

int
MotifStateTree::_get_connection_end(
        MotifStateTreeNodeOP const & node,
        String const & bp_name) {

    int node_end_index = -1;

    if (bp_name != "") {
        auto ei = node->data()->get_end_index(bp_name);

        if (!node->available_pos(ei)) {
            throw MotifStateTreeException(
                    "cannot add connection with " + std::to_string(node->index()) + " and "
                            "end name " + bp_name + " as the end is blocked");
        }

        if (ei == node->data()->block_end_add()) {
            throw MotifStateTreeException(
                    "cannot add connection with " + std::to_string(node->index()) + " and "
                            "end name " + bp_name + " as the end is blocked");
        }

        if (connections_.in_connection(node->index(), bp_name)) {
            throw MotifStateTreeException(
                    "cannot add connection with " + std::to_string(node->index()) +
                    " and end name " + bp_name + " as this end is "
                            "already in a connection");
        }
        node_end_index = ei;

    } else {
        auto node_indexes = node->available_children_pos();
        node_indexes.erase(node_indexes.begin(), node_indexes.begin() + 1);

        if (node_indexes.size() > 1) {
            throw MotifStateTreeException(
                    "cannot connect nodes " + std::to_string(node->index()) + " its unclear "
                            " which ends to attach as there is more then one possibility");
        }

        if (node_indexes.size() == 0) {
            throw MotifStateTreeException(
                    "cannot connect nodes " + std::to_string(node->index()) + " there are "
                            "no ends free ends to attach too");
        }

        auto node_index_name = node->data()->end_name(node_indexes[0]);
        if (connections_.in_connection(node->index(), node_index_name)) {
            throw MotifStateTreeException(
                    "cannot connect nodes " + std::to_string(node->index()) + " there are "
                            "no ends free ends to attach too");
        }

        node_end_index = node_indexes[0];

    }

    return node_end_index;

}



//add functions ////////////////////////////////////////////////////////////////////////////////////


int
MotifStateTree::add_state(
        motif::MotifStateOP const & state,
        int parent_index,
        int parent_end_index) {

    for (auto const & n : tree_) {
        if (n->data()->uuid() == state->uuid()) {
            throw MotifStateTreeException(
                    "cannot add state: " + state->name() + " to tree as its uuid is " +
                    "already present in the tree");
        }
    }

    auto parent = _get_parent(parent_index);

    if (parent == nullptr) {
        auto n_data = std::make_shared<MSNodeData>(state);
        return tree_.add_data(n_data, (int) state->end_states().size(), -1, -1);
    }

    auto avail_pos = _get_available_parent_end_pos(parent, parent_end_index);

    for (auto const & p : avail_pos) {
        auto n_data = std::make_shared<MSNodeData>(state);

        get_aligned_motif_state(parent->data()->cur_state->end_states()[p],
                                n_data->cur_state,
                                n_data->ref_state);

        if (sterics_ && _steric_clash(n_data)) { continue; }

        return tree_.add_data(n_data, (int) state->end_states().size(), parent->index(), p);

    }

    return -1;
}

int
MotifStateTree::add_state(
        motif::MotifStateOP const & state,
        int parent_index,
        String const & parent_end_name) {

    auto parent = _get_parent(parent_index);
    auto parent_end_index = _get_parent_index_from_name(parent, parent_end_name);
    return add_state(state, parent_index, parent_end_index);
}

int
MotifStateTree::add_mst(
        MotifStateTreeOP const & mst,
        int parent_index,
        int parent_end_index) {

    auto parent = _get_parent(parent_index);
    auto avail_pos = _get_available_parent_end_pos(parent, parent_end_index);

    if (avail_pos.size() == 0 && tree_.size() != 0) {
        throw MotifStateTreeException(
                "cannot add motif state tree to current tree, parent has no available "
                        "ends to add to");
    } else if (avail_pos.size() == 0) {
        avail_pos.push_back(-1);
    }

    int i = -1;
    int j = 0;
    auto index_dict = std::map<int, int>();
    for (auto const & n : *mst) {
        i++;

        if (i == 0) {
            j = add_state(n->data()->ref_state, parent_index, avail_pos[0]);
        } else {
            int ind = index_dict[n->parent_index()];
            int pei = n->parent_end_index();
            j = add_state(n->data()->ref_state, ind, pei);
        }

        index_dict[n->index()] = j;
        if (j == -1) {
            throw MotifStateTreeException("could not add motif state tree to this tree");
        }
    }

    return j;


}

int
MotifStateTree::add_mst(
        MotifStateTreeOP const & mst,
        int parent_index,
        String const & parent_end_name) {

    auto parent = _get_parent(parent_index);
    auto parent_end_index = _get_parent_index_from_name(parent, parent_end_name);
    return add_mst(mst, parent_index, parent_end_index);
}

void
MotifStateTree::add_connection(
        int i,
        int j,
        String const & i_bp_name,
        String const & j_bp_name) {

    auto node_i = MotifStateTreeNodeOP(nullptr);
    auto node_j = MotifStateTreeNodeOP(nullptr);

    try { node_i = tree_.get_node(i); }
    catch (data_structure::tree::TreeException) {
        throw MotifTreeException(
                "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
                std::to_string(i) + " does not exist");
    }

    try { node_j = tree_.get_node(j); }
    catch (data_structure::tree::TreeException) {
        throw MotifTreeException(
                "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
                std::to_string(j) + " does not exist");
    }

    auto node_i_ei = _get_connection_end(node_i, i_bp_name);
    auto node_j_ei = _get_connection_end(node_j, j_bp_name);

    auto node_i_end_name = node_i->data()->end_name(node_i_ei);
    auto node_j_end_name = node_j->data()->end_name(node_j_ei);


    connections_.add_connection(i, j, node_i_end_name, node_j_end_name);
}

void
MotifStateTree::replace_state(
        int i,
        motif::MotifStateOP const & new_state) {

    auto n = tree_.get_node(i);
    if (new_state->end_states().size() != n->data()->ref_state->end_states().size()) {
        throw MotifStateTreeException(
                "attempted to replace a state with a different number of ends");
    }

    auto old_state = n->data()->ref_state;
    n->data()->ref_state = new_state;
    n->data()->cur_state = std::make_shared<motif::MotifState>(*new_state);
    n->data()->uuid(old_state->uuid());

    for (int i = 0; i < new_state->end_states().size(); i++) {
        if (connections_.in_connection(n->index(), old_state->end_names()[i])) {
            connections_.update_connection_name(n->index(),
                                                old_state->end_names()[i],
                                                new_state->end_names()[i]);
        }

    }


    queue_.push(n);
    MotifStateTreeNodeOP current, parent;
    int pei;
    while (!queue_.empty()) {
        current = queue_.front();
        queue_.pop();

        parent = current->parent();
        if (parent == nullptr) { continue; }
        pei = current->parent_end_index();
        aligner_.get_aligned_motif_state(parent->data()->cur_state->end_states()[pei],
                                         current->data()->cur_state,
                                         current->data()->ref_state);

        for (auto const & c : current->children()) {
            if (c != nullptr) { queue_.push(c); }
        }

    }
}



//removal functions ////////////////////////////////////////////////////////////////////////////////



void
MotifStateTree::remove_node(
        int i) {

    if (i == -1) {
        i = last_node()->index();
    }

    try {
        auto n = get_node(i);
        tree_.remove_node(n);
        connections_.remove_connections_to(i);
    }
    catch (MotifStateTreeException) {
        throw MotifStateTreeException(
                "cannot remove node with index: " + std::to_string(i) + " as it does not exist");
    }
}

void
MotifStateTree::remove_node_level(
        int level) {

    if (level == -1) { level = tree_.level(); }

    auto remove = std::vector<MotifStateTreeNodeOP>();
    for (auto const & n : tree_) {
        if (n->level() >= level) {
            remove.push_back(n);
        }
    }

    std::reverse(remove.begin(), remove.end());

    for (auto const & n : remove) {
        remove_node(n->index());
    }

    // find new highest node index
    int max_index = 0;
    for (auto const & n : tree_) {
        if (n->index() > max_index) { max_index = n->index(); }
    }
    tree_.index(max_index);

}


//outputing functions //////////////////////////////////////////////////////////////////////////////


String
MotifStateTree::topology_to_str() {
    String s;

    for (auto const & n : tree_) {
        s += n->data()->ref_state->name() + "," + n->data()->ref_state->end_names()[0] + ",";
        s += n->data()->ref_state->end_ids()[0] + "," + std::to_string(n->parent_index()) + ",";
        s += std::to_string(n->parent_end_index()) + " ";
    }
    s += "|";
    for (auto const & c : connections_) {
        s += c->to_str() + " ";
    }

    return s;
}


MotifTreeOP
MotifStateTree::to_motif_tree() {

    auto mt = std::make_shared<MotifTree>();
    //mt->set_option_value("sterics", options_.get_bool("sterics"));
    mt->set_option_value("sterics", false);
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for (auto const & n : tree_) {
        i++;
        motif::MotifOP m;
        /*if(n->data()->ref_state->name() != "") {
            m = resources::Manager::instance().motif(n->data()->ref_state->name(), n->data()->ref_state->end_ids()[0]);
        }
        else {
            m = resources::Manager::instance().motif("", n->data()->ref_state->end_ids()[0]);
        }*/
        m = resources::Manager::instance().motif(n->data()->name(), "",
                                                 n->data()->end_name(0));
        m->new_res_uuids();

        if (i == 0) {
            align_motif(n->data()->cur_state->end_states()[0],
                        m->ends()[0], m);
            mt->add_motif(m);
            continue;
        }

        parent_index = n->parent_index();
        parent_end_index = n->parent_end_index();
        if (parent_end_index == -1) {
            throw MotifStateTreeException("cannot convert to motif tree");
        }

        j = mt->add_motif(m, parent_index, parent_end_index);
        if (j == -1) {
            throw MotifStateTreeException("failed to add motif in to_motif_tree");
        }


    }

    for (auto const & c : connections_) {
        mt->add_connection(c->i(), c->j(), c->name_i(), c->name_j());
    }
    mt->set_option_value("sterics", options_.get_bool("sterics"));

    return mt;
}


//option functions /////////////////////////////////////////////////////////////////////////////////


void
MotifStateTree::setup_options() {
    options_.add_option("sterics", true, base::OptionType::BOOL);
    options_.add_option("clash_radius", 2.9f, base::OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
MotifStateTree::update_var_options() {
    sterics_ = options_.get_bool("sterics");
    clash_radius_ = options_.get_int("clash_radius");
}

}















