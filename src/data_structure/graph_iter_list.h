//
// Created by Joseph Yesselman on 1/19/18.
//

#ifndef RNAMAKE_NEW_GRAPH_ITER_LIST_H
#define RNAMAKE_NEW_GRAPH_ITER_LIST_H

#include <queue>

#include <base/types.hpp>
#include <base/assertions.h>
#include <base/vector_container.h>
#include <data_structure/graph_base.h>
#include <data_structure/graph_adjacency_list.h>

namespace data_structure {

template<typename DataType, typename AdjacencyListType>
class IterList {
public:
    struct VisitedNode {
        std::shared_ptr<VisitedNode> parent;
        Index index;

        inline
        VisitedNode(
                std::shared_ptr<VisitedNode> nparent,
                Index nindex):
                parent(nparent),
                index(nindex) {}

        bool
        index_in_path(
                Index i) {
            if(index == i) { return true; }
            auto current = parent;
            while(current != nullptr) {
                if(current->index == i) { return true; }
                current = current->parent;
            }
            return false;
        }

        int
        path_length() {
            auto path_length = 1;
            auto current = parent;
            while(current != nullptr) {
                path_length += 1;
                current = current->parent;
            }
            return path_length;
        }
    };

    typedef std::shared_ptr<VisitedNode> VisitedNodeOP;

public:
    IterList():
            iter_list_(std::vector<Node<DataType> *>()),
    open_(std::queue<Index>()),
    seen_(std::map<Index, int>()) {
        iter_list_.reserve(100);
    }

public:
    typedef typename std::vector<Node<DataType> const *>::const_iterator const_iterator;
    typedef typename std::vector<Node<DataType> *>::iterator iterator;

    iterator begin() { return iter_list_.begin(); }
    iterator end()   { return iter_list_.end(); }

    const_iterator begin() const noexcept { return iter_list_.begin(); }
    const_iterator end() const noexcept   { return iter_list_.end(); }

public:
    virtual
    void
    transversal(
            AdjacencyListType & adj_list,
            Index start_n) {
        iter_list_.resize(0);

        seen_ = std::map<Index, int>();
        auto neighbors = std::vector<Index>();
        neighbors.reserve(10);

        open_.push(start_n);
        seen_[start_n] = 1;

        while(open_.size() > 0) {
            auto current = open_.front();
            open_.pop();
            Node<DataType> * p = &adj_list.get_node(current);
            iter_list_.push_back(p);
            get_neighbors(current, adj_list, neighbors);
            for(auto const & n : neighbors) { open_.push(n); }
        }

        // ensure you find all nodes not connected or in another sub graph
        for(auto const & kv : adj_list) {
            if(seen_.find(kv.first)  != seen_.end()) { continue; }
            auto new_start = kv.first;
            open_.push(new_start);
            seen_[new_start] = 1;

            while(open_.size() > 0) {
                auto current = open_.front();
                open_.pop();
                Node<DataType> * p = &adj_list.get_node(current);
                iter_list_.push_back(p);
                get_neighbors(current, adj_list, neighbors);
                for(auto const & n : neighbors) { open_.push(n); }
            }
        }
    }

    virtual
    void
    path_transversal(
            AdjacencyListType & adj_list,
            Index start_n,
            Index end_n) {
        iter_list_.resize(0);

        seen_ = std::map<Index, int>();
        auto neighbors = std::vector<Index>();
        neighbors.reserve(10);

        auto start_vn = std::make_shared<VisitedNode>(nullptr, start_n);
        auto current_round = std::vector<VisitedNodeOP>{start_vn};
        auto next_round = std::vector<VisitedNodeOP>();
        auto end_vn = VisitedNodeOP(nullptr);

        while(current_round.size() > 0) {
            for(auto vn : current_round) {
                get_neighbors_path(vn, adj_list, neighbors);
                for(auto const & ni : neighbors) {
                    auto new_vn = std::make_shared<VisitedNode>(vn, ni);
                    if(ni == end_n ) {
                        if(end_vn == nullptr) {
                            end_vn = new_vn;
                        }
                        else if(end_vn->path_length() > new_vn->path_length()) {
                            end_vn = new_vn;
                        }
                    }
                    else {
                        next_round.push_back(new_vn);
                    }
                }
            }
            current_round = next_round;
            next_round = std::vector<VisitedNodeOP>();
        }

        if(end_vn == nullptr) {
            throw GraphException(
                    "there is no path between nodes: " + std::to_string(start_n) +  " and " +
                    std::to_string(end_n) + " !");
        }

        auto path_back = std::vector<Index>();
        while(end_vn != nullptr) {
            path_back.push_back(end_vn->index);
            end_vn = end_vn->parent;
        }
        int pos = path_back.size()-1;
        for(int i = pos; i >= 0; i--) {
            iter_list_.push_back(&adj_list.get_node(path_back[i]));
        }

    }

protected:
    virtual
    void
    get_neighbors(
            Index ni,
            AdjacencyListType & adj_list,
            std::vector<Index> & neighbors) {
        neighbors.resize(0);
        auto & edges = adj_list.get_node_edges(ni);
        for(auto const & e : edges) {
            if(e == nullptr) { continue; }
            if(seen_.find(e->partner(ni)) != seen_.end()) { continue; }
            neighbors.push_back(e->partner(ni));
            seen_[e->partner(ni)];
        }
    }

    virtual
    void
    get_neighbors_path(
            VisitedNodeOP vn,
            AdjacencyListType & adj_list,
            std::vector<Index> & neighbors) {
        neighbors.resize(0);
        auto & edges = adj_list.get_node_edges(vn->index);
        for(auto const & e : edges) {
            if(e == nullptr) { continue; }
            if(vn->index_in_path(e->partner(vn->index))) { continue; }
            neighbors.push_back(e->partner(vn->index));
        }
    }

protected:
    std::vector<Node<DataType> *> iter_list_;
    std::queue<Index> open_;
    std::map<Index, int> seen_;

};

template<typename DataType, typename AdjacencyListType>
class DirectedIterList : public IterList<DataType, AdjacencyListType> {
public:
    typedef IterList<DataType, AdjacencyListType> BaseClass;
    typedef typename BaseClass::VisitedNodeOP VisitedNodeOP;

public:
    DirectedIterList():
            BaseClass() {}

public:
    virtual
    void
    transversal(
            AdjacencyListType & adj_list,
            Index start_n) {
        this->iter_list_.resize(0);

        this->seen_ = std::map<Index, int>();
        auto neighbors = std::vector<Index>();
        neighbors.reserve(10);

        this->open_.push(start_n);
        this->seen_[start_n] = 1;

        while(this->open_.size() > 0) {
            auto current = this->open_.front();
            this->open_.pop();
            Node<DataType> * p = &adj_list.get_node(current);
            this->iter_list_.push_back(p);
            get_neighbors(current, adj_list, neighbors);
            for(auto const & n : neighbors) { this->open_.push(n); }
        }

        // ensure you find all nodes not connected or in another sub graph
        bool found = 1;
        while(found) {
            found = 0;
            for (auto const & kv : adj_list) {
                if (this->seen_.find(kv.first) != this->seen_.end()) { continue; }
                if (adj_list.has_parent(kv.first)) { continue; }
                found = 1;
                auto new_start = kv.first;
                this->open_.push(new_start);
                this->seen_[new_start] = 1;

                while (this->open_.size() > 0) {
                    auto current = this->open_.front();
                    this->open_.pop();
                    Node<DataType> *p = &adj_list.get_node(current);
                    this->iter_list_.push_back(p);
                    get_neighbors(current, adj_list, neighbors);
                    for (auto const & n : neighbors) { this->open_.push(n); }
                }
            }
        }
    }

    void
    sub_graph_transversal(
            AdjacencyListType & adj_list,
            Index start_n,
            Index end_n) {
        BaseClass::path_transversal(adj_list, start_n, end_n);
        auto neighbors = std::vector<Index>();
        neighbors.reserve(10);
        int pos = 0;
        while(pos < this->iter_list_.size()) {
            get_neighbors(this->iter_list_[pos]->index(), adj_list, neighbors);
            for(auto const & ni : neighbors) {
                auto found = false;
                for(auto const & n : this->iter_list_) {
                    if(n->index() == ni) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    this->iter_list_.push_back(&adj_list.get_node(ni));
                }

            }

            pos += 1;
        }
    }

protected:
    virtual
    void
    get_neighbors(
            Index ni,
            AdjacencyListType & adj_list,
            std::vector<Index> & neighbors) {
        neighbors.resize(0);
        auto & edges = adj_list.get_node_edges(ni);
        for(auto const & e : edges) {
            if(e == nullptr) { continue; }
            auto partner = e->partner(ni);
            if(adj_list.has_parent(partner) && adj_list.get_parent_index(partner) == ni) {
                if(this->seen_.find(e->partner(ni)) != this->seen_.end()) { continue; }
                neighbors.push_back(e->partner(ni));
                this->seen_[e->partner(ni)];
            }
        }
    }

    virtual
    void
    get_neighbors_path(
            VisitedNodeOP vn,
            AdjacencyListType & adj_list,
            std::vector<Index> & neighbors) {
        neighbors.resize(0);
        auto & edges = adj_list.get_node_edges(vn->index);
        for(auto const & e : edges) {
            if(e == nullptr) { continue; }
            auto partner = e->partner(vn->index);
            if(adj_list.has_parent(partner) && adj_list.get_parent_index(partner) == vn->index) {
                if(vn->index_in_path(e->partner(vn->index))) { continue; }
                neighbors.push_back(e->partner(vn->index));
            }
        }
    }

};



}

#endif //RNAMAKE_NEW_GRAPH_ITER_LIST_H
