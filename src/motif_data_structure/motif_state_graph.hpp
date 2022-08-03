//
//  motif_state_graph.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef motif_state_graph_hpp
#define motif_state_graph_hpp

#include "base/option.h"
#include "data_structure/graph/graph.h"
#include "motif_data_structure/motif_graph.h"
#include "motif_data_structure/motif_state_node.hpp"
#include <motif/motif_state_aligner.h>
#include <stdio.h>

namespace motif_data_structure {

class MotifStateGraphException : public std::runtime_error {
public:
  MotifStateGraphException(String const &message)
      : std::runtime_error(message) {}
};

class MotifStateGraph {
public:
  MotifStateGraph();

  MotifStateGraph(MotifGraphOP const &);

  MotifStateGraph(MotifStateGraph const &);

  ~MotifStateGraph() {}

public: // iterators
  typedef typename data_structure::graph::GraphNodeOPs<MSNodeDataOP>::iterator
      iterator;
  typedef
      typename data_structure::graph::GraphNodeOPs<MSNodeDataOP>::const_iterator
          const_iterator;

  iterator begin() {
    _update_align_list();
    return _align_list.begin();
  }

  iterator end() { return _align_list.end(); }

private:
  void _setup_from_mg(MotifGraphOP const &);

public:
  size_t size() { return _graph.size(); }

private: // add function helpers
  data_structure::graph::GraphNodeOP<MSNodeDataOP> _get_parent(String const &,
                                                               int);

  Indexes _get_available_parent_end_pos(
      data_structure::graph::GraphNodeOP<MSNodeDataOP> const &, int);

  int _get_parent_index_from_name(
      data_structure::graph::GraphNodeOP<MSNodeDataOP> const &, String const &);

  int _get_connection_end(
      data_structure::graph::GraphNodeOP<MSNodeDataOP> const &, String const &);

  inline int _steric_clash(MSNodeDataOP const &new_data) {
    float dist;
    for (auto const &n : _graph.nodes()) {
      for (auto const &b1 : n->data()->cur_state->beads()) {
        for (auto const &b2 : new_data->cur_state->beads()) {
          dist = b1.distance(b2);
          if (dist < _clash_radius) {
            return 1;
          }
        }
      }
    }
    return 0;
  }

public: // add functions
  int add_state(motif::MotifStateOP const &state, int parent_index = -1,
                int parent_end_index = -1, int orphan = 0);

  int add_state(motif::MotifStateOP const &state, int parent_index,
                String const &parent_end_name);

  void add_connection(int, int, String const &, String const &);

  void replace_state(int i, motif::MotifStateOP const &);

  void remove_state(int pos = -1);

public: // remove functions
  void remove_level(int level);

public: // graph wrappers
  inline data_structure::graph::GraphNodeOP<MSNodeDataOP> last_node() {
    return _graph.last_node();
  }

  inline data_structure::graph::GraphNodeOP<MSNodeDataOP> const &
  get_node(int i) const {
    return _graph.get_node(i);
  }

  inline data_structure::graph::GraphNodeOP<MSNodeDataOP> const
  get_node(util::Uuid const &uuid) const {
    for (auto const &n : _graph) {
      if (n->data()->uuid() == uuid) {
        return n;
      }
    }
    throw MotifStateGraphException(
        "cannot get node with uuid no motif has it in this tree");
  }

  inline data_structure::graph::GraphNodeOP<MSNodeDataOP> const
  get_node(String const &m_name) const {
    auto node = data_structure::graph::GraphNodeOP<MSNodeDataOP>(nullptr);
    for (auto const &n : _graph) {
      if (n->data()->name() == m_name) {
        if (node != nullptr) {
          throw MotifStateGraphException(
              "cannot get node with name: " + m_name +
              " there is more then one motif "
              "with this name");
        }

        node = n;
      }
    }

    if (node == nullptr) {
      throw MotifStateGraphException("cannot get node with name: " + m_name +
                                     " there is no motif in the tree with "
                                     "this name");
    }

    return node;
  }

  void increase_level() { return _graph.increase_level(); }

public: // motif graph wrappers
  MotifGraphOP to_motif_graph();

private: // misc functions
  void _update_align_list();

  void _align_states(int pos = -1);

public: // getters
  data_structure::graph::GraphNodeOPs<MSNodeDataOP> const
  unaligned_nodes() const;

public: // Options Wrappers
  inline float get_int_option(String const &name) {
    return _options.get_int(name);
  }

  inline float get_float_option(String const &name) {
    return _options.get_float(name);
  }

  inline String get_string_option(String const &name) {
    return _options.get_string(name);
  }

  inline bool get_bool_option(String const &name) {
    return _options.get_bool(name);
  }

  template <typename T>
  void set_option_value(String const &name, T const &val) {
    _options.set_value(name, val);
    update_var_options();
  }

private:
  void setup_options();

  void update_var_options();

private:
  data_structure::graph::GraphStatic<MSNodeDataOP> _graph;
  data_structure::graph::GraphNodeOPs<MSNodeDataOP> _align_list;
  motif::MotifStateAligner _aligner;
  base::Options _options;
  std::map<int, int> _aligned;
  int _update_align_list;
  // options
  float _clash_radius;
  bool _sterics;
};

typedef std::shared_ptr<MotifStateGraph> MotifStateGraphOP;

} // namespace motif_data_structure

#endif /* motif_state_graph_hpp */
