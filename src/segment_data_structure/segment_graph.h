//
// Created by Joseph Yesselman on 1/14/18.
//

#ifndef RNAMAKE_NEW_SEGMENT_GRAPH_H
#define RNAMAKE_NEW_SEGMENT_GRAPH_H

#include <type_traits>

#include <data_structure/graph/graph.h>
#include <resource_management/resource_manager.h>
#include <structure/all_atom/aligner.hpp>
//#include <structure/aligner.h>

using namespace data_structure::graph;
namespace segment_data_structure {

template <typename SegmentType, typename AlignerType> class SegmentGraph {
public:
  typedef std::shared_ptr<SegmentType> SegmentTypeOP;

public: // construction ////////////////////////////////////////////////////
  SegmentGraph()
      : _aligner(AlignerType()),
        _graph(FixedEdgeDirectedGraph<SegmentTypeOP>()), _needs_update(false),
        _first_update(INT_MAX) {}

  SegmentGraph(SegmentGraph const &sg)
      : _aligner(AlignerType()), _graph(sg._graph), _needs_update(false),
        _first_update(INT_MAX) {
    _update_default_transveral();
  }

  ~SegmentGraph() = default;

public: // iteration ///////////////////////////////////////////////////////
  typedef typename Indexes::iterator iterator;
  typedef typename Indexes::const_iterator const_iterator;

  iterator begin() noexcept { return _path.begin(); }
  iterator end() noexcept { return _path.end(); }

  const_iterator begin() const noexcept { return _path.begin(); }
  const_iterator end() const noexcept { return _path.end(); }

public: // operators ///////////////////////////////////////////////////////
  const SegmentType &operator[](Index ni) { return *_graph.get_node_data(ni); }

public:
  [[nodiscard]] inline size_t get_num_segments() const {
    return _graph.get_num_nodes();
  }

  [[nodiscard]] inline Connections const &
  get_segment_connections(Index ni) const {
    return _graph.get_node_edges(ni);
  }

public:
  [[nodiscard]] inline bool has_parent(Index ni) const {
    return _graph.has_parent(ni);
  }

  [[nodiscard]] inline Index get_parent_index(Index ni) const {
    return _graph.get_parent_index(ni);
  }

  [[nodiscard]] inline Index get_parent_end_index(Index ni) const {
    return _graph.get_parent_end_index(ni);
  }

  [[nodiscard]] inline bool are_segments_connected(Index n1, Index n2) const {
    return _graph.edge_between_nodes(n1, n2);
  }

public:
  Index add(SegmentTypeOP &seg) {
    Index ni = _graph.add_node(seg, seg->get_num_ends());
    _update_default_transveral();
    return ni;
  }

  Index add(SegmentTypeOP &seg, Index parent_index,
            String const &parent_end_name) {
    const auto &parent = _graph.get_node_data(parent_index);
    auto parent_end_index = parent->get_end_index(parent_end_name);
    _aligner.align(*parent, *seg, parent_end_index);
    //structure::all_atom::align_segment(*parent, *seg, parent_end_index);
    auto ni = _graph.add_node(seg, seg->get_num_ends(), 0,
                              ConnectionPoint{parent_index, parent_end_index});
    _update_default_transveral();
    return ni;
  }

public:
  void add_connection(const ConnectionPoint &nie1,
                      const ConnectionPoint &nie2) {
    return _graph.add_edge(nie1, nie2);
  }

public:
  void remove(Index ni) {
    _graph.remove_node(ni);
    auto roots = _graph.get_root_indexes();
    if (roots.size() > 0) {
      _graph.setup_transversal(roots[0]);
    }
  }

  void replace(Index pos, SegmentTypeOP &seg, bool copy = true) {
    if (!copy) {
      _graph.get_node_data(pos) = seg;
    } else {
      _graph.get_node_data(pos) = SegmentTypeOP(seg);
    }
    _needs_update = true;
    if (_first_update > pos) {
      _first_update = pos;
    }
    _update_alignments(_first_update);
  }

public:
  [[nodiscard]] const String &get_end_name(Index n_index,
                                           Index end_index) const {
    return _graph.get_node_data(n_index)->get_end(end_index).get_name();
  }

  [[nodiscard]] const FixedEdgeDirectedGraph<SegmentTypeOP> get_graph() {
    return _graph;
  }

  [[nodiscard]] SegmentTypeOP get_node_data(Index pos) {
    return _graph.get_node_data(pos);
  }

private:
  void _update_default_transveral() {
    auto roots = _graph.get_root_indexes();
    if (roots.size() > 0) {
      _graph.setup_transversal(roots[0]);
    }
    _path = _graph.get_index_path();
  }

  void _update_alignments(int start) {
    auto parent_index = 0;
    auto parent_end_index = 0;
    auto needs_update = false;
    for (const auto &n : _graph) {
      if (!has_parent(n->get_index())) {
        needs_update = false;
        continue;
      }
      if (start == n->get_index()) {
        needs_update = true;
      }
      if (!needs_update) {
        continue;
      }
      parent_index = get_parent_index(n->get_index());
      parent_end_index = get_parent_end_index(n->get_index());

      _aligner.align(*_graph.get_node_data(parent_index),
                     *_graph.get_node_data(n->get_index()), parent_end_index);
    }
  }

private:
  AlignerType _aligner;
  FixedEdgeDirectedGraph<SegmentTypeOP> _graph;
  bool _needs_update;
  int _first_update;
  Indexes _path;
};

typedef SegmentGraph<structure::all_atom::Segment, structure::all_atom::Aligner>
    SegmentGraphAllAtom;

// template <typename SegmentType, typename AlignerType>
// using SegmentGraphOP = std::shared_ptr<SegmentGraph<SegmentType,
//  AlignerType>>;

// typedefs
//////////////////////////////////////////////////////////////////////////////////

// conversion functions
/////////////////////////////////////////////////////////////////////////////////

/*
namespace __helpers {

template <typename SegmentType, typename AlignerType>
void add_motif(Index i, SegmentGraph<SegmentType, AlignerType> const &g,
               SegmentGraph<SegmentType, AlignerType> &new_g,
               std::map<Index, Index> &index_convert) {

  auto pos = -1;
  auto &seg = g.get_segment(i);

  if (g.has_parent(i)) {
    auto pi = g.get_parent_index(i);
    auto pei = g.get_parent_end_index(i);
    auto new_pi = index_convert[pi];
    pos =
        new_g.add_segment(seg, new_pi, new_g.get_segment_end_name(new_pi, pei));
  }

  else {
    pos = new_g.add_segment(seg);
  }

  index_convert[i] = pos;
};

template <typename SegmentType, typename AlignerType>
base::VectorContainerOP<data_structure::Edge>
get_extra_connections(SegmentGraph<SegmentType, AlignerType> const &g,
                      std::map<String, int> &seen_connections) {

  auto new_edges = std::vector<data_structure::Edge>();

  for (auto const &n : g) {
    auto &connections = g.get_motif_connections(n->index());
    for (auto const &c : connections) {
      if (c == nullptr) {
        continue;
      }
      auto key = c->to_str();
      if (seen_connections.find(key) != seen_connections.end()) {
        continue;
      }

      new_edges.push_back(*c);

      seen_connections[key] = 1;
    }
  }

  return std::make_shared<base::VectorContainer<data_structure::Edge>>(
      new_edges);
};

} // namespace __helpers

// template <typename SegmentType, typename AlignerType>
// SegmentGraphOP<SegmentType, AlignerType>
// convert_ideal_helices_to_basepair_steps(
//         SegmentGraph<SegmentType, AlignerType> const & g,
//         resources::ResourceManager const & rm) {
//     typedef data_structure::NodeIndexandEdge NodeIndexandEdge;
//     auto new_g = std::make_shared<SegmentGraph<SegmentType, AlignerType>>();
//     auto index_convert = std::map<Index, Index>();
//     auto args = StringStringMap{{"name", "HELIX.IDEAL"}};
//     auto aligner = AlignerType();
//     auto seen_connections = std::map<String, int>();
//     for(auto const & n : g) {
//         if(g.has_parent(n->index())) {
//             auto pi = g.get_parent_index(n->index());
//             auto pei = g.get_parent_end_index(n->index());
//
//             auto key1 = data_structure::Edge(n->index(), pi, 0,
//             pei).to_str(); auto key2 = data_structure::Edge(pi, n->index(),
//             pei, 0).to_str();
//
//             seen_connections[key1] = 1; seen_connections[key2] = 1;
//
//         }
//
//         // not a helix or not idealized
//         if(n->data().get_segment_type() != util::SegmentType::HELIX) {
//             __helpers::add_motif(n->index(), g, *new_g, index_convert);
//             continue;
//         }
//         else if(n->data().get_name_str().substr(0, 5) != "HELIX") {
//             __helpers::add_motif(n->index(), g, *new_g, index_convert);
//             continue;
//         }
//
//         // ideal helix to break up
//         auto num_res = n->data().get_num_residues();
//         auto num_bp_steps = num_res / 2;
//
//         auto start = rm.get_segment(args);
//         auto pos = -1;
//         auto org_pei = 1;
//
//         if(g.has_parent(n->index())) {
//             auto pi = g.get_parent_index(n->index());
//             auto pei = g.get_parent_end_index(n->index());
//             auto new_pi = index_convert[pi];
//             pos = new_g->add_segment(*start, new_pi,
//             new_g->get_segment_end_name(new_pi, pei));
//
//         }
//         else { // no parent but make sure its starting from the correct
//         orientation
//             aligner.align(n->data().get_end(0), *start);
//             pos = new_g->add_segment(*start);
//         }
//
//         for(int i = 0; i < num_bp_steps-2; i++) {
//             auto ideal_bp_step = rm.get_segment(args);
//             pos = new_g->add_segment(*ideal_bp_step, pos,
//             start->get_end_name(1));
//         }
//
//         index_convert[n->index()] = pos;
//
//     }
//
//     for(auto const & n : g) {
//         auto & connections = g.get_motif_connections(n->index());
//         for(auto const & c : connections) {
//             if(c == nullptr) { continue; }
//             auto key = c->to_str();
//             if(seen_connections.find(key) != seen_connections.end()) {
//             continue; }
//
//             auto new_ni = index_convert[c->node_i];
//             auto new_nj = index_convert[c->node_j];
//
//             if(c->node_i == 0) { new_ni = 0; }
//             if(c->node_j == 0) { new_nj = 0; }
//
//             new_g->add_connection(data_structure::NodeIndexandEdge{new_ni,
//             c->edge_i},
//                                   data_structure::NodeIndexandEdge{new_nj,
//                                   c->edge_j});
//
//             seen_connections[key] = 1;
//
//         }
//     }
//
//     return new_g;
// }

template <typename SegmentType, typename AlignerType>
secondary_structure::SegmentGraphOP
get_secondary_structure_graph(SegmentGraph<SegmentType, AlignerType> const &g) {

  if (std::is_same<SegmentType, secondary_structure::Segment>::value) {
    LOGE << "attempting to generate a secondary structure graph from an "
            "existing one";
    throw std::runtime_error("attempting to generate a secondary structure "
                             "graph from an existing one");
  }

  auto ss_g = std::make_shared<secondary_structure::SegmentGraph>();
  auto seen_connections = std::map<String, int>();

  for (auto const &n : g) {
    auto seg = n->data().get_secondary_structure();

    if (g.has_parent(n->index())) {
      auto pi = g.get_parent_index(n->index());
      auto pei = g.get_parent_end_index(n->index());

      auto key1 = data_structure::Edge(n->index(), pi, 0, pei).to_str();
      auto key2 = data_structure::Edge(pi, n->index(), pei, 0).to_str();
      seen_connections[key1] = 1;
      seen_connections[key2] = 1;

      ss_g->add_segment(*seg, pi, ss_g->get_segment_end_name(pi, pei));

    }

    else {
      ss_g->add_segment(*seg);
    }
  }

  auto extra_connections =
      __helpers::get_extra_connections(g, seen_connections);
  for (auto const &e : *extra_connections) {
    ss_g->add_connection(data_structure::NodeIndexandEdge{e.node_i, e.edge_i},
                         data_structure::NodeIndexandEdge{e.node_j, e.edge_j});
  }

  return ss_g;
}

} // namespace segment_data_structure

namespace structure {

// inline
// SegmentGraphOP
// convert_ideal_helices_to_basepair_steps(
//         SegmentGraph const & sg,
//         resources::ResourceManager const & rm) {
//     return
//     segment_data_structure::convert_ideal_helices_to_basepair_steps<Segment,
//     Aligner>(sg, rm);
// }

inline secondary_structure::SegmentGraphOP
get_secondary_structure_graph(SegmentGraph const &sg) {
  return segment_data_structure::get_secondary_structure_graph<Segment,
                                                               Aligner>(sg);
};
     */
} // namespace segment_data_structure

#endif // RNAMAKE_NEW_SEGMENT_GRAPH_H
