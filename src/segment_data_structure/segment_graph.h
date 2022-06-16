//
// Created by Joseph Yesselman on 1/14/18.
//

#ifndef RNAMAKE_NEW_SEGMENT_GRAPH_H
#define RNAMAKE_NEW_SEGMENT_GRAPH_H

#include <data_structure/graph.h>
#include <structure/aligner.h>
//#include <resources/resource_manager.h>
#include <secondary_structure/aligner.h>
#include <secondary_structure/segment.h>

namespace segment_data_structure {

template <typename SegmentType, typename AlignerType> class SegmentGraph {
public:
  SegmentGraph()
      : aligner_(AlignerType()),
        graph_(data_structure::FixedEdgeDirectedGraph<SegmentType>()),
        needs_update_(false), first_update_(INT_MAX) {}

  SegmentGraph(SegmentGraph const &sg)
      : aligner_(AlignerType()), graph_(sg.graph_), needs_update_(false),
        first_update_(INT_MAX) {
    _update_default_transveral();
  }

public:
  typedef typename data_structure::FixedEdgeDirectedGraph<
      SegmentType>::const_iterator const_iterator;
  typedef typename data_structure::FixedEdgeDirectedGraph<SegmentType>::iterator
      iterator;

  iterator begin() {
    if (needs_update_) {
      _update_alignments(first_update_);
      needs_update_ = false;
      first_update_ = INT_MAX;
    }
    return graph_.begin();
  }
  iterator end() { return graph_.end(); }

  const_iterator begin() const noexcept {
    if (needs_update_) {
      _update_alignments(first_update_);
      needs_update_ = false;
      first_update_ = INT_MAX;
    }

    return graph_.begin();
  }
  const_iterator end() const noexcept { return graph_.end(); }

public:
  inline size_t get_num_segments() const { return graph_.get_num_nodes(); }

  inline SegmentType const &get_segment(Index ni) const {

    if (needs_update_) {
      _update_alignments(first_update_);
      needs_update_ = false;
      first_update_ = INT_MAX;
    }

    return graph_.get_node_data(ni);
  }

  inline SegmentType &get_segment(Index ni) {
    if (needs_update_) {
      _update_alignments(first_update_);
      needs_update_ = false;
      first_update_ = INT_MAX;
    }
    return graph_.get_node_data(ni);
  }

  inline std::vector<data_structure::Edge const *> const &
  get_segment_connections(Index ni) const {
    return graph_.get_node_edges(ni);
  }

public:
  inline bool has_parent(Index ni) const { return graph_.has_parent(ni); }

  inline Index get_parent_index(Index ni) const {
    return graph_.get_parent_index(ni);
  }

  inline Index get_parent_end_index(Index ni) const {
    return graph_.get_parent_end_index(ni);
  }

  inline bool are_motifs_connected(Index n1, Index n2) const {
    return graph_.edge_between_nodes(n1, n2);
  }

  inline std::vector<data_structure::Edge const *> const &
  get_motif_connections(Index ni) const {
    return graph_.get_node_edges(ni);
  }

public:
  Index add_segment(SegmentType const &seg) {
    auto ni = graph_.add_node(seg, seg.get_num_ends());
    _update_default_transveral();
    return ni;
  }

  Index add_segment(SegmentType const &seg, Index parent_index,
                    String const &parent_end_name) {
    auto seg_copy = seg;
    auto &parent = graph_.get_node_data(parent_index);
    auto parent_end_index = parent.get_end_index(parent_end_name);
    aligner_.align(parent.get_end(parent_end_index), seg_copy);
    auto ni = graph_.add_node(
        seg_copy, seg.get_num_ends(), 0,
        data_structure::NodeIndexandEdge{parent_index, parent_end_index});
    _update_default_transveral();
    return ni;
  }

public:
  void add_connection(data_structure::NodeIndexandEdge const &nie1,
                      data_structure::NodeIndexandEdge const &nie2) {
    return graph_.add_edge(nie1, nie2);
  }

public:
  void remove_segment(Index ni) {
    graph_.remove_node(ni);
    auto roots = graph_.get_root_indexes();
    if (roots.size() > 0) {
      graph_.setup_transversal(roots[0]);
    }
  }

  void replace_segment(Index pos, SegmentType const &seg, bool copy = true) {

    if (!copy) {
      get_segment(pos) = seg;
    } else {
      get_segment(pos) = SegmentType(seg);
    }

    needs_update_ = true;
    if (first_update_ > pos) {
      first_update_ = pos;
    }
  }

public:
  String get_segment_end_name(Index n_index, Index end_index) const {
    auto temp = graph_.get_node_data(n_index).get_end(end_index).get_name_str();

    std::cout << "Inside get_segment_end_nake: " << temp << std::endl;

    return temp;
  }

public:
  void write_nodes_to_pdbs(String const &name) const {
    if (needs_update_) {
      _update_alignments(first_update_);
      needs_update_ = false;
      first_update_ = INT_MAX;
    }

    for (auto const &n : graph_) {
      n->data().write_pdb(name + "." + std::to_string(n->index()) + ".pdb");
    }
  }

protected:
  void _update_default_transveral() {
    auto roots = graph_.get_root_indexes();
    if (roots.size() > 0) {
      graph_.setup_transversal(roots[0]);
    }
  }

  void _update_alignments(int start) const {

    auto parent_index = 0;
    auto parent_end_index = 0;
    auto needs_update = false;
    for (auto &n : graph_) {
      if (!has_parent(n->index())) {
        needs_update = false;
        continue;
      }

      if (start == n->index()) {
        needs_update = true;
      }

      if (!needs_update) {
        continue;
      }

      parent_index = get_parent_index(n->index());
      parent_end_index = get_parent_end_index(n->index());

      aligner_.align(
          graph_.get_node_data(parent_index).get_end(parent_end_index),
          n->data());
    }
  }

protected:
  AlignerType aligner_;
  // lots of mutatables to allow lazy updates
  mutable data_structure::FixedEdgeDirectedGraph<SegmentType> graph_;
  mutable bool needs_update_;
  mutable int first_update_;
};

template <typename SegmentType, typename AlignerType>
using SegmentGraphOP = std::shared_ptr<SegmentGraph<SegmentType, AlignerType>>;

} // namespace segment_data_structure

// typedefs
//////////////////////////////////////////////////////////////////////////////////

namespace structure {

typedef segment_data_structure::SegmentGraph<Segment, Aligner> SegmentGraph;
typedef std::shared_ptr<SegmentGraph> SegmentGraphOP;

} // namespace structure

namespace secondary_structure {

typedef segment_data_structure::SegmentGraph<Segment, Aligner> SegmentGraph;
typedef std::shared_ptr<SegmentGraph> SegmentGraphOP;

} // namespace secondary_structure

// conversion functions
/////////////////////////////////////////////////////////////////////////////////

namespace segment_data_structure {

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

} // namespace structure

#endif // RNAMAKE_NEW_SEGMENT_GRAPH_H
