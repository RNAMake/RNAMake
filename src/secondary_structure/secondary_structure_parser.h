//
//  secondary_structure_parser.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_parser__
#define __RNAMake__secondary_structure_parser__

#include <cassert>
#include <map>
#include <stdio.h>

// RNAMake
#include "base/exception.hpp"
#include "data_structure/graph/graph.h"
#include "secondary_structure/basepair.h"
#include "secondary_structure/motif.h"
#include "secondary_structure/pose.h"
#include "secondary_structure/residue.h"
#include "secondary_structure/structure.h"

namespace secondary_structure {

enum NodeType { UNPAIRED = 0, PAIRED = 1 };

struct NodeData {
  NodeData() {}

  NodeData(ResidueOPs const &nresidues, NodeType const &ntype)
      : residues(nresidues), type(ntype) {}

  ResidueOPs residues;
  NodeType type;
};

class SecondaryStructureChainGraph {
public:
  SecondaryStructureChainGraph()
      : _graph(data_structure::graph::GraphStatic<NodeData>()) {}

  ~SecondaryStructureChainGraph() {}

public:
  typedef
      typename data_structure::graph::GraphStatic<NodeData>::iterator iterator;
  typedef typename data_structure::graph::GraphStatic<NodeData>::const_iterator
      const_iterator;

  iterator begin() { return _graph.begin(); }
  iterator end() { return _graph.end(); }

  const_iterator begin() const { return _graph.begin(); }
  const_iterator end() const { return _graph.end(); }

public:
  size_t size() { return _graph.size(); }

  data_structure::graph::GraphNodeOPs<NodeData> const &nodes() {
    return _graph.nodes();
  }

public:
  int add_chain(NodeData const &data, int parent_index = -1, int orphan = 0) {

    auto parent = _graph.last_node();
    if (parent_index != -1) {
      parent = _graph.get_node(parent_index);
    }
    if (parent == nullptr) {
      return _graph.add_data(data, -1, -1, -1, 3);
    }

    return _graph.add_data(data, parent_index, 1, 0, 3, orphan);
  }

  int get_node_by_res(ResidueOP const &res) {
    for (auto const &n : _graph) {
      for (auto const &r : n->data().residues) {
        if (r->uuid() == res->uuid()) {
          return n->index();
        }
      }
    }
    return -1;
  }

  void pair_res(int n_i, int n_j) {

    assert(graph_.get_node(n_i)->data().type == NodeType::PAIRED &&
           "unpaired node is being paired");
    assert(graph_.get_node(n_j)->data().type == NodeType::PAIRED &&
           "unpaired node is being paired");
    _graph.connect(n_i, n_j, 2, 2);
  }

private:
  data_structure::graph::GraphStatic<NodeData> _graph;
};

typedef std::shared_ptr<data_structure::graph::GraphNode<NodeData>> SSNodeOP;
typedef std::shared_ptr<SecondaryStructureChainGraph>
    SecondaryStructureChainGraphOP;

class Parser {
public:
  Parser()
      : _seen(std::map<SSNodeOP, int>()), _structure(StructureOP()),
        _residues(ResidueOPs()), _pairs(BasepairOPs()) {}

  ~Parser() {}

public:
  SecondaryStructureChainGraphOP parse(String const &, String const &);

  MotifOPs parse_to_motifs(String const &, String const &);

  MotifOP parse_to_motif(String const &, String const &);

  PoseOP parse_to_pose(String const &, String const &);

  void reset() {
    _structure = StructureOP();
    _residues = ResidueOPs();
    _pairs = BasepairOPs();
  }

private:
  void _add_unpaired_residues_to_graph(SecondaryStructureChainGraphOP &,
                                       ResidueOPs const &, int);

  void _add_paired_res_to_graph(SecondaryStructureChainGraphOP &,
                                ResidueOP const &, int);

  BasepairOP _get_previous_pair(ResidueOP const &);

  ResidueOP _previous_res(ResidueOP const &r) {

    int i = (int)(std::find(_residues.begin(), _residues.end(), r) -
                  _residues.begin());
    if (i == 0) {
      return nullptr;
    } else {
      return _residues[i - 1];
    }
  }

  int _start_of_chain(ResidueOP const &r) {

    for (auto const &c : _structure->chains()) {
      if (c->first() == r) {
        return 1;
      }
    }
    return 0;
  }

  ResidueOP _get_bracket_pair(ResidueOP const &r_start) {

    int bracket_count = 0;
    int start = 0;
    for (auto const &r : _residues) {
      if (r_start == r && !start) {
        start = 1;
      } else if (start) {
      } else {
        continue;
      }

      if (r->dot_bracket() == "(") {
        bracket_count += 1;
      } else if (r->dot_bracket() == ")") {
        bracket_count -= 1;
        if (bracket_count == 0) {
          return r;
        }
      }
    }

    throw Exception("cannot find pair in _get_bracket_pair");
  }

  MotifOP _generate_motif(SSNodeOP const &);

  SSNodeOP _walk_nodes(SSNodeOP const &);

  MotifOP _build_motif(StructureOP const &);

private:
  MotifOPs _parse_to_motifs(SecondaryStructureChainGraphOP);

private:
  BasepairOPs _pairs;
  ResidueOPs _residues;
  StructureOP _structure;
  std::map<SSNodeOP, int> _seen;
  ChainOP _chain;
};

} // namespace secondary_structure
#endif /* defined(__RNAMake__secondary_structure_parser__) */
