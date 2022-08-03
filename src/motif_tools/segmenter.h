//
//  segmenter.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__segmenter__
#define __RNAMake__segmenter__

#include <base/types.hpp>
#include <map>
#include <motif/motif.h>
#include <motif/motif_factory.h>
#include <queue>
#include <stdio.h>
#include <structure/residue.h>

namespace motif_tools {

/*
 * Exception for Segmenter
 */
class SegmenterException : public std::runtime_error {
public:
  /**
   * Standard constructor for SegmenterException
   * @param   message   Error message for segmenter
   */
  SegmenterException(String const &message) : std::runtime_error(message) {}
};

struct Pair {
  Pair(structure::ResidueOP const &nres1, structure::ResidueOP const &nres2,
       int ndist)
      : res1(nres1), res2(nres2), dist(ndist) {}

  int inline contains(structure::ResidueOP const &res) const {

    if (res == res1 || res == res2) {
      return 1;
    } else {
      return 0;
    }
  }

  structure::ResidueOP res1, res2;
  int dist;
};

struct Segments {
  inline Segments(motif::MotifOP const &nremoved,
                  motif::MotifOP const &nremaining)
      : removed(nremoved), remaining(nremaining) {}

  motif::MotifOP removed, remaining;
};

typedef std::shared_ptr<Pair> PairOP;
typedef std::vector<PairOP> PairOPs;
typedef std::shared_ptr<Segments> SegmentsOP;

struct PairSearchNode {
  PairSearchNode(PairOPs const &npairs) : pairs(npairs), score(0) {
    for (auto const &p : pairs) {
      score += p->dist;
    }
  }

  ~PairSearchNode() {}

  inline int contains(structure::ResidueOP const &res) const {

    for (auto const &p : pairs) {
      if (p->contains(res)) {
        return 1;
      }
    }
    return 0;
  }

  inline int contains_pair(PairOP const &pair) {

    for (auto const &p : pairs) {
      if (pair == p) {
        return 1;
      }
    }
    return 0;
  }

  PairOPs pairs;
  int score;
};

typedef std::vector<PairSearchNode> PairSearchNodes;

struct PairSearchNodeCompare {
  bool operator()(PairSearchNode const &node1, PairSearchNode const &node2) {

    if (node1.score > node2.score) {
      return true;
    } else {
      return false;
    }
  }
};

typedef std::priority_queue<PairSearchNode, PairSearchNodes,
                            PairSearchNodeCompare>
    PairSearchNodePriortyQueue;

class PairSearch {
public:
  PairSearch()
      : _queue(PairSearchNodePriortyQueue()), _solutions(PairSearchNodes()),
        _values(std::map<structure::ResidueOP, float>()) {}

  ~PairSearch() {}

public:
  PairSearchNodes const &search(structure::ResidueOPs const &res,
                                PairOPs const &pairs,
                                PairOPs const &end_pairs) {
    _res = res;
    _get_default_values(res, pairs, end_pairs);

    auto all_pairs = PairOPs();
    for (auto const &p : pairs) {
      all_pairs.push_back(p);
    }
    for (auto const &p : end_pairs) {
      all_pairs.push_back(p);
    }

    for (auto const &p : pairs) {
      auto n = PairSearchNode(PairOPs{p});
      n.score += _get_estimated_score(n);
      _queue.push(n);
    }

    auto current = PairSearchNode(PairOPs{});
    while (!_queue.empty()) {
      current = _queue.top();
      _queue.pop();

      for (auto p : all_pairs) {
        if (current.contains(p->res1) || current.contains(p->res2)) {
          continue;
        }
        auto pairs = current.pairs;
        pairs.push_back(p);
        auto node = PairSearchNode(pairs);
        auto estimated_score = _get_estimated_score(node);
        if (estimated_score == 0) {
          _solutions.push_back(node);
          if (_solutions.size() > 10) {
            return _solutions;
          }
        } else {
          node.score += _get_estimated_score(node);
          _queue.push(node);
        }
      }
    }

    return _solutions;
  }

private:
  void _get_default_values(structure::ResidueOPs const &res,
                           PairOPs const &pairs, PairOPs const &end_pairs) {

    auto all_pairs = PairOPs();
    for (auto const &p : pairs) {
      all_pairs.push_back(p);
    }
    for (auto const &p : end_pairs) {
      all_pairs.push_back(p);
    }

    for (auto const &r : res) {
      auto total = 0.0, count = 0.0;
      for (auto const &p : all_pairs) {
        if (p->contains(r)) {
          total += p->dist;
          count += 1;
        }
      }
      _values[r] = (total / count) / 2;
    }
  }

  int _get_estimated_score(PairSearchNode const &n) {
    int score = 0;
    for (auto const &r : _res) {
      if (n.contains(r)) {
        continue;
      }
      score += _values[r];
    }
    return score;
  }

private:
  PairSearchNodePriortyQueue _queue;
  PairSearchNodes _solutions;
  std::map<structure::ResidueOP, float> _values;
  structure::ResidueOPs _res;
};

class Segmenter {
public:
  Segmenter() : mf_(motif::MotifFactory()) {}

  ~Segmenter() {}

public:
  SegmentsOP apply(structure::RNAStructureOP const &,
                   structure::BasepairOPs const &);

private:
  structure::ChainOP _get_subchain(structure::RNAStructureOP const &,
                                   PairOP const &);

  void _get_pairs(structure::RNAStructureOP const &,
                  structure::ResidueOPs const &);

  SegmentsOP _get_segments(structure::RNAStructureOP const &,
                           structure::ResidueOPs &,
                           structure::BasepairOPs const &,
                           structure::ResidueOPs const &,
                           structure::BasepairOPs const &);

private:
  motif::MotifFactory _mf;
  PairOPs _pairs, _end_pairs;
};

} // namespace motif_tools

#endif /* defined(__RNAMake__segmenter__) */
