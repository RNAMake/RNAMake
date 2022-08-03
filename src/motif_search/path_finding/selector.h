//
// Created by Joseph Yesselman on 3/18/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SELECTOR_H
#define RNAMAKE_NEW_PATH_FINDING_SELECTOR_H

#include <data_structure/graph/graph.h>
#include <motif/motif_state.h>
#include <motif/motif_state_ensemble.h>
#include <resources/motif_state_sqlite_library.h>

namespace motif_search {
namespace path_finding {

class SelectorException : public std::runtime_error {
public:
  SelectorException(String const &message) : std::runtime_error(message) {}
};

struct SelectorNodeData {
  inline SelectorNodeData(String const &n_name,
                          motif::MotifStateOPs const &n_motif_states,
                          int n_type)
      : name(n_name), motif_states(n_motif_states), type(n_type) {}

  String name;
  motif::MotifStateOPs motif_states;
  int type;
};

typedef std::shared_ptr<SelectorNodeData> SelectorNodeDataOP;

class Selector {
public:
  Selector() : _pos(0), _max(0), _parent_type(0) {}

  virtual ~Selector() {}

  virtual Selector *clone() const { return new Selector(*this); }

public:
  virtual void add(String const &lib_name) {
    auto type = _graph.size();
    auto motif_states = motif::MotifStateOPs();
    auto ms_lib = resources::MotifStateSqliteLibrary(lib_name);
    ms_lib.load_all();
    for (auto const &ms : ms_lib) {
      motif_states.push_back(ms);
    }
    auto d = std::make_shared<SelectorNodeData>(lib_name, motif_states, type);
    _add(d);
  }

  virtual void add(String const &lib_name, motif::MotifStateEnsembleOP mse) {
    auto type = _graph.size();
    auto motif_states = motif::MotifStateOPs();
    for (auto const &mem : mse->members()) {
      motif_states.push_back(mem->motif_state);
    }
    auto d = std::make_shared<SelectorNodeData>(lib_name, motif_states, type);
    _add(d);
  }

  virtual void add(motif::MotifOP motif) {
    auto type = _graph.size();
    auto motif_states = motif::MotifStateOPs{motif->get_state()};
    auto d =
        std::make_shared<SelectorNodeData>(motif->name(), motif_states, type);
    _add(d);
  }

public:
  void connect(String const &, String const &);

public:
  void start(int parent_type) const {
    _parent_type = parent_type;
    if (parent_type == -1) {
      _max = 0;
    } else {
      _max = _graph.get_node(parent_type)->connections().size() - 1;
    }
    _pos = -1;
  }

  SelectorNodeDataOP next() const {
    if (finished()) {
      throw SelectorException(
          "cannot enumerate anymore selector is done enumerating");
    }
    _pos += 1;
    if (_parent_type == -1) {
      return _graph.get_node(0)->data();
    } else {
      auto c = _graph.get_node(_parent_type)->connections()[_pos];
      return c->partner(_parent_type)->data();
    }
  }

  bool finished() const { return _pos == _max; }

  size_t size() const { return _graph.size(); }

protected:
  void _add(SelectorNodeDataOP d) {
    for (auto const &n : _graph) {
      if (d->name == n->data()->name) {
        throw SelectorException("cannot have two nodes with the same name");
      }
    }
    _graph.add_data(d, -1, 1);
  }

protected:
  data_structure::graph::GraphDynamic<SelectorNodeDataOP> _graph;
  mutable int _pos, _max, _parent_type;
};

class RoundRobinSelector : public Selector {
public:
  RoundRobinSelector() : Selector() {}

  ~RoundRobinSelector() {}

  Selector *clone() const { return new RoundRobinSelector(*this); }

public:
  virtual void add(String const &lib_name) {
    Selector::add(lib_name);
    _add_new_connections();
  }

  virtual void add(String const &lib_name, motif::MotifStateEnsembleOP mse) {
    Selector::add(lib_name, mse);
    _add_new_connections();
  }

  virtual void add(motif::MotifOP motif) {
    Selector::add(motif);
    _add_new_connections();
  }

private:
  void _add_new_connections() {
    auto i = (int)_graph.size() - 1;
    for (int j = 0; j < _graph.size(); j++) {
      if (j == i) {
        continue;
      }
      _graph.connect(j, i);
    }
  }
};

typedef std::shared_ptr<Selector> SelectorOP;

SelectorOP default_selector();

} // namespace path_finding
} // namespace motif_search

#endif // RNAMAKE_NEW_PATH_FINDING_SELECTOR_H
