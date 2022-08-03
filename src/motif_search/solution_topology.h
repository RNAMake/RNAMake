//
// Created by Joseph Yesselman on 3/23/19.
//

#ifndef RNAMAKE_NEW_SOLUTION_TOPOLOGY_H
#define RNAMAKE_NEW_SOLUTION_TOPOLOGY_H

#include <data_structure/graph.h>
#include <motif_data_structure/motif_state_ensemble_graph.h>
#include <resources/motif_state_sqlite_library.h>

namespace motif_search {

class SolutionTopologyTemplate {
public:
  enum class NodeType { LIBRARY, MOTIF_STATE, ENSEMBLE };

  class Node {
  public:
    inline Node(String const &lib_name)
        : _type(NodeType::LIBRARY), _lib_name(lib_name) {}

    inline Node(motif::MotifStateOP ms)
        : _type(NodeType::MOTIF_STATE), _motif_state(ms) {}

    inline Node(motif::MotifStateEnsembleOP mse)
        : _type(NodeType::ENSEMBLE), _ensemble(mse) {}

  public:
    inline NodeType get_type() const { return _type; }

    inline String const &get_lib_name() const {
      if (_type != NodeType::LIBRARY) {
        throw std::runtime_error("cannot get_lib_name(), wrong type");
      }
      return _lib_name;
    }

    inline motif::MotifStateOP get_motif_state() const {
      if (_type != NodeType::MOTIF_STATE) {
        throw std::runtime_error("cannot get_motif_state(), wrong type");
      }
      return _motif_state;
    }

    inline motif::MotifStateEnsembleOP get_motif_state_ensemble() const {
      if (_type != NodeType::ENSEMBLE) {
        throw std::runtime_error(
            "cannot get_motif_state_ensemble(), wrong type");
      }
      return _ensemble;
    }

  private:
    NodeType _type;
    String _lib_name;
    motif::MotifStateOP _motif_state;
    motif::MotifStateEnsembleOP _ensemble;
  };

public:
  SolutionTopologyTemplate()
      : _g(data_structure::FixedEdgeDirectedGraph<Node>()) {
    _motif_type_ends = StringIntMap();
    _motif_type_ends["ideal_helices"] = 2;
    _motif_type_ends["flex_helices"] = 2;
    _motif_type_ends["twoway"] = 2;
    _motif_type_ends["unique_twoway"] = 2;
  }

  ~SolutionTopologyTemplate() {}

public: // iterators
  typedef typename data_structure::FixedEdgeDirectedGraph<Node>::const_iterator
      const_iterator;
  typedef
      typename data_structure::FixedEdgeDirectedGraph<Node>::iterator iterator;

  iterator begin() { return _g.begin(); }
  iterator end() { return _g.end(); }

  const_iterator begin() const noexcept { return _g.begin(); }
  const_iterator end() const noexcept { return _g.end(); }

public:
  void add_library(String const &lib_name) {
    auto num_ends = _get_ends_for_motif_type(lib_name);
    auto n = Node(lib_name);
    _g.add_node(n, num_ends);
    _update_default_transveral();
  }

  void add_library(String const &lib_name,
                   data_structure::NodeIndexandEdge const &parent_nie) {
    auto num_ends = _get_ends_for_motif_type(lib_name);
    auto n = Node(lib_name);
    _g.add_node(n, num_ends, 0, parent_nie);
    _update_default_transveral();
  }

public:
  void add_motif_state(motif::MotifStateOP ms) {
    auto n = Node(ms);
    _g.add_node(n, ms->end_names().size());
    _update_default_transveral();
  }

  void add_motif_state(motif::MotifStateOP ms,
                       data_structure::NodeIndexandEdge const &parent_nie) {
    auto n = Node(ms);
    _g.add_node(n, ms->end_names().size(), 0, parent_nie);
    _update_default_transveral();
  }

public:
  void add_ensemble(motif::MotifStateEnsembleOP mse) {
    auto n = Node(mse);
    _g.add_node(n, mse->num_end_states());
    _update_default_transveral();
  }

  void add_ensemble(motif::MotifStateEnsembleOP mse,
                    data_structure::NodeIndexandEdge const &parent_nie) {
    auto n = Node(mse);
    _g.add_node(n, mse->num_end_states(), 0, parent_nie);
    _update_default_transveral();
  }

public:
  inline bool has_parent(Index ni) const { return _g.has_parent(ni); }

  inline Index get_parent_index(Index ni) const {
    return _g.get_parent_index(ni);
  }

  inline Index get_parent_end_index(Index ni) const {
    return _g.get_parent_end_index(ni);
  }

private:
  void _update_default_transveral() {
    auto roots = _g.get_root_indexes();
    if (roots.size() > 0) {
      _g.setup_transversal(roots[0]);
    }
  }

  int _get_ends_for_motif_type(String const &motif_type) {
    if (_motif_type_ends.find(motif_type) == _motif_type_ends.end()) {
      throw std::runtime_error("motif type: " + motif_type +
                               " not accepted in SolutionTopology");
    }
    return _motif_type_ends[motif_type];
  }

private:
  data_structure::FixedEdgeDirectedGraph<Node> _g;
  StringIntMap _motif_type_ends;
};

typedef std::shared_ptr<SolutionTopologyTemplate> SolutionTopologyTemplateOP;

class SolutionToplogy;
typedef std::shared_ptr<SolutionToplogy> SolutionToplogyOP;

class SolutionToplogyFactory {
public:
  SolutionToplogyFactory() {
    setup_options();
    update_var_options();
  }

  ~SolutionToplogyFactory() {}

public:
  SolutionToplogyOP
  generate_toplogy(SolutionTopologyTemplate const &sol_template) {

    auto mseg =
        std::make_shared<motif_data_structure::MotifStateEnsembleOPGraph>();
    for (auto const &n : sol_template) {
      auto mse = motif::MotifStateEnsembleOP();
      if (n->data().get_type() == SolutionTopologyTemplate::NodeType::LIBRARY) {
        auto lib = _get_library(n->data().get_lib_name());
        mse = _parse_library_into_ensemble(lib);
      } else if (n->data().get_type() ==
                 SolutionTopologyTemplate::NodeType::ENSEMBLE) {
        mse = n->data().get_motif_state_ensemble();
      } else if (n->data().get_type() ==
                 SolutionTopologyTemplate::NodeType::MOTIF_STATE) {
        mse = std::make_shared<motif::MotifStateEnsemble>(
            n->data().get_motif_state());
      }

      else {
        throw std::runtime_error("not implemented");
      }

      if (sol_template.has_parent(n->index())) {
        auto parent_nie = data_structure::NodeIndexandEdge{
            sol_template.get_parent_index(n->index()),
            sol_template.get_parent_end_index(n->index())};
        mseg->add_ensemble(mse, parent_nie);
      } else {
        mseg->add_ensemble(mse);
      }
    }

    return std::make_shared<SolutionToplogy>(mseg);
  }

private:
  resources::MotifStateSqliteLibraryOP _get_library(String const &lib_name) {
    if (_libraries.find(lib_name) == _libraries.end()) {
      _libraries[lib_name] =
          std::make_shared<resources::MotifStateSqliteLibrary>(lib_name);
      _libraries[lib_name]->load_all();
    }
    return _libraries[lib_name];
  }

  motif::MotifStateEnsembleOP
  _parse_library_into_ensemble(resources::MotifStateSqliteLibraryOP library) {
    auto motif_states = motif::MotifStateOPs();
    auto energies = Reals();

    auto is_helix_lib = false;
    if (library->get_name() == "flex_helices" ||
        library->get_name() == "ideal_helices" ||
        library->get_name() == "ideal_helices_min") {
      is_helix_lib = true;
    }

    for (auto const &ms : *library) {
      if (is_helix_lib && (ms->size() > _parameters.max_helix_size ||
                           ms->size() < _parameters.min_helix_size)) {
        continue;
      }
      ms->new_uuids();
      motif_states.push_back(std::make_shared<motif::MotifState>(*ms));
      energies.push_back(1);
    }
    if (motif_states.size() == 0) {
      LOG_ERROR << library->get_name() << " has no motifs! ";
      if (is_helix_lib) {
        LOG_ERROR << " this is a helix library did you set the min and max "
                     "helix incorrectly??";
      }
      exit(0);
    }
    return std::make_shared<motif::MotifStateEnsemble>(motif_states, energies);
  }

private:
  struct Parameters {
    int max_helix_size, min_helix_size;
  };

  void setup_options() {
    _options.add_option("max_helix_size", 99, base::OptionType::INT);
    _options.add_option("min_helix_size", 6, base::OptionType::INT);
    _options.lock_option_adding();
  }

  void update_var_options() {
    _parameters.max_helix_size = _options.get_int("max_helix_size");
    _parameters.min_helix_size = _options.get_int("min_helix_size");
  }

public: // option wrappers
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
  base::Options _options;
  Parameters _parameters;
  std::map<String, resources::MotifStateSqliteLibraryOP> _libraries;
  std::map<String, motif::MotifStateEnsembleOP> _ensembles;
};

class SolutionToplogy {
public:
  SolutionToplogy(motif_data_structure::MotifStateEnsembleOPGraphOP mseg) {
    _mseg = mseg;
    _rng = util::RandomNumberGenerator();
    _solution_nie = _mseg->get_leafs();
  }

public:
  typedef
      typename motif_data_structure::MotifStateEnsembleOPGraph::const_iterator
          const_iterator;
  typedef typename motif_data_structure::MotifStateEnsembleOPGraph::iterator
      iterator;

  iterator begin() { return _mseg->begin(); }
  iterator end() { return _mseg->end(); }

  const_iterator begin() const noexcept { return _mseg->begin(); }
  const_iterator end() const noexcept { return _mseg->end(); }

public:
  // TODO come up with better system to distingiush these too options
  // This one is for Path finding
  motif_data_structure::MotifStateGraphOP
  initialize_solution(structure::BasepairStateOP bp_state) {

    auto ms = std::make_shared<motif::MotifState>(
        "start", Strings{"start", "start"}, Strings{"", ""},
        structure::BasepairStateOPs{bp_state, bp_state}, math::Vector3s(), 0, 0,
        0);

    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);
    msg->add_state(ms);

    for (auto const &n : *_mseg) {
      auto ms = get_motif_state(n->index());
      if (_mseg->has_parent(n->index())) {
        msg->add_state(ms, _mseg->get_parent_index(n->index()) + 1,
                       _mseg->get_parent_end_index(n->index()));
      } else {
        msg->add_state(ms);
      }
    }
    return msg;
  }

  // this one is for MC
  motif_data_structure::MotifStateGraphOP
  initialize_solution_no_start(structure::BasepairStateOP bp_state) {

    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);
    // msg->add_state(ms);

    for (auto const &n : *_mseg) {
      auto ms = get_motif_state(n->index());
      if (_mseg->has_parent(n->index())) {
        msg->add_state(ms, _mseg->get_parent_index(n->index()),
                       _mseg->get_parent_end_index(n->index()));
      } else {
        msg->add_state(ms);
      }
    }
    return msg;
  }

  motif::MotifStateOP get_motif_state(Index pos) {
    _max_member = _mseg->get_ensemble(pos)->size();
    _mem_pos = _rng.randrange(_max_member);
    return _mseg->get_ensemble(pos)->get_member(_mem_pos)->motif_state;
  }

  std::vector<data_structure::NodeIndexandEdge> const &get_solution_nie() {
    return _solution_nie;
  }

  inline size_t size() { return _mseg->size(); }

  inline size_t get_ensemble_size(Index pos) {
    return _mseg->get_ensemble(pos)->size();
  }

private:
  motif_data_structure::MotifStateEnsembleOPGraphOP _mseg;
  util::RandomNumberGenerator _rng;
  int _max_member, _mem_pos;
  std::vector<data_structure::NodeIndexandEdge> _solution_nie;
};

} // namespace motif_search

#endif // RNAMAKE_NEW_SOLUTION_TOPOLOGY_H
