//
// Created by Joseph Yesselman on 2019-04-21.
//

#ifndef RNAMAKE_NEW_THERMO_FLUC_STERICS_H
#define RNAMAKE_NEW_THERMO_FLUC_STERICS_H

#include <motif_data_structure/motif_state_graph.hpp>

namespace thermo_fluctuation {
namespace graph {
namespace sterics {

class Sterics {
public:
  Sterics() {}

  virtual ~Sterics() {}

  virtual Sterics *clone() const = 0;

public:
  virtual bool clash(motif_data_structure::MotifStateGraphOP) = 0;
};

typedef std::shared_ptr<Sterics> StericsOP;

class NoSterics : public Sterics {
public:
  NoSterics() {}

  ~NoSterics() {}

  Sterics *clone() const { return new NoSterics(*this); }

public:
  bool clash(motif_data_structure::MotifStateGraphOP msg) { return false; }
};

class SelectiveSterics : public Sterics {
public:
  SelectiveSterics(Indexes const node_indices_1, Indexes const node_indices_2,
                   float steric_radius)
      : node_indices_1_(node_indices_1), node_indices_2_(node_indices_2),
        steric_radius_(steric_radius) {}

  ~SelectiveSterics() {}

  Sterics *clone() const { return new SelectiveSterics(*this); }

public:
  bool clash(motif_data_structure::MotifStateGraphOP msg) {
    for (auto const &i : node_indices_1_) {
      for (auto const &j : node_indices_2_) {
        for (auto const &b2 : msg->get_node(i)->data()->cur_state->beads()) {
          for (auto const &b1 : msg->get_node(j)->data()->cur_state->beads()) {
            if (b1.distance(b2) < steric_radius_) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }

private:
  Indexes node_indices_1_, node_indices_2_;
  float steric_radius_;
};

} // namespace sterics
} // namespace graph
} // namespace thermo_fluctuation

#endif // RNAMAKE_NEW_THERMO_FLUC_STERICS_H
