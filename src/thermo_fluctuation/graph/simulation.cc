//
// Created by Joseph Yesselman on 2019-04-09.
//

#include <thermo_fluctuation/graph/simulation.h>

namespace thermo_fluctuation {
namespace graph {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// option functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulation::setup_options() {
  options_.add_option("temperature", 298.15f, base::OptionType::FLOAT);
  options_.add_option("steric_radius", 2.2f, base::OptionType::FLOAT);
  options_.add_option("cutoff", 4.5f, base::OptionType::FLOAT);

  options_.lock_option_adding();
  update_var_options();
}

void Simulation::update_var_options() {
  parameters_.temperature = get_float_option("temperature");
  parameters_.steric_radius = get_float_option("steric_radius");
  parameters_.cutoff = get_float_option("cutoff");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulation::setup(
    motif_data_structure::MotifStateEnsembleGraph const &mseg,
    data_structure::NodeIndexandEdge const &start,
    data_structure::NodeIndexandEdge const &end) {
  // setup sampler and generate initial conformation
  sampler_ = std::make_shared<Sampler>(mseg);
  sampler_->set_temperature(parameters_.temperature);
  msg_ = sampler_->get_initial_state();

  start_ = start;
  end_ = end;

  bool target_an_aligned_end = false;
  if (msg_->get_node(end_.node_index)->data()->block_end_add() ==
      end_.edge_index) {
    target_an_aligned_end = true;
  }
  scorer_->setup(target_an_aligned_end);
}

bool Simulation::next() {
  auto count = 0;
  auto done = false;
  while (!done) {
    if (sampler_->next(msg_) == 0) {
      count += 1;
      if (count == 100) {
        LOG_ERROR << "sampler is stuck unsure what happened";
        exit(0);
      }
    }
    count = 0;
    done = true;
  }

  if (sterics_->clash(msg_)) {
    return false;
  }

  end_state_1_ = msg_->get_node(start_.node_index)
                     ->data()
                     ->get_end_state(start_.edge_index);
  end_state_2_ =
      msg_->get_node(end_.node_index)->data()->get_end_state(end_.edge_index);
  score_ = scorer_->score(*end_state_1_, *end_state_2_);
  if (score_ < parameters_.cutoff) {
    return true;
  }

  return false;
}

} // namespace graph
} // namespace thermo_fluctuation
