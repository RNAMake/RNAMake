//
// Created by Joseph Yesselman on 4/25/18.
//

#include <thermo_fluctuation/graph/sampler.h>

namespace thermo_fluctuation {
namespace graph {

motif_data_structure::MotifStateGraphOP
Sampler::get_initial_state() {
    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);

    for (auto const & n : mseg_) {
        max_sizes_.push_back(n->data().size());
        mem_pos_ = rng_.randrange(max_sizes_.back());
        if (n->data().size() != 1) {
            indexes_.push_back(n->index());
        }
        auto ms = n->data().get_member(mem_pos_)->motif_state;
        ms->new_uuids();
        energies_.push_back(n->data().get_member(mem_pos_)->energy);
        states_.push_back(mem_pos_);
        if (mseg_.has_parent(n->index())) {
            msg->add_state(ms, mseg_.get_parent_index(n->index()), mseg_.get_parent_end_index(n->index()));
        } else {
            msg->add_state(ms);
        }
    }
    return msg;
}


int
Sampler::next(
        motif_data_structure::MotifStateGraphOP msg) {
    // select index from nodes that have more than 1 member
    node_num_ = indexes_[rng_.randrange(indexes_.size())];
    mem_pos_ = rng_.randrange(max_sizes_[node_num_]);

    member_ = mseg_.get_ensemble(node_num_).get_member(mem_pos_);

    if (mc_.accept(energies_[node_num_], member_->energy)) {
        _update(msg);
        return 1;
    } else {
        return 0;
    }

}

void
Sampler::_update(
        motif_data_structure::MotifStateGraphOP msg) {
    last_state_pos_ = states_[node_num_];
    states_[node_num_] = mem_pos_;
    energies_[node_num_] = member_->energy;
    last_state_ = msg->get_node(node_num_)->data()->ref_state;
    last_num_ = node_num_;
    msg->replace_state(node_num_, member_->motif_state);
}

}
}