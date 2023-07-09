//
// Created by Joseph Yesselman on 4/25/18.
//

#ifndef TEST_THERMO_FLUC_GRAPH_SAMPLER_H
#define TEST_THERMO_FLUC_GRAPH_SAMPLER_H

#include <util/monte_carlo.h>
#include <motif_data_structure/motif_state_ensemble_graph.h>

namespace thermo_fluctuation {
namespace graph {

class Sampler {
public:
    Sampler(
            motif_data_structure::MotifStateEnsembleGraph const & mseg) :
            mseg_(mseg) {
        set_temperature(298.15f);
    }

public:
    motif_data_structure::MotifStateGraphOP
    get_initial_state();

    int
    next(
            motif_data_structure::MotifStateGraphOP);

    void
    undo();

public:

    inline
    void
    set_temperature(
            float temp) {
        //set MonteCarlo temperature to kBT, Boltzmann constant in pN.A/K
        mc_.set_temperature(temp * 1.3806488e-1);

    }

private:
    void
    _update(
            motif_data_structure::MotifStateGraphOP);

private:
    util::MonteCarlo mc_;
    util::RandomNumberGenerator rng_;
    motif_data_structure::MotifStateEnsembleGraph mseg_;
    motif::MotifStateEnsembleMemberOP member_;
    Ints max_sizes_, indexes_, states_;
    Floats energies_;
    int node_num_, mem_pos_, accept_;
    //hold last move to undo
    int last_state_pos_, last_num_;
    motif::MotifStateOP last_state_;
};

typedef std::shared_ptr<Sampler> SamplerOP;

}
}

#endif //TEST_THERMO_FLUC_GRAPH_SAMPLER_H
