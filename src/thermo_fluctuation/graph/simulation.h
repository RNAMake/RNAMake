//
// Created by Joseph Yesselman on 2019-04-09.
//

#ifndef RNAMAKE_NEW_THERMO_FLUC_GRAPH_SIMULATION_H
#define RNAMAKE_NEW_THERMO_FLUC_GRAPH_SIMULATION_H

#include <thermo_fluctuation/graph/sampler.h>
#include <thermo_fluctuation/graph/scorer.h>
#include <thermo_fluctuation/graph/sterics.h>

namespace thermo_fluctuation {
namespace graph {


class Simulation {
public:
    Simulation(
            ScorerOP scorer,
            sterics::StericsOP sterics_):
            scorer_(scorer->clone()),
            sterics_(sterics_->clone()) {
        options_ = base::Options();
        setup_options();
    }

    ~Simulation() {}

public:

    void
    setup(
            motif_data_structure::MotifStateEnsembleGraph const &,
            data_structure::NodeIndexandEdge const &,
            data_structure::NodeIndexandEdge const &);

    bool
    next();

public: // setters


public: // outputs
    void
    write_pdbs(
            String const & name = "nodes") {
        msg_->to_motif_graph()->write_pdbs(name);
    }

    String
    get_pdb_str() {
        return msg_->to_motif_graph()->pdb_str();
    }


public: //option wrappers

    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }

    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }

    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }

    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }


    template<typename T>
    void
    set_option_value(
            String const & name,
            T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }

protected:
    void
    setup_options();

    void
    update_var_options();


private:
    struct Parameters {
        float temperature, steric_radius, cutoff;
    };

private:
    motif_data_structure::MotifStateGraphOP msg_;
    data_structure::NodeIndexandEdge start_, end_;
    structure::BasepairStateOP end_state_1_, end_state_2_;
    ScorerOP scorer_;
    sterics::StericsOP sterics_;
    SamplerOP sampler_;
    Parameters parameters_;
    base::Options options_;
};

typedef std::shared_ptr<Simulation> SimulationOP;

}
}

#endif //RNAMAKE_NEW_THERMO_FLUC_GRAPH_SIMULATION_H