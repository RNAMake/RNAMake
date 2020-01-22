//
//  thermo_simulation.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef thermo_simulation_hpp
#define thermo_simulation_hpp

#include <stdio.h>

#include <stdio.h>

#include <base/application.hpp>
#include <data_structure/graph.h>
#include <motif_data_structure/motif_graph.h>
#include <motif_data_structure/motif_state_ensemble_graph.h>
#include <resources/resource_manager.h>
#include <thermo_fluctuation/graph/simulation.h>

struct ConnectionTemplate {
    data_structure::NodeIndexandEdge start;
    data_structure::NodeIndexandEdge end;

    inline
    bool
    operator == (
            ConnectionTemplate const & c) const {
        return (c.start == start && c.end == end);
    }
};

class ThermoSimulationApp : public base::Application {
public:
    struct Parameters {
        String mg_file, log_level, extra_sequences, score_file;
        int steps, n;
        bool movie;
    };


public:
    ThermoSimulationApp():
            base::Application(),
            rm_(resources::Manager::instance()){}
    
    ~ThermoSimulationApp() {}
    
public:
    
    void
    setup_options();
    
    void
    parse_command_line(
        int,
        const char **);
    
    void
    run();

private:

    float
    _calc_initial_score(
            thermo_fluctuation::graph::ScorerOP,
            motif_data_structure::MotifGraphOP,
            data_structure::NodeIndexandEdge const &,
            data_structure::NodeIndexandEdge const &);

    int
    _calc_num_hits(
            thermo_fluctuation::graph::SimulationOP,
            motif_data_structure::MotifStateEnsembleGraph const &,
            data_structure::NodeIndexandEdge const &,
            data_structure::NodeIndexandEdge const &,
            int);

    ConnectionTemplate
    _guess_connection(
            motif_data_structure::MotifGraphOP);


    struct EnsembleConversionResults {
        inline
        EnsembleConversionResults(
                motif_data_structure::MotifStateEnsembleGraphOP n_mseg,
                std::map<int, int> const & n_index_hash):
                mseg(n_mseg),
                index_hash(n_index_hash) {}

        motif_data_structure::MotifStateEnsembleGraphOP mseg;
        std::map<int, int> index_hash;
    };

    typedef std::shared_ptr<EnsembleConversionResults> EnsembleConversionResultsOP;

    EnsembleConversionResultsOP
    _get_mseg(
            motif_data_structure::MotifGraphOP);

    void
    _parse_extra_sequences(
            String const &);
    
    
private:
    resources::Manager & rm_;
    Parameters parameters_;
    std::map<int, Strings> extra_sequences_;
        
};



#endif /* thermo_simulation_hpp */
