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
        String mg_file, log_level;
        int steps, n;
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
    
    
private:
    resources::Manager & rm_;
    Parameters parameters_;
        
};



#endif /* thermo_simulation_hpp */
