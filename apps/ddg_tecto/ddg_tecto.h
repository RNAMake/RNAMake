//
// Created by Joseph Yesselman on 2019-04-19.
//

#ifndef RNAMAKE_NEW_DDG_TECTO_H
#define RNAMAKE_NEW_DDG_TECTO_H

#include <ddg_tecto/tecto_graph_setup.h>
#include <base/application.hpp>
#include <resources/resource_manager.h>
#include <thermo_fluctuation/graph/simulation.h>

class SimulateTectosGraphApp : public base::Application {
public:
    SimulateTectosGraphApp():
        base::Application(),
        _rm(resources::Manager::instance()) {}

    ~SimulateTectosGraphApp() {}

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
    String
    _convert_sequence_to_RNA(
            String const &);

    secondary_structure::PoseOP
    _get_ss_pose(
            String const &,
            String const &);

    TectoGraphInfoOP
    _get_tecto_graph_info(
            String const &,
            secondary_structure::Pose const &,
            secondary_structure::Pose const &);

private:

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

    data_structure::NodeIndexandEdge
    _get_updated_node_index_and_edge(
            data_structure::NodeIndexandEdge const &,
            std::map<int, int> const &);

private:
    struct Parameters {
        String fseq, fss, cseq, css, log_level, setup;
        int steps, n;
    };

    std::map<std::string, std::string> parameters_map = {{"fseq", std::string() },
                             {"fss", std::string()},
                             {"cseq", std::string()},
                             {"css", std::string()},
                             {"log_level", std::string()},
                             {"setup", std::string()}};

private:
    resources::Manager & _rm;
    Parameters _parameters;


};


#endif //RNAMAKE_NEW_DDG_TECTO_H
