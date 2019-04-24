//
// Created by Joseph Yesselman on 2019-04-19.
//

#include <ddg_tecto/ddg_tecto.h>
#include <ddg_tecto/tecto_graph_setup.h>
#include <base/backtrace.hpp>
#include <math/stats.h>
#include <secondary_structure/secondary_structure_parser.h>



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SimulateTectosGraphApp::setup_options() {
    add_option("fseq", "", base::OptionType::STRING, true);
    add_option("fss",  "", base::OptionType::STRING, true);
    add_option("cseq", "", base::OptionType::STRING, true);
    add_option("css",  "", base::OptionType::STRING, true);

    // optional
    add_option("steps", 1000000, base::OptionType::INT);
    add_option("n", 1, base::OptionType::INT);
    add_option("log_level", "info", base::OptionType::STRING);
    add_option("setup", "default", base::OptionType::STRING);
    add_option("sterics", "selective_sterics", base::OptionType::STRING);
    add_option("scoring", "old_frame_scorer", base::OptionType::STRING);
    add_option("logging", "none", base::OptionType::STRING);

}

void
SimulateTectosGraphApp::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);
    _parameters.fseq = get_string_option("fseq");
    _parameters.fss  = get_string_option("fss");
    _parameters.cseq = get_string_option("cseq");
    _parameters.css  = get_string_option("css");
    _parameters.log_level = get_string_option("log_level");
    _parameters.setup = get_string_option("setup");
    _parameters.steps = get_int_option("steps");
    _parameters.n = get_int_option("n");

    // setup logging
    auto log_level = base::log_level_from_str(_parameters.log_level);
    base::init_logging(log_level);

    // report command line options
    LOG_INFO << "log level set to: " << _parameters.log_level;
    LOG_INFO << "flow peice -> sequence: " << _parameters.fseq << " structure: " << _parameters.fss;
    LOG_INFO << "chip peice -> sequence: " << _parameters.cseq << " structure: " << _parameters.css;
}

void
SimulateTectosGraphApp::run() {
    using namespace thermo_fluctuation::graph;

    auto flow_p = _get_ss_pose(_parameters.fseq, _parameters.fss);
    auto chip_p = _get_ss_pose(_parameters.cseq, _parameters.css);
    auto graph_info = _get_tecto_graph_info(_parameters.setup, *flow_p, *chip_p);
    auto mseg_info = _get_mseg(graph_info->mg);
    auto start = _get_updated_node_index_and_edge(graph_info->start, mseg_info->index_hash);
    auto end   = _get_updated_node_index_and_edge(graph_info->end, mseg_info->index_hash);

    auto scorer = std::make_shared<OldFrameScorer>();

    auto steric_node_indices_1 = Ints{end.node_index, end.node_index-1};
    auto steric_node_indices_2 = Ints{start.node_index};

    auto steric_checker = std::make_shared<sterics::SelectiveSterics>(steric_node_indices_1, steric_node_indices_2, 2.2f);
    auto sim = std::make_shared<Simulation>(scorer, steric_checker);

    auto hits = std::vector<double>();
    for(int j = 0; j < _parameters.n; j++) {
        int count = 0;
        sim->setup(*mseg_info->mseg, start, end);
        for (int i = 0; i < _parameters.steps; i++) {
            count += sim->next();
        }
        hits.push_back((double)count);
    }
    LOG_INFO << "RUNS: " << _parameters.n << " AVG: " << math::mean(hits) << " STDEV: " << math::stdev(hits);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

secondary_structure::PoseOP
SimulateTectosGraphApp::_get_ss_pose(
        String const & sequence,
        String const & structure) {
    auto parser = secondary_structure::Parser();
    return parser.parse_to_pose(_convert_sequence_to_RNA(sequence), structure);

}

String
SimulateTectosGraphApp::_convert_sequence_to_RNA(
        String const & sequence) {
    auto seq_rna = String("");
    for(auto const & e : sequence) {
        if(e == 'T') { seq_rna += 'U'; }
        else         { seq_rna += e;   }
    }
    return seq_rna;
}

TectoGraphInfoOP
SimulateTectosGraphApp::_get_tecto_graph_info(
        String const & setup_name,
        secondary_structure::Pose const & flow_p,
        secondary_structure::Pose const & chip_p) {
    auto setup_procedure = TectoGraphSetupOP(nullptr);
    if(setup_name == "default") {
        setup_procedure = TectoGraphSetupOP(std::make_shared<DefaultSetup>()->clone());
    }
    else {
        LOG_ERROR << "invalid setup procedure: " << setup_name; exit(0);
    }
    return setup_procedure->get_graph(flow_p, chip_p);

}

SimulateTectosGraphApp::EnsembleConversionResultsOP
SimulateTectosGraphApp::_get_mseg(
        motif_data_structure::MotifGraphOP mg) {

    auto index_hash = std::map<int, int>();
    auto mseg = std::make_shared<motif_data_structure::MotifStateEnsembleGraph>();
    int i = 0, j = 0;
    for(auto const & n : *mg) {
        // build / get motif state ensemble
        auto mse = motif::MotifStateEnsembleOP(nullptr);

        if(n->data()->mtype() == util::MotifType::HELIX) {
            // is not a basepair step
            if(n->data()->residues().size() > 4) {
                LOG_ERROR << "supplied a helix motif: "+n->data()->name() +" that is not a basepair step, this is no supported";
                exit(0);
            }
            try {
                mse = _rm.motif_state_ensemble(n->data()->end_ids()[0]);
            }
            catch(resources::ResourceManagerException const & e) {
                LOG_ERROR << "cannot find motif state ensemble for basepair with id: " + n->data()->end_ids()[0] <<
                          "check to make sure its a Watson-Crick basepair step";
                exit(0);
            }
        }
        else {
            mse = std::make_shared<motif::MotifStateEnsemble>(n->data()->get_state());
        }
        if(i == 0) {
            j = mseg->add_ensemble(*mse);
        }
        else {
            int pi = index_hash[n->parent_index()];
            j = mseg->add_ensemble(*mse, data_structure::NodeIndexandEdge{pi, n->parent_end_index()});
        }
        index_hash[n->index()] = j;
        i++;
    }
    return std::make_shared<EnsembleConversionResults>(mseg, index_hash);
}

data_structure::NodeIndexandEdge
SimulateTectosGraphApp::_get_updated_node_index_and_edge(
        data_structure::NodeIndexandEdge const & org_nie,
        std::map<int, int> const & index_hash) {
    return data_structure::NodeIndexandEdge{index_hash.at(org_nie.node_index), org_nie.edge_index};
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    //load tectos
    auto tecto_dir = String(base::base_dir()+"/rnamake/lib/RNAMake/apps/simulate_tectos");
    resources::Manager::instance().add_motif(tecto_dir+"/resources/GAAA_tetraloop");
    resources::Manager::instance().add_motif(tecto_dir+"/resources/GGAA_tetraloop");

    auto app = SimulateTectosGraphApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}