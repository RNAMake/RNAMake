// Made by Chris Jurich 08/25/2020

#ifndef __BUILD_SQLITE_LIBRARIES_H__
#define __BUILD_SQLITE_LIBRARIES_H__

#include <set>
#include <unordered_set>
#include <map>
#include <dirent.h>
#include <sstream>
#include <functional>
#include <filesystem>

#include <base/file_io.h>
#include <base/env_manager.h>
#include <base/log.h>
#include <base/types.h>
#include <base/settings.h>
#include <base/application.hpp>
#include <base/global_constants.h>
#include <base/backtrace.h>
#include <base/exception.h>

#include <resources/motif_sqlite_library.h>
#include <resources/sqlite_library.h>
#include <motif/motif_factory.h>
#include <structure/basepair_state.h>
#include <motif_data_structure/motif_tree.h>

// external
#include <CLI/CLI.hpp>
#include <sqlite_modern/sqlite_modern_cpp.h>

enum DIRECTORIES {
    MOTIF_LIBRARIES,
    MOTIFS,
    MOTIF_STATE_LIBRARIES,
    MOTIF_ENSEMBLE_LIBRARIES
};

// this structure represents a build option
struct BuildOpt {
    std::function<void()> handle;
    bool selected = false;
    std::vector<BuildOpt*> children = {nullptr};
    std::vector<DIRECTORIES> output_dirs;

    void
    select() {
        selected = true;
        for(auto& child : children) {
            if (child != nullptr) {
                child->select();
            }
        }
    }
};

class BuildSqliteLibraries : public base::Application {

public:
    BuildSqliteLibraries() :
    lib_names_(resources::MotifSqliteLibrary::get_libnames()),
    app_("BuildSqliteLibraries"),
    build_opts_({
            // working methods
            {"ideal_helices",                  BuildOpt{[this] { _build_ideal_helices(); },false,{nullptr},{DIRECTORIES::MOTIFS,DIRECTORIES::MOTIF_LIBRARIES}}},
            {"basic_libraries",                BuildOpt{[this] { _build_basic_libraries(); },false,{nullptr},{DIRECTORIES::MOTIF_LIBRARIES}}},
            {"helix_ensembles",                BuildOpt{[this] { _build_helix_ensembles(); },false,{nullptr},{DIRECTORIES::MOTIF_LIBRARIES,DIRECTORIES::MOTIF_ENSEMBLE_LIBRARIES}}},
            {"new_bp_steps",                   BuildOpt{[this] { _build_new_bp_steps(); },false,{nullptr},{DIRECTORIES::MOTIF_LIBRARIES}}},
            {"motif_state_libraries",          BuildOpt{[this] { _build_motif_state_libraries(); },false,{nullptr},{DIRECTORIES::MOTIF_STATE_LIBRARIES}}},
            {"unique_twoway_library",          BuildOpt{[this] { _build_unique_twoway_library(); },false,{nullptr},{DIRECTORIES::MOTIF_LIBRARIES}}},
            {"ss_and_seq_libraries",           BuildOpt{[this] { _build_ss_and_seq_libraries(); },false,{nullptr},{DIRECTORIES::MOTIF_ENSEMBLE_LIBRARIES}}},

                // methods that don't work or don't have working python versions
                //{"trimmed_ideal_helix_library",    BuildOpt{[this] { _build_trimmed_ideal_helix_library();},false,{nullptr}}}
                //{"existing_motif_library",         BuildOpt{[this] { build_existing_motif_library(); },false}},
                //{"flex_helix_library",             BuildOpt{[this] { build_flex_helix_library(); },false}},
             {"le_helix_libraries",                   BuildOpt{[this] { _build_le_helix_lib(); },false,{nullptr},{DIRECTORIES::MOTIF_LIBRARIES}}},
                //{"motif_ensemble_state_libraries", BuildOpt{[this] { _build_motif_ensemble_state_libraries(); },false}},
    })
    {

    }
    //_build_ideal_helices ->
    // _build_basic_libraries & _build_flex_helix_library ->
    // _build_helix_ensembles & _build_new_bp_steps & _build_unique_twoway_library
    // -> _build_motif_state_libraries & _build_motif_ensemble_state_libraries
    //

public:

public:
    void
    run() override ;

public:
    void
    setup_options() override ;

public: // librrary building methods

    void
	_build_ideal_helices();

    void
	_build_basic_libraries();

    void
	_build_existing_motif_library();

    void
	_build_flex_helix_library();

    void
	_build_helix_ensembles();

    void
	_build_new_bp_steps();

    void
	_build_le_helix_lib();

    void
	_build_motif_state_libraries();

    void
	_build_unique_twoway_library();

    void
	_build_ss_and_seq_libraries();

    void
	_build_motif_ensemble_state_libraries();

    void
	_build_trimmed_ideal_helix_library();

public: // public helper methods
    void
    _directory_setup() {

        if(!std::filesystem::exists(parameters_.output_dir)) {
            std::filesystem::create_directory(parameters_.output_dir);
        } else {
            LOGW<<"WARNING: The directory "<<parameters_.output_dir<<
            " already exists. The contained files WILL be overwritten.";
        }

        // should I just delete all files here?
        for(auto& required : required_dirs_) {
            String path = parameters_.output_dir;
            switch (required) {
                case MOTIF_LIBRARIES: {
                    path += "/motif_libraries_new/";
                    break;
                }
                case MOTIFS: {
                    path += "/motifs/";
                    break;
                }
                case MOTIF_STATE_LIBRARIES: {
                    path += "/motif_state_libraries_new/";
                    break;
                }
                case MOTIF_ENSEMBLE_LIBRARIES: {
                    path += "/motif_ensemble_libraries_new/";
                    break;
                }
            }
            if(std::filesystem::exists(path)) {
                LOGW<<"WARNING: The directory "<<path<<" is not empty. Files WILL be overwritten.";
            } else {
                std::filesystem::create_directory(path);
            }
        }
    }

    base::LogLevel
    log_level() const {
        return base::log_level_from_str(parameters_.log_level);
    }

    void
    _generate_method_dependency_tree() {
        // gotta build that tree!!
        // flex_hleix_library and motif_
        build_opts_.at("ideal_helices").children = {&build_opts_.at("basic_libraries"),
                //&build_opts_.at("flex_helix_library"),
                                                    &build_opts_.at("ss_and_seq_libraries")
        };
        build_opts_.at("basic_libraries").children = {&build_opts_.at("helix_ensembles"),
                                                      &build_opts_.at("new_bp_steps"),
                                                      &build_opts_.at("unique_twoway_library")};
        //build_opts_["flex_helix_library"].children = build_opts_["basic_libraries"].children;
        build_opts_.at("helix_ensembles").children = {&build_opts_.at("motif_state_libraries"),
                                                      &build_opts_.at("le_helix_libraries")
                                                      //&build_opts_["motif_ensemble_state_libraries"]
        };
        build_opts_.at("new_bp_steps").children = build_opts_.at("helix_ensembles").children;
        build_opts_.at("unique_twoway_library").children = build_opts_.at("helix_ensembles").children;
    }

public: // public group members
    CLI::App app_;

private:
    StringStringMap
    _get_valid_dirs(String const&);

private:
    struct Parameters {
        bool  all = false;
        String input_dir;
        String output_dir;
        String log_level;
    };

    StringStringMap lib_names_;
    std::map<String,motif::MotifOP> motif_map_;
    std::map<String,BuildOpt> build_opts_;
    Parameters parameters_;
    std::set<DIRECTORIES> required_dirs_;
    const Strings motif_keys_{"data","name","end_name","end_id","id"};
    const Strings ensemble_keys_{"data","name","id"};
    const Strings method_order_{
                "ss_and_seq_libraries",
                "motif_state_libraries",
    //            "motif_ensemlble_state_libraries", // these are commented out for now... they will eventually be re-implemented CJ 09/20
                "le_helix_libraries",
                "helix_ensembles",
                "new_bp_steps",
                "unique_twoway_library",
                "basic_libraries",
     //           "flex_helix_library",     // these are commented out for now... they will eventually be re-implemented CJ 09/20
                "ideal_helices"
    };

};
#endif // __BUILD_SQLITE_LIBRARIES_H__