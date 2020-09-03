// Made by Chris Jurich 08/25/2020

#ifndef __BUILD_SQLITE_LIBRARIES_H__
#define __BUILD_SQLITE_LIBRARIES_H__

#include <set>
#include <unordered_set>
#include <map>
#include <dirent.h>
#include <sstream>
#include <functional>

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

#include <CLI/CLI.hpp>

struct BuildOpt {
    std::function<void()> handle;
    bool selected = false;
};

class BuildSqliteLibraries : public base::Application {

public:
    BuildSqliteLibraries() :
    lib_names_(resources::MotifSqliteLibrary::get_libnames()),
    app_("BuildSqliteLibraries"),
    build_opts_({
            // working methods
            {"helix_ensembles",                BuildOpt{[this] { _build_helix_ensembles(); },false}},
            {"ideal_helices",                  BuildOpt{[this] { _build_ideal_helices(); },false}},
            {"basic_libraries",                BuildOpt{[this] { _build_basic_libraries(); },false}},
            {"helix_ensembles",                BuildOpt{[this] { _build_helix_ensembles(); },false}},
            {"new_bp_steps",                   BuildOpt{[this] { _build_new_bp_steps(); },false}},
            {"motif_state_libraries",          BuildOpt{[this] { _build_motif_state_libraries(); },false}},
            {"unique_twoway_library",          BuildOpt{[this] { _build_unique_twoway_library(); },false}},
            {"ss_and_seq_libraries",           BuildOpt{[this] { _build_ss_and_seq_libraries(); },false}},
            {"trimmed_ideal_helix_library",    BuildOpt{[this] { _build_trimmed_ideal_helix_library();},false}}
           // methods that don't work or don't have working python versions
                //{"existing_motif_library",         BuildOpt{[this] { build_existing_motif_library(); },false}},
                //{"flex_helix_library",             BuildOpt{[this] { build_flex_helix_library(); },false}},
                //{"le_helix_lib",                   BuildOpt{[this] { build_le_helix_lib(); },false}},
                //{"motif_ensemble_state_libraries", BuildOpt{[this] { _build_motif_ensemble_state_libraries(); },false}},
    })
    {

    }

public:
    void
    run() override ;

public:
    void
    setup_options() override ;

public:
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

    void
    _directory_setup() {
        for( const auto& dir : {"/motif_librariesV2/","motif_state_librariesV2/","/motif_ensemble_librariesV2/"}) {
            const auto& full_dir = base::resources_path() + dir ;
            if(!base::is_dir(full_dir)) {
                const auto status = mkdir(full_dir.c_str(),0);

                if(status) {
                    LOGF<<"Unable to create directory: "<<full_dir;
                    LOGF<<"Exiting";
                }
            }
        }
    }

public:
    base::LogLevel
    log_level() const {
        return base::log_level_from_str(parameters_.log_level);
    }

private:
    struct Parameters {
        bool  all = false;
        String input_dir;
        String log_level;
    };

public:
    CLI::App app_;

private:
    StringStringMap
    _get_valid_dirs(String const&);

private:
    StringStringMap lib_names_;
    std::map<String,motif::MotifOP> motif_map_;
    std::map<String,BuildOpt> build_opts_;
    Parameters parameters_;

private:
    const Strings motif_keys_{"data","name","end_name","end_id","id"};
    const Strings ensemble_keys_{"data","name","id"};
};
#endif // __BUILD_SQLITE_LIBRARIES_H__