
// Made by Chris Jurich 08/25/2020



#ifndef __BUILD_SQLITE_LIBRARIES_H__
#define __BUILD_SQLITE_LIBRARIES_H__

#include <set>
#include <unordered_set>
#include <map>
#include <dirent.h>


#include <base/file_io.h>
#include <base/env_manager.h>
#include <base/log.h>
#include <base/types.h>
#include <base/settings.h>
#include <base/application.hpp>

#include <resources/motif_sqlite_library.h>
#include <resources/sqlite_library.h>
#include <motif/motif_factory.h>


#include <argparse/argparse.hpp>

class BuildSqliteLibraries : public base::Application {

public:
    BuildSqliteLibraries() :
    lib_names_(resources::MotifSqliteLibrary::get_libnames())
    {}

public:
    void
    run() override ;

    void
    setup_options() override ;

    void
    build_helix_ensembles();

    void
    build_ideal_helices();

    void
    build_trimmed_ideal_helix_library();

    void
    build_new_bp_steps();

private:
    StringStringMap
    _get_valid_dirs(String const&);

private:
    argparse::ArgumentParser options; 
    StringStringMap lib_names_; 
    std::map<String,motif::MotifOP> motif_map_;

};
#endif // __BUILD_SQLITE_LIBRARIES_H__
