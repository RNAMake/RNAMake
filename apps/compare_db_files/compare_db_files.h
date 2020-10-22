#ifndef RNAMAKE_COMPARE_DB_FILES
#define RNAMAKE_COMPARE_DB_FILES

#include <filesystem>
#include <vector>

#include <base/types.h>
#include <base/application.hpp>
#include <base/log.h>
#include <resources/motif_sqlite_library.h>

#include <CLI/CLI.hpp>

class CompareDBFiles : public base::Application {
public:
    CompareDBFiles() :
                app_("compare_db_files") {

    }

    void
    setup_options() override {
        app_.add_option("MOTIF-DIR",parameters_.dir_1)
            ->check(CLI::ExistingDirectory)
            ->required();

        app_.add_option("--log_level",parameters_.log_level,"level for global logging")
                ->check(CLI::IsMember(std::set<String>{"debug","error","fatal","info","verbose","warn"}))
                ->default_val("info")
                ->group("Core Inputs");


    }

    void
    run() override ;

private:
    std::vector<std::filesystem::path> dir_1_files_;
private:
    struct Parameters {
        std::filesystem::path dir_1;
        String log_level;
    };

    Parameters parameters_;
private:
    std::vector<std::filesystem::path>
    _get_paths(std::filesystem::path const& path ) const {
        auto db_files = std::vector<std::filesystem::path>();

        for(auto& file : std::filesystem::recursive_directory_iterator( path,
                            std::filesystem::directory_options::skip_permission_denied))  {
            const auto& fp = file.path();
            if(fp.extension() == ".db") {
                db_files.emplace_back(fp.string().substr(path.string().size()));
            }
        }

        return db_files;
    }

private:
    static
    std::vector<std::filesystem::path>
    _get_files()  {
        auto files = std::vector<std::filesystem::path>();
        for (const auto &pr : resources::MotifSqliteLibrary::get_libnames()) {
            files.emplace_back(pr.second);
        }
        return files;
    }

public:
    CLI::App app_ ;

};

#endif // RNAMAKE_COMPARE_DB_FILES