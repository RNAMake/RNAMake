#include <build_sqlite_libraries/build_sqlite_libraries.h>

#include <utility>
#include <motif/motif_ensemble.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper structs brought over from python
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MotifEnsembles = std::vector<motif::MotifEnsemble>;

struct MotifandEnd{
    motif::MotifOP motif;
    int index;
    MotifandEnd(motif::MotifOP  motif, int index) : motif(std::move(motif)), index(index) {

    }
};

using MotifandEnds = std::vector<MotifandEnd>;

class SSandSeqCluster {
    public:
        String id;
        MotifandEnds motif_and_ends;
    public:
        SSandSeqCluster(String  id) :
        id(std::move(id)), motif_and_ends(MotifandEnds{}) {}

    public:
        bool
        motif_matches_end(motif::MotifOP const & motif, int end_index) {
            if(id == motif->end_ids()[end_index]) {
                motif_and_ends.emplace_back(motif,end_index); return true; } else {
                return false;
            }
        }

    public:
        bool
        motif_matches(motif::MotifOP const& motif) {
            const auto end_id_len = motif->end_ids().size();
            for(auto ii = 0; ii<end_id_len; ++ii) {
                if(motif_matches_end(motif,ii)) {
                    return true;
                }
            }
            return false;
        }

};

using SSandSeqClusters = std::vector<SSandSeqCluster>;

struct MotifCluster {
    String end_id;
    MotifandEnds motif_and_ends;

    MotifCluster(motif::MotifOP const& motif, structure::BasepairStateOP state) : state(std::move(state))
    {
        motifs.push_back(motif);
    }

    structure::BasepairStateOP state; 
    motif::MotifOPs motifs;
};

using MotifClusters = std::vector<MotifCluster>;

MotifClusters
cluster_motifs(resources::MotifSqliteLibrary::iterator && start,
                    resources::MotifSqliteLibrary::iterator const & end,
                    float max_distance=1.5f) {
    
    if(start == end) {
        return MotifClusters{};
    }
    
    auto clusters = MotifClusters{};
    clusters.emplace_back(*start, (*start)->ends()[0]->state());

    for(; start != end; ++start) {

        auto found(0);

        for(auto& c : clusters) {
            const auto dist = c.state->diff((*start)->ends()[1]->state());
            if(dist < max_distance) {
                ++found;
                c.motifs.push_back((*start));
                break;
            }
        }
        if(!found) {
            clusters.emplace_back(MotifCluster{(*start),(*start)->ends()[1]->state()});
        }
    }

    return clusters;
}

MotifClusters
cluster_motifs(motif::MotifOPs const & motifs, float max_distance=1.5f) {

    if(motifs.empty()) {
        LOGW << "Empty list of motifs encountered";
        return MotifClusters{};
    }

    auto clusters = MotifClusters{};
    clusters.emplace_back(motifs[0], motifs[0]->ends()[0]->state());

    const auto num_motifs = motifs.size();
    for(auto ii = 0; ii<num_motifs; ++ii) {
        auto found(0);
        const auto& curr_motif = motifs[ii];
        for(auto& c : clusters) {
            const auto dist = c.state->diff((curr_motif)->ends()[1]->state());
            if(dist < max_distance) {
                ++found;
                c.motifs.push_back((curr_motif));
                break;
            }
        }

        if(!found) {
            clusters.emplace_back(MotifCluster{(curr_motif),(curr_motif)->ends()[1]->state()});
        }
    }

    return clusters;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// method definitions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
StringStringMap
BuildSqliteLibraries::_get_valid_dirs(String const& base_dir) {
    
    auto paths = StringStringMap{};

#if __has_include(<filesystem>)

    if(!std::filesystem::exists(base_dir)) {
        LOGW<<"WARNING: The directory "<<base_dir<<" does not exist.";
        return paths;
    }

    for(auto& dir : std::filesystem::directory_iterator(base_dir)) {
        const String dir_name = std::filesystem::path(dir).string(); ;
        if(std::filesystem::exists(dir_name + "/ref_frames.dat")) {
            const auto last_slash = dir_name.find_last_of('/') + 1;
            const auto motif_name = dir_name.substr(last_slash);
            paths[motif_name] = dir_name;
        }
    }

#else
    DIR *dir = opendir(base_dir.c_str());

    struct dirent *entry = readdir(dir);

    while (entry != nullptr)
    {
        if (entry->d_type == DT_DIR) {
            const String dir_name = base_dir +"/"+ entry->d_name;
            if(base::file_exists(String{dir_name} + String{"/ref_frames.dat"})) {
                const auto last_slash = dir_name.find_last_of('/') + 1;
                const auto motif_name = dir_name.substr(last_slash);
                paths[motif_name] = dir_name;
            }
            
        }
        entry = readdir(dir);
    }

    closedir(dir);
#endif
    return paths;
}

void
BuildSqliteLibraries::_build_ideal_helices() {

    auto paths = _get_valid_dirs(parameters_.input_dir + "/helices/");

    auto mf = motif::MotifFactory{};
    auto data = std::vector<Strings>{};
    auto reversed_data = std::vector<Strings>{};
    auto motif_ct(0); 

    for(auto& [folder, path] : paths) {
        auto tokens = base::split_str_by_delimiter(folder,".");

        auto motif = mf.motif_from_file(path);
        
        const auto last_dot = folder.find_last_of('.') + 1;
        const auto num = folder.substr(last_dot);
        motif->name(
            "HELIX.IDEAL" + (num == "0" ? "" : "." + num)
                );

        auto aligned_motif = mf.align_motif_to_common_frame(motif,0); 
        
        auto entry = Strings{
                aligned_motif->to_str(), aligned_motif->name(), aligned_motif->end_name(0), aligned_motif->end_ids()[0],std::to_string(motif_ct)
        };

        resources::sqlite3_escape(entry);

        data.push_back(entry);

        auto reversed_motif = mf.can_align_motif_to_end(aligned_motif,1);
        auto reversed_aligned_motif = mf.align_motif_to_common_frame(reversed_motif,0);
        mf._setup_secondary_structure(reversed_aligned_motif);
        auto reversed_entry = Strings{
            reversed_aligned_motif->to_str(),reversed_aligned_motif->name(),reversed_aligned_motif->end_name(0),reversed_aligned_motif->end_ids()[0],std::to_string(motif_ct)
        };

        resources::sqlite3_escape(reversed_entry);

        reversed_data.push_back(reversed_entry);

        ++motif_ct;
    }
    
    const auto path = parameters_.output_dir +"/motif_libraries_new/ideal_helices.db";

    resources::build_sqlite_library(path,data,motif_keys_,"id");

    const auto reversed_path = parameters_.output_dir +"/motif_libraries_new/ideal_helices_reversed.db";
    resources::build_sqlite_library(reversed_path,reversed_data,motif_keys_,"id");

    auto sql_lib = resources::MotifSqliteLibrary{"ideal_helices"};
    auto motif = sql_lib.get("HELIX.IDEAL.3");

    auto outfile = std::ofstream(parameters_.output_dir + "/motifs/base.motif");
    outfile<<motif->to_str();
    outfile.close();
}

void
BuildSqliteLibraries::_build_unique_twoway_library() {
    auto valid_dirs = _get_valid_dirs(parameters_.input_dir + "/two_ways/");
    auto motifs = motif::MotifOPs{};
    auto mf = motif::MotifFactory{};
    for(auto&& [motif_name,path] : valid_dirs) {
        if(motif_name.find("TWOWAY") != std::string::npos) {
            try {
                motifs.push_back(
                        mf.motif_from_file(path)
                );
            } catch(std::runtime_error& error ) {
                LOGW<<error.what()<<" in file "<<path;
            }
        }
    }
    auto clusters = cluster_motifs(motifs,9.0f);

    auto data = std::vector<Strings>{};
    auto mes_data = std::vector<Strings>{};
    auto outfile = std::ofstream("sim_list_new");
    auto count(0);

    for(auto& c : clusters) {
        auto lowest = c.motifs[0];
        for(const auto& m : c.motifs) {
            if (lowest->score() > m->score()) {
                lowest = m;
            }
        }

        if(lowest->name() == "TWOWAY.1GID.6" ||
            lowest->name() == "TWOWAY.1GID.2" ||
            lowest->name() == "TWOWAY.2GDI.4" ||
            lowest->name() == "TWOWAY.2VQE.18" ) {
            continue;
        }

        ++count;

        outfile<<lowest->name()<<","<<lowest->ends()[0]->name()<<"|";

        for(const auto& m : c.motifs) {
            outfile<<lowest->name()<<","<<lowest->ends()[0]->name()<<"|";
        }
        outfile<<"\n";

        auto entry = Strings{
                lowest->to_str(), lowest->name(), lowest->ends()[0]->name(), lowest->end_ids()[0], std::to_string(count)
        };

        resources::sqlite3_escape(entry);

        data.push_back(entry);
    }

    outfile.close();
    const auto path = parameters_.output_dir + "/motif_libraries_new/unique_twoway.db";
    resources::build_sqlite_library(path, data, motif_keys_,"id");
}

void
BuildSqliteLibraries::_build_new_bp_steps() {
    auto motif_lib = resources::MotifSqliteLibrary{"bp_steps"};
    motif_lib.load_all();

    auto motifs = motif::MotifOPs{};
    for(auto& motif : motif_lib)  {
        auto tokens = base::split_str_by_delimiter(motif->name(),".");
        if(tokens.size() == 1) {
            motifs.push_back(motif); // this branch is never used... pretty sure this logic fro plython doesn't hold up CJ 09/20
        }
    }
    
    auto ii(0);
    auto motif_ct(0);
    auto unique = std::unordered_set<String>();
    auto mf = motif::MotifFactory{}; 

    auto motif_data = std::vector<Strings>{};

    for(auto& motif : motifs) {
        auto tokens = base::split_str_by_delimiter(motif->name(),"=");
        
        const auto& ends = motif->end_ids();

        if(unique.find(ends[0]) != unique.end() || 
            unique.find(ends[1]) != unique.end() ) {
            continue;
        }
        
        const auto old_name = motif->name();
        unique.insert(ends[0]); 
        
        if (unique.find(ends[1]) == unique.end()) {
            unique.insert(ends[1]);
        }
    
        motif->name("BP." + std::to_string(ii));
        
        auto motif_aligned = mf.can_align_motif_to_end(motif,1);
        motif_aligned = mf.can_align_motif_to_end(motif_aligned,1);

        auto entry = Strings{
                motif->to_str(), motif->name(), motif->end_name(0), motif->end_ids()[0], std::to_string(motif_ct++)
        };

        resources::sqlite3_escape(entry);

        motif_data.push_back(entry);

        auto aligned_entry = Strings{
                motif_aligned->to_str(), motif_aligned->name(), motif_aligned->end_name(0), motif_aligned->end_ids()[0],std::to_string(motif_ct++)
        };

        resources::sqlite3_escape(aligned_entry) ;

        motif_data.push_back(aligned_entry);
        
        ++ii;
     }

    const auto path = parameters_.output_dir +"/motif_libraries_new/new_bp_steps.db";

    resources::build_sqlite_library(path, motif_data, motif_keys_, "id");

}

void
BuildSqliteLibraries::_build_flex_helix_library() {
    // DETAIL this one isn't implemented bc the files don't exist --CJ
    throw base::RNAMakeImplementationExcepetion("the method: \"BuildSqliteLibraries::_build_flex_helix_library()\" is not fully implemented. Do not call it ");
}

void
BuildSqliteLibraries::_build_le_helix_lib() {

    auto paths = _get_valid_dirs(parameters_.input_dir + "/helices/");

    auto mf = motif::MotifFactory{};
    auto data = std::vector<Strings>{};
    auto reversed_data = std::vector<Strings>{};
    auto motif_ct(0);

    for (auto&[folder, path] : paths) {
        auto tokens = base::split_str_by_delimiter(folder, ".");

        if (tokens[1] != "LE") {
            continue;
        }

        auto motif = mf.motif_from_file(path);

        auto successes = std::vector<std::pair<motif::MotifOP, int>>{};

        for (auto ii = 0; ii < motif->ends().size(); ++ii) {
            auto aligned_motif = mf.align_motif_to_common_frame(motif, ii);

            if (aligned_motif != nullptr) {
                successes.emplace_back(aligned_motif, ii);
            }
        }

        for (auto & succ : successes) {
            const auto &motif = succ.first;
            const auto &index = succ.second;
            const auto aligned_motif = mf.align_motif_to_common_frame(motif, index);
            auto entry = Strings{aligned_motif->to_str(),
                                 aligned_motif->name(),
                                 aligned_motif->ends()[0]->name(),
                                 aligned_motif->end_ids()[0],
                                 std::to_string(motif_ct++)};
            resources::sqlite3_escape(entry);

            data.emplace_back(entry);
        }


    }

    const auto path = parameters_.output_dir + "/motif_libraries_new/le_helices.db";
    resources::build_sqlite_library(path, data, motif_keys_, "id" );
}

    void
BuildSqliteLibraries::_build_motif_ensemble_state_libraries() {

    throw base::RNAMakeImplementationExcepetion("the method: \"BuildSqliteLibraries::_build_motif_ensemble_state_library()\" is not fully implemented. Do not call it ");

    // not really sure how to implement this one... couldn't find a counterpart for MotifEnsembleSqliteLIbrary -CJ -09/20
    //       for libname in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
//    print libname
//    me_lib = sqlite_library.MotifEnsembleSqliteLibrary(libname)
//    me_lib.load_all()
//
//    keys = ['data', 'name', 'id']
//    data = []
//
//    for i, me in enumerate(me_lib.all()):
//        mse = me.get_state()
//        #print len(me.members)
//        data.append([mse.to_str(), mse.id, i])
//
//    path = settings.RESOURCES_PATH +"/motif_state_ensemble_libraries/"+libname+".db"
//    sqlite_library.build_sqlite_library(path, data, keys, 'id')

}

void
BuildSqliteLibraries::_build_existing_motif_library() {
    throw base::RNAMakeImplementationExcepetion("the method: \"BuildSqliteLibraries::_build_existing_motif_library()\" is not fully implemented. Do not call it ");
    if(!base::file_exists(base::resources_path() + "/existing.motifs")) {
        LOGW<<"The file \"existing.motifs\"... unable to build_existing_motif_library";
    }
    const auto lines = base::get_lines_from_file("existing.motifs");
    auto motifs = motif::MotifOPs{};

    for(const auto& line : lines) {
    }
}

void
BuildSqliteLibraries::_build_basic_libraries() {
    const auto motif_types = util::MotifTypes {
        util::MotifType::TWOWAY,
        util::MotifType::NWAY,
        util::MotifType::HAIRPIN,
        util::MotifType::TCONTACT
    };

    const auto bad_keys = std::unordered_set<String>{
            "TWOWAY.2GDI.4-X20-X45",
            "TWOWAY.1S72.46-02097-02647",
            "TWOWAY.2GDI.6-Y20-Y45"
    };

    auto file_lines = Strings{};
    if(std::filesystem::exists(parameters_.input_dir + "/motifs_extra_bps.csv")) {
        file_lines = base::get_lines_from_file("motifs_extra_bps.csv");
        file_lines.erase(file_lines.begin(),file_lines.begin()+1);
    } else {
        LOGW<<"WARNING: The file "<<parameters_.input_dir <<
        "/motifs_extra_bps.csv does not exist. This will limit what can be done in BuildSqliteLibraries::_build_basic_libraries().";
    }

    auto has_extra_bps = std::map<String,int>{};

    for(const auto& line : file_lines ) {
        auto tokens = base::split_str_by_delimiter(line,",");
        if(tokens.size() > 1 && std::stof(tokens[2]) > 0 ) {
            has_extra_bps[tokens[0]] = 1;
        }
    }

    auto folders = StringStringMap{};
    //TODO change this.. a bit hacky but it works for now CJ 09/20
    auto path = std::filesystem::recursive_directory_iterator(parameters_.input_dir,
                                                               std::filesystem::directory_options::skip_permission_denied);

    for(auto& dir : std::filesystem::directory_iterator(parameters_.input_dir,std::filesystem::directory_options::skip_permission_denied)) {
        if(!std::filesystem::is_directory(dir)) {
            continue;
        }
        auto entries = _get_valid_dirs(dir.path().string());
        folders.insert(entries.begin(),entries.end());
    }

    auto mf = motif::MotifFactory{};
    for(const auto& type : motif_types) {
        const auto type_str = util::type_to_str(type);

        auto successes = std::vector<std::pair<motif::MotifOP,int>>{};

        for(const auto& [dir, file] : folders) {
            if(type_str == "TCONTACT" && dir.find("TC.") == std::string::npos) {
                continue;
             } else if (dir.find(type_str) == std::string::npos && type_str != "TCONTACT") {
                continue;
            }

            if(type == util::MotifType::TWOWAY &&
                    has_extra_bps.find(dir) != has_extra_bps.end()) {
                    continue;
            }

            auto motif = motif::MotifOP{};
            try {
                motif = mf.motif_from_file(file);
            } catch (std::runtime_error& error) {
                LOGW<<"Error: "<<error.what()<<". Couldn't build motif for file "<<file<<". Skipping.";
                continue;
            }

            const auto& num_ends = motif->ends().size();

            if( (type != util::MotifType::HAIRPIN && num_ends == 1 ) || num_ends == 0 ) {
                continue;
            }


            for(auto ii = 0; ii<num_ends; ++ii ) {
                auto motif_added = mf.can_align_motif_to_end(motif,ii);

                if(motif_added == nullptr)  {
                    continue;
                }

                const auto key = String{
                    motif_added->name() + "-" + motif_added->ends()[0]->name()
                };

                if(bad_keys.find(key) != bad_keys.end()) {
                    continue;
                }

                successes.emplace_back(motif_added,ii);
            }
        }

        auto data = std::vector<Strings>{};


        auto motif_ct(0);

        for(const auto& entry : successes) {

            auto motif = entry.first;
            auto index = entry.second;

            auto m_added = mf.align_motif_to_common_frame(motif,index);

            if(type == util::MotifType::TWOWAY) {
                auto remove = structure::BasepairOPs{};
                for(const auto& bp : m_added->basepairs()) {

                    if(std::find(m_added->ends().begin(),m_added->ends().end(),bp) == m_added->ends().end()
                            && bp->bp_type() == "cW-W") {

                        remove.push_back(bp);
                    }
                }

                m_added->remove_bad_bps(remove);

                try {
                    mf._setup_secondary_structure(m_added);
                } catch(std::runtime_error& error) {
                    LOGW<<"Error: "<<error.what();
                    continue;
                }

            }

            auto data_entry = Strings{
                    m_added->to_str(), m_added->name(), m_added->ends()[0]->name(), m_added->end_ids()[0],std::to_string(motif_ct++)
            };
            resources::sqlite3_escape(data_entry);

            data.push_back(data_entry);

        }

        if(data.empty())  {
            LOGW<<"No data generated for motif type: "<<util::type_to_str(type);
            continue;
        }

        auto motif_name = util::type_to_str(type);
        std::transform(motif_name.begin(), motif_name.end(), motif_name.begin(),tolower);
        const auto path = parameters_.output_dir + "/motif_libraries_new/" + motif_name + ".db";
        resources::build_sqlite_library(path,data,motif_keys_,"id");
    }

}

void
BuildSqliteLibraries::_build_motif_state_libraries() {
    const auto path = std::filesystem::path(parameters_.output_dir + "/motif_state_libraries_new/");

    const auto total = lib_names_.size();
    auto index(1);
    for(const auto& pair : lib_names_) {
        LOGI<<"Beginning motif db file "<<index<<" of "<<total<<": "<<pair.first;
        auto data = std::vector<Strings>{};

        auto motif_lib = resources::MotifSqliteLibrary{pair.first};
        motif_lib.load_all();
        auto ii(0);
        for(const auto& motif : motif_lib) {
            const auto& motif_state = motif->get_state();

            auto data_entry = Strings{
                motif_state->to_str(),motif_state->name(),motif_state->end_names()[0],motif_state->end_ids()[0],std::to_string(ii++)
            };
            resources::sqlite3_escape(data_entry);
            data.push_back(data_entry);
        }
        LOGI<<"\tFound "<<data.size()<<" entries. Beginning to write to db file...";

        resources::build_sqlite_library(path.string() + pair.first + ".db", data, motif_keys_, "id");

        LOGI<<"\tWrote "<<data.size()<<" entries to file "<<path.string() + pair.first + ".db";
        LOGI<<"Finishing motif db file "<<index++<<" of "<<total<<": "<<pair.first;
    }
}

void
BuildSqliteLibraries::_build_ss_and_seq_libraries() {
    const auto libbames = Strings{"twoway","tcontact","hairpin","nway"};

    for(const auto& libname : libbames) {
        auto motif_lib = resources::MotifSqliteLibrary{libname};
        motif_lib.load_all();

        auto clusters = SSandSeqClusters{};
        auto motif_ensembles = MotifEnsembles{};
        auto mes_names = Strings{};
        auto motifs = motif::MotifOPs{};

        for(const auto& motif : motif_lib) {
            auto matched(false);

            for(auto& c : clusters) {
                if(c.motif_matches_end(motif,0)) {
                    matched = true;
                }
            }

            if(!matched) {
                clusters.emplace_back(motif->end_ids()[0]);
                clusters.rbegin()->motif_matches_end(motif,0);
            }
        }

        auto data = std::vector<Strings>{};

        const auto num_clusters = clusters.size();

        for(auto ii = 0; ii<num_clusters; ++ii) {
            const auto &cluster = clusters[ii];
            auto all_motifs = motif::MotifOPs{};

            for (const auto &motif_end : cluster.motif_and_ends) {
                all_motifs.push_back(motif_end.motif);
            }

            if (libname != "twoway") {
                auto energies = std::vector<float>(all_motifs.size(), 1);
                auto motif_ensemble = motif::MotifEnsemble{cluster.id,all_motifs, energies};

                motif_ensembles.push_back(motif_ensemble);
                motifs.push_back(motif_ensemble.members()[0]->motif);
                motifs.rbegin()->get()->name(motif_ensemble.id());
                mes_names.push_back(motif_ensemble.id());

                auto data_entry = Strings{
                        motifs[0]->to_str(), motif_ensemble.id(), std::to_string(ii)
                };

                resources::sqlite3_escape(data_entry);
                data.push_back(data_entry);
                continue;
            }

            auto m_clusters = cluster_motifs(all_motifs);
            auto clustered_motifs = motif::MotifOPs{};
            auto energies = std::vector<float>{};

            const auto num_sub_clusters = m_clusters.size();
            for (auto jj = 0; jj < num_sub_clusters; ++jj) {
                const auto &c_motifs = m_clusters[jj];
                clustered_motifs.push_back(c_motifs.motifs[0]);
                const auto pop = float(c_motifs.motifs.size()) / float(all_motifs.size());
                energies.push_back(-base::constants::kBT * std::log(pop));
            }

            auto motif_ensemble = motif::MotifEnsemble{cluster.id, clustered_motifs, energies};

            auto data_entry = Strings{
                    motif_ensemble.to_str(), motif_ensemble.id(), std::to_string(ii)
            };

            resources::sqlite3_escape(data_entry);
            data.push_back(data_entry);
        }
        const auto path = parameters_.output_dir + "/motif_ensemble_libraries_new/" + libname + ".db";
        resources::build_sqlite_library(path,data,ensemble_keys_,"id");
    }
}

void
BuildSqliteLibraries::_build_trimmed_ideal_helix_library() {
    auto helix_folders = _get_valid_dirs(parameters_.input_dir + "/motifs/helices/" ) ;
    const auto lib_name =  lib_names_.find("ideal_helices");
    if( lib_name == lib_names_.end()) {
        LOGW <<"\"ideal_helices\" not found in library names. No ideal helices built";
        return;
    }
    auto sql_lib = resources::MotifSqliteLibrary{"ideal_helices"};
    auto mf = motif::MotifFactory{};
    auto motif_ct(0);
    auto data = std::vector<Strings>{};

    for(auto&& folder_dir : helix_folders ) {
        if(folder_dir.first.find("IDEAL") == std::string::npos &&
                folder_dir.first.find("ideal") == std::string::npos) {
            continue;
        }

        auto m = mf.motif_from_file(folder_dir.second);
        auto entry = Strings{
                m->to_str(), m->name(), m->end_name(0), m->end_ids()[0], std::to_string(motif_ct++)
        };

        resources::sqlite3_escape(entry);

        data.push_back(entry);
    }

    const auto path = parameters_.output_dir  + "/motif_state_libraries_new/" + lib_name->first + "_min.db";
    resources::build_sqlite_library(path,data,motif_keys_,"id");
}

void
BuildSqliteLibraries::_build_helix_ensembles() {
    const auto hel_dir = std::filesystem::path(parameters_.input_dir + "/helices/");
    if( ! std::filesystem::exists(hel_dir) ) {
        LOGW<<"WARNING: directory "<<hel_dir.string()<<" does not exist. Cannot build helix ensembles. Exiting...";
    }

    const auto helix_paths = _get_valid_dirs(hel_dir.string());

    auto clusters = SSandSeqClusters{};
    constexpr auto kB = 1.3806488e-1;
    constexpr auto kBT = kB * 298.15;

    auto mf = motif::MotifFactory{};

    LOGI<<"Beginning cluster construction. Found "<<helix_paths.size()<<" files to use...";

    for(const auto& [dir,path] : helix_paths ) {
        const auto tokens = base::split_str_by_delimiter(dir,".");

        if(tokens.size() > 1 && (tokens[1] == "IDEAL" || tokens[1] == "LE")) {
            continue;
        }
        auto motif = mf.motif_from_file(path);
        const auto num_bps = motif->basepairs().size();
        for(auto ii = 0; ii<num_bps - 1; ++ii) {
            const auto& basepairs = motif->basepairs();
            if(basepairs[ii]->bp_type() != "cW-W" || basepairs[ii+1]->bp_type() != "cW-W") {
                continue;
            }
            auto bps = structure::BasepairOPs{basepairs[ii],basepairs[ii+1]};

            auto res = structure::ResidueOPs{};

            for(const auto& bp : bps)  {
                for(const auto& r : bp->residues()) {
                    if(std::find(res.begin(),res.end(),r) == res.end()) {
                        res.push_back(r);
                    }
                }
            }

            if(res.size() != 4) {
                continue;
            }

            auto m_bps = mf.motif_from_bps(bps);

            if(m_bps->end_ids().size() != 2 || m_bps->end_ids()[0].size() != 11) {
                continue;
            }

            auto matched(false);

            for( auto& c : clusters) {
                if(c.motif_matches(m_bps)) {
                   matched = true;
                }
            }

            if( not matched) {
                clusters.emplace_back(m_bps->end_ids()[0]);
                clusters.back().motif_matches(m_bps);
                if (m_bps->end_ids()[0] == m_bps->end_ids()[1]) {
                    clusters.emplace_back(m_bps->end_ids()[1]);
                    clusters.back().motif_matches(m_bps);
                }
            }

        }
    }

    LOGI<<"Finished cluster construction for "<<helix_paths.size()<<" files.";

    auto mes_data = std::vector<Strings>{};
    auto all_mes_data = std::vector<String>{};
    auto all_count(0);
    auto motif_data = std::vector<Strings>{};
    auto unique_data = std::vector<Strings>{};
    auto unique = std::unordered_set<String>{};
    auto bp_count(1);

    auto outfile = std::fstream("bp_step_info.csv");
    outfile<<"name,bp_type,pop,d\n";

    auto count(1);

    const auto num_clusters = clusters.size();

    LOGI<<"Beginning cluster calculations. There are "<<num_clusters<<" clusters ...";

    for(auto ii=0; ii<num_clusters; ++ii) {

        const auto &cluster = clusters[ii];
        auto &motif = cluster.motif_and_ends[0].motif;

        if (unique.find(motif->end_ids()[0]) != unique.end()
            || unique.find(motif->end_ids()[1]) != unique.end()) {
            continue;
        }

        unique.insert(motif->end_ids()[0]);
        unique.insert(motif->end_ids()[1]);

        ++all_count;

        auto tokens = base::split_str_by_delimiter(cluster.id, "_");

        const auto num_motif_and_ends = cluster.motif_and_ends.size();
        auto aligned_motifs = motif::MotifOPs{};

        for (auto jj = 0; jj < num_motif_and_ends; ++jj) {

            const auto &motif_and_end = cluster.motif_and_ends[jj];
            const auto &m_and_e_motif = motif_and_end.motif;
            const auto &m_and_e_index = motif_and_end.index;

            auto motif_aligned = motif::MotifOP{};
            try {
                motif_aligned = mf.can_align_motif_to_end(m_and_e_motif, m_and_e_index);
            } catch (std::runtime_error &error) {
                LOGW << "Error: " << error.what();
                continue;
            }

            if (motif_aligned == nullptr) {
                continue;
            }

            motif_aligned = mf.align_motif_to_common_frame(motif_aligned, m_and_e_index);
            auto fail(false);

            for (const auto &bp : motif_aligned->basepairs()) {

                const auto vec1 = bp->res1()->get_atom("C2'")->coords() - bp->res1()->get_atom("O4'")->coords();
                const auto vec2 = bp->res2()->get_atom("C2'")->coords() - bp->res2()->get_atom("O4'")->coords();

                if (vec1.dot(vec2) > 2.f) {
                    fail = true;
                    break;
                }
            }

            if (fail) {
                continue;
            }
            aligned_motifs.push_back(motif_aligned);
        }

        auto motif_clusters = cluster_motifs(aligned_motifs, 0.65f);

        auto clustered_motifs = motif::MotifOPs{};
        auto energies = std::vector<float>{};
        auto ss = std::stringstream{};
        ss<<tokens[0][0]<<tokens[2][1]<<"="<<tokens[0][1]<<tokens[2][0];
        String dir_name = ss.str();

        const auto num_motif_clusters = motif_clusters.size();
        for (auto jj = 0; jj < num_motif_clusters; ++jj) {
            auto cluster_motifs = motif_clusters[jj].motifs;
            auto curr_motif = cluster_motifs[0];
            curr_motif->mtype(util::MotifType::HELIX);
            curr_motif->name("BP." + std::to_string(ii) + "." + std::to_string(jj));

            auto motif_entry = Strings{
                    curr_motif->to_str(), curr_motif->name(), curr_motif->ends()[0]->name(), cluster.id,
                    std::to_string(count++)
            };

            resources::sqlite3_escape(motif_entry);

            motif_data.push_back(motif_entry);

            clustered_motifs.push_back(curr_motif);

            const auto pop = float(cluster_motifs.size()) / float(aligned_motifs.size());

            outfile << curr_motif->name() << "," << curr_motif->end_ids()[0] << "," << pop << ","
                    << curr_motif->ends()[1]->d().to_str();
            outfile << "," << curr_motif->ends()[1]->r() << "\n";

            energies.push_back(-kBT * std::log(pop));
        }


        //auto motif_ensemble = MotifEnsemble{};

        if(clustered_motifs.empty()) {
            continue;
        }
        auto motif_ensemble = motif::MotifEnsemble{cluster.id, clustered_motifs, energies};
        auto e_motif = motif_ensemble.members()[0]->motif;
        auto pop = std::exp(motif_ensemble.members()[0]->energy / -kBT);
        e_motif->name("BP." + std::to_string(ii));
        e_motif->to_pdb("BP." + std::to_string(ii) + ".pdb");

        outfile << e_motif->name() << "," << e_motif->end_ids()[0] << "," << std::to_string(pop) << ","
                << e_motif->ends()[1]->d().to_str();
        outfile << "," << e_motif->ends()[1]->r().to_str() << "\n";

        auto mes_entry = Strings{
                motif_ensemble.to_str(),motif_ensemble.id(),std::to_string(bp_count)
        };

        resources::sqlite3_escape(mes_entry);

        mes_data.push_back(mes_entry);

        auto unique_data_entry = Strings{
                e_motif->to_str(),e_motif->name(),e_motif->ends()[0]->name(), motif_ensemble.id(), std::to_string(bp_count++)
        };

        resources::sqlite3_escape(unique_data_entry);

        unique_data.push_back(unique_data_entry);

        clustered_motifs.clear();
        energies.clear();

        for (const auto &member : motif_ensemble.members()) {

            auto motif_aligned = mf.can_align_motif_to_end(member->motif, 1);
            motif_aligned = mf.can_align_motif_to_end(motif_aligned, 1);
            clustered_motifs.push_back(motif_aligned);
            energies.push_back(member->energy);
            pop = std::exp(member->energy / -kBT);

            if (e_motif->end_ids()[0] != e_motif->end_ids()[1]) {
                outfile << motif_aligned->name() << "," << motif_aligned->end_ids()[0] << "," << pop << ",";
                outfile << motif_aligned->ends()[1]->d().to_str() << "," << motif_aligned->ends()[1]->r().to_str()
                        << "\n";
            }

            auto motif_data_entry = Strings{
                    motif_aligned->to_str(),motif_aligned->name(),motif_aligned->ends()[0]->name(),motif_aligned->end_ids()[0],std::to_string(count++)
            };

            resources::sqlite3_escape(motif_data_entry);

            motif_data.push_back(motif_data_entry);
        }

        motif_ensemble = motif::MotifEnsemble{clustered_motifs[0]->end_ids()[0],clustered_motifs,energies};

        if(motif->end_ids()[0] != motif->end_ids()[1]) {
            pop = std::exp((*motif_ensemble.members().rbegin())->energy/-kBT);
            auto mes_data_entry = Strings{
                motif_ensemble.to_str(),motif_ensemble.id(),std::to_string(bp_count)
            };

            resources::sqlite3_escape(mes_data_entry);

            *motif = *motif_ensemble.members()[0]->motif; // not sure this is the right way to do it?? CJ 09/20

            motif->name(String{"BP."} + std::to_string(clusters.size() - 1));
            outfile<<motif->name()<<","<<motif->end_ids()[0]<<",";
            outfile<<motif->ends()[1]->d().to_str()<<","<<motif->ends()[1]->r().to_str()<<std::endl;
        }

        auto unique_entry = Strings{
            motif->to_str(),motif->name(),motif->ends()[0]->name(),motif_ensemble.id(),std::to_string(bp_count++)
        };

        resources::sqlite3_escape(unique_entry);

        unique_data.push_back(unique_entry);

    }

    LOGI<<"Finished cluster calculations for "<<num_clusters;

    String path = parameters_.output_dir + "/motif_ensemble_libraries_new/bp_steps.db";
    resources::build_sqlite_library(path, mes_data, ensemble_keys_, "id");

    path = parameters_.output_dir + "/motif_libraries_new/bp_steps.db";
    resources::build_sqlite_library(path, motif_data, motif_keys_, "id");

    path = parameters_.output_dir + "/motif_libraries_new/new_bp_steps.db";
    resources::build_sqlite_library(path, unique_data, motif_keys_, "id");
}

void
BuildSqliteLibraries::setup_options()   {

    app_.add_flag("--all",parameters_.all,"flag to build all available libraries")
            ->group("Core Inputs");

    app_.add_option("--in_dir",parameters_.input_dir,"directory that contains input data folders")
            ->default_val(base::resources_path())
            ->check(CLI::ExistingDirectory)
            ->group("Core Inputs");

    app_.add_option("--out_dir",parameters_.output_dir,"directory where the files will be published to")
            ->required()
            ->group("Core Inputs");

    app_.add_option("--log_level",parameters_.log_level,"level for global logging")
            ->check(CLI::IsMember(std::set<String>{"debug","error","fatal","info","verbose","warn"}))
            ->default_val("info")
            ->group("Core Inputs");

    for(auto& pair : build_opts_ )  {
        app_.add_flag("--" + pair.first,pair.second.selected)
                ->group("Build Options");
    }

}

void
BuildSqliteLibraries::run()   {

    base::init_logging(log_level());

    _generate_method_dependency_tree();

    for(auto& pairs : build_opts_) {
        if(parameters_.all || pairs.second.selected) {
            pairs.second.select();
            for(auto& dir_option : pairs.second.output_dirs) {
                required_dirs_.insert(dir_option);
            }
        }
    }

    _directory_setup();

    const auto total = std::count_if(build_opts_.cbegin(),build_opts_.cend(), [] (const auto& pair) {
        return pair.second.selected;
    });

    if(total == 0)  {
       LOGW<<"WARNING: No build options were selected, nothing will be done. Run ./build_sqlite_libraries -h to see options";
       LOGI<<"Exiting...";
       exit(1);
    }

    auto index(1);

    for(const auto& method : method_order_)  {

        auto& build = build_opts_.at(method);
        if(build.selected) {
            LOGI << "Beginning build " << index << " of " << total << ": " << method << " ...";
            build.handle();
            LOGI << "Finishing build " << index++ << " of " << total << ": " << method;
        }
    }

}

int main(int argc, const char** argv) {
    /* NOTES: existing.motifs does not exist... did not ressurrect this. Should ? Therefore, this part does not match
     * */
    std::set_terminate(base::save_backtrace);

    auto env_mgr = base::EnvManager(Strings{"RNAMAKE"});    
    env_mgr.set_envs();
    auto app = BuildSqliteLibraries{};

    app.setup_options();
    CLI11_PARSE(app.app_,argc,argv);
    app.run();
}
