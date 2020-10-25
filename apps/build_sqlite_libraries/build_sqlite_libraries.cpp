#include <build_sqlite_libraries/build_sqlite_libraries.h>

#include <utility>
#include <motif/motif_ensemble.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper structs brought over from python
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define FIX
using MotifEnsembles = std::vector<motif::MotifEnsemble>;

struct MotifandEnd{
    motif::MotifOP motif;
    int end_index;
    MotifandEnd(motif::MotifOP m, int ei) : motif(std::move(m)), end_index(ei) {

    }
};

using MotifandEnds = std::vector<MotifandEnd>;

class SSandSeqCluster {
    public:
        String end_id;
        MotifandEnds motif_and_ends;
    public:
        SSandSeqCluster(String  end_id) :
        end_id(std::move(end_id)), motif_and_ends(MotifandEnds{}) {}

    public:
        bool
        motif_matches_end(motif::MotifOP const & m, int ei) {
            if(end_id == m->end_ids()[ei]) {
                motif_and_ends.emplace_back(m, ei);
                return true;
            } else {
                return false;
            }
        }

    public:
        bool
        motif_matches(motif::MotifOP const& m) {
            for(auto i = 0; i<m->end_ids().size(); ++i) {
                auto r = motif_matches_end(m,i) ;
                if(r) {
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
    clusters.emplace_back(*start, (*start)->ends()[1]->state());
    ++start;
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
    auto clusters = MotifClusters{};
    if(motifs.empty()) {
        LOGW << "Empty list of motifs encountered";
        return clusters;
    }

    clusters.emplace_back(motifs[0], motifs[0]->ends()[1]->state());

    const auto num_motifs = motifs.size();
    for(auto i = 0; i<num_motifs; ++i) {
        if (i == 0) {
            continue;
        }
        auto found(0);
        const auto& m = motifs[i];
        for(auto& c : clusters) {
            const auto dist = c.state->diff((m)->ends()[1]->state());
            if(dist < max_distance) {
                ++found;
                c.motifs.push_back((m));
                break;
            }
        }

        if(found == 0)  {
            clusters.emplace_back(MotifCluster{(m),(m)->ends()[1]->state()});
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
BuildSqliteLibraries::_build_avg_helix_library() {
    auto lines = base::get_lines_from_file(base::resources_path() + "/setup/average_helices.dat");
    auto mf = motif::MotifFactory{};
    auto avg_helices = motif::MotifOPs{};

    for(auto i = 0; i < lines.size(); ++i) {
        try {
            auto m = std::make_shared<motif::Motif>(lines[i],
                                                    structure::ResidueTypeSetManager::getInstance().residue_type_set());
            avg_helices.push_back(std::move(m));
        } catch( ... ){
            LOGW<<"WARNING: Error occurred while trying to create Motif...";
            continue;
        }
    }

    auto succeses = std::vector<std::pair<motif::MotifOP, int>>{};
    for(auto& m : avg_helices) {
        auto m_added = mf.can_align_motif_to_end(m, 0);

        if (m_added == nullptr) {
            continue;
        } else {
            succeses.emplace_back(m_added,0);
        }
    }

    auto data = std::vector<Strings>{};

    for(auto i = 0 ; i < succeses.size(); ++i) {
        auto m = succeses[i].first;
        auto ei = succeses[i].second;

        auto m_added = mf.align_motif_to_common_frame(m,ei);

        auto data_entry = Strings{
            m_added->to_str(), m_added->name(), m_added->ends()[0]->name(), m_added->end_ids()[0], std::to_string(i)
        };

        resources::sqlite3_escape(data_entry);

        data.push_back(data_entry);
    }

    const auto path = parameters_.output_dir + "/motif_libraries_new/avg_helices.db";
    resources::build_sqlite_library(path, data, motif_keys_, "id");

}

void
BuildSqliteLibraries::_build_ideal_helices() {

    auto paths = _get_valid_dirs(parameters_.input_dir + "/helices/");
    auto mlib = motif::MotifOPs{};
    auto mf = motif::MotifFactory();
    for(int ii = 0 ; ii < 21; ++ii) {
        auto motif = mf.motif_from_file(parameters_.input_dir + "/helices/HELIX.LE." + std::to_string(ii));
        if(motif == nullptr)  {
            LOGW<<"Construction of LE "<<ii<<" failed";
        } else {
            mlib.push_back(motif);
        }

    }
    auto data = std::vector<Strings>{};
    auto rdata = std::vector<Strings>{};
    auto count(0);
    for(auto& m : mlib) {
        auto spl = base::split_str_by_delimiter(m->name(),".");
        if(spl.back() == "0")  {
            m->name("HELIX.IDEAL");
        } else {
            m->name("HELIX.IDEAL." + spl.back());
        }

        auto aligned_motif = mf.align_motif_to_common_frame(m,0);
        
        auto entry = Strings{
                aligned_motif->to_str(), aligned_motif->name(), aligned_motif->ends()[0]->name(), aligned_motif->end_ids()[0],std::to_string(count)
        };

        resources::sqlite3_escape(entry);
        data.push_back(entry);

        aligned_motif = mf.can_align_motif_to_end(m,1);
        aligned_motif = mf.align_motif_to_common_frame(aligned_motif,0);
        mf._setup_secondary_structure(aligned_motif);

        auto reversed_entry = Strings{
                aligned_motif->to_str(),aligned_motif->name(),aligned_motif->ends()[0]->name(),aligned_motif->end_ids()[0],std::to_string(count)
        };

        resources::sqlite3_escape(reversed_entry);
        rdata.push_back(reversed_entry);
        ++count;
    }

    const auto path = parameters_.output_dir +"/motif_libraries_new/ideal_helices.db";
    resources::build_sqlite_library(path,data,motif_keys_,"id");
    const auto reversed_path = parameters_.output_dir +"/motif_libraries_new/ideal_helices_reversed.db";
    resources::build_sqlite_library(reversed_path,rdata,motif_keys_,"id");
    auto sql_lib = resources::MotifSqliteLibrary{1,
                             parameters_.output_dir + resources::MotifSqliteLibrary::get_libnames().at("ideal_helices")};
    auto motif = sql_lib.get("HELIX.IDEAL.3");

    auto outfile = std::ofstream(parameters_.output_dir + "/motifs/base.motif");
    outfile<<motif->to_str();
    outfile.close();
}

void
BuildSqliteLibraries::_build_unique_twoway_library() {
    auto valid_dirs = _get_valid_dirs(parameters_.input_dir + "/two_ways/");
    auto mlib = resources::MotifSqliteLibrary(1,
                              parameters_.output_dir + resources::MotifSqliteLibrary::get_libnames().at("twoway"));
    mlib.load_all();
    auto motifs = motif::MotifOPs{};
    for(auto& m : mlib) {
        motifs.push_back(m);
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
    auto mlib = resources::MotifSqliteLibrary{1,
                     parameters_.output_dir + resources::MotifSqliteLibrary::get_libnames().at("bp_steps")};
    mlib.load_all();

    auto motifs = motif::MotifOPs{};
    auto motif_data = std::vector<Strings>{};
    auto count(0);
    for(auto& m : mlib)  {
        auto spl = base::split_str_by_delimiter(m->name(),".");
        if(spl.size() == 1) {
            motifs.push_back(m); // this branch is never used... pretty sure this logic fro plython doesn't hold up CJ 09/20
        }
    }
    
    auto i(0);
    auto unique = std::unordered_set<String>();
    auto unique_m = std::vector<Strings>();
    auto mf = motif::MotifFactory{};

    for(auto& m : motifs) {
        auto name_spl = base::split_str_by_delimiter(m->name(),"=");
        

        if(unique.find(m->end_ids()[0]) != unique.end() ) {
            continue;
        }
        if (unique.find(m->end_ids()[1]) != unique.end() ) {
            continue;
        }
        
        const auto old_name = m->name();
        unique.insert(m->end_ids()[0]);
        
        if (unique.find(m->end_ids()[1]) == unique.end()) {
            unique.insert(m->end_ids()[1]);
        }
    
        m->name("BP." + std::to_string(i));
        
        auto m_a = mf.can_align_motif_to_end(m,1);
        m_a = mf.can_align_motif_to_end(m_a,1);

        auto entry = Strings{
                m->to_str(), m->name(), m->ends()[0]->name(), m->end_ids()[0], std::to_string(count)
        };
        ++count;
        resources::sqlite3_escape(entry);
        motif_data.push_back(entry);

        auto aligned_entry = Strings{
                m_a->to_str(), m_a->name(), m_a->ends()[0]->name(), m_a->end_ids()[0],std::to_string(count)
        };
        ++count;
        resources::sqlite3_escape(aligned_entry) ;
        unique_m.push_back(aligned_entry);

        ++i;
     }

    const auto path = parameters_.output_dir +"/motif_libraries_new/new_bp_steps.db";

    resources::build_sqlite_library(path, motif_data, motif_keys_, "id");

}

void
BuildSqliteLibraries::_build_flex_helix_library() {
    LOGW<<"WARNING: BuildSqliteLibraries::_build_flex_helix_library() causes the program to crash "
          "from an exception being thrown motif::Motif().. not continuing...";
    return;
    auto helix_files = std::vector<std::filesystem::path>();

    for(const auto& dir : std::filesystem::directory_iterator(base::resources_path() + "/setup/flex_helices/",
                                                              std::filesystem::directory_options::skip_permission_denied)) {
        const auto& path = dir.path();
        if(path.extension() == ".dat" ) {
            helix_files.push_back(path);
        }
    }

    auto flex_helices = motif::MotifOPs{};

    for(const auto& fname : helix_files) {
        const auto name = base::filename(fname);
        const auto length = base::split_str_by_delimiter(name,"_")[1];
        const auto lines = base::get_lines_from_file(fname);
        for(int i = 0; i < lines.size(); ++i) {
            try {
                auto motif = std::make_shared<motif::Motif>(
                        lines[i],
                        structure::ResidueTypeSetManager::getInstance().residue_type_set()
                );
                flex_helices.push_back(motif);
            } catch(std::runtime_error& error) {
                LOGW<<"WARNING: Unable to make motif from file in"
                      " BuildSqliteLibraries::_build_flex_helix_library. Error: "<<error.what()<<"...";
                continue;
            }
        }
    }

    auto successes = std::vector<std::pair<motif::MotifOP,int>>{};
    auto mf = motif::MotifFactory{};
    for(auto& m : flex_helices) {
        auto m_added = mf.can_align_motif_to_end(m, 0);

        if(m_added == nullptr) {
            continue;
        } else {
            successes.emplace_back(m_added, 0);
        }
    }

    auto data = std::vector<Strings>{};

    for(auto i = 0; i < successes.size(); ++i) {
        auto& m = successes[i].first;
        auto& ei = successes[i].second;
        auto m_added = mf.align_motif_to_common_frame(m, ei);

        auto entry = Strings {
            m_added->to_str(), m_added->name(), m_added->ends()[0]->name(), m_added->end_ids()[0],std::to_string(i)
        };

        resources::sqlite3_escape(entry);
        data.push_back(entry);
    }

    const auto path = parameters_.output_dir + "/motif_libraries_new/flex_helices.db";
    resources::build_sqlite_library(path, data, motif_keys_, "id");

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
    const auto existing_path = base::resources_path() + "/setup/existing.motifs";
    if(!base::file_exists(existing_path)) {
        LOGW<<"The file \"existing.motifs\"does not exist... unable to build_existing_motif_library";
        return;
    }
    const auto lines = base::get_lines_from_file(existing_path);
    auto motifs = motif::MotifOPs{};

    for(const auto& line : lines) {
        try {
            auto m = std::make_shared<motif::Motif>(line,
                                                    structure::ResidueTypeSetManager::getInstance().residue_type_set());
            motifs.push_back(std::move(m));
        } catch(std::runtime_error& error) {
            LOGW<<"WARNING: Could not create motif from line. Error: "<<error.what()<<"...";
        }
    }

    auto successes = std::vector<std::pair<motif::MotifOP, int>>{};
    auto mf = motif::MotifFactory{};
    for(auto& m : motifs) {
        if(m->ends().size() == 2) {
            m->mtype(util::MotifType::TWOWAY);
        } else {
            m->mtype(util::MotifType::NWAY);
        }

        for(auto ei = 0; ei < m->ends().size(); ++ei) {
            auto m_added = mf.can_align_motif_to_end(m, ei);
            if(m_added == nullptr){
                continue;
            }

            successes.emplace_back(m_added, ei);
        }
    }
    auto data = std::vector<Strings>{};
    auto removals(0);
    for(auto i = 0; i < successes.size(); ++i) {
        auto& s = successes[i];
        auto m = s.first;
        auto ei = s.second;

        auto m_added = mf.align_motif_to_common_frame(m, ei);

        auto remove = secondary_structure::BasepairOPs{};

        for(auto it = m_added->UNSAFE_basepairs().begin() ;
            it != m_added->UNSAFE_basepairs().end(); ) {
            auto bp = *it;
            auto in_end(false);
            for(int i = 0; i < m_added->ends().size(); ++i) {
                if(*m_added->ends()[i] == *bp) {
                    in_end = true;
                    break;
                }
            }

            if(!in_end && bp->bp_type() == "cW-W") {
                ++removals;
                it = m_added->UNSAFE_basepairs().erase(it);
            } else {
                ++it;
            }

        }
        mf._setup_secondary_structure(m_added);

        auto entry = Strings{
            m_added->to_str(), m_added->name(), m_added->ends()[0]->name(), m_added->end_ids()[0], std::to_string(i)
        };

        resources::sqlite3_escape(entry);
        data.push_back(entry);
    }
    const auto path = parameters_.output_dir + "/motif_libraries_new/existing.db";
    resources::build_sqlite_library(path, data, motif_keys_, "id");
}

void
BuildSqliteLibraries::_build_basic_libraries() {
    const auto types = util::MotifTypes {
        util::MotifType::TWOWAY,
        util::MotifType::NWAY,
        util::MotifType::HAIRPIN,
        util::MotifType::TCONTACT
    };
    const auto type_dirs = std::map<util::MotifType,std::filesystem::path>{
            {util::MotifType::TWOWAY, "/two_ways/"},
            {util::MotifType::NWAY, "/junctions/"},
            {util::MotifType::HAIRPIN, "/hairpins/"},
            {util::MotifType::TCONTACT, "/tertiary_contacts/"}
    };

    const auto bad_keys = std::unordered_set<String>{
            "TWOWAY.2GDI.4-X20-X45",
            "TWOWAY.1S72.46-02097-02647",
            "TWOWAY.2GDI.6-Y20-Y45",
    };

    auto lines = Strings{};
    char* rnamake = std::getenv("RNAMAKE");

    if(rnamake == nullptr) {
        LOGW<<"WARNING: The variable RNAMAKE is not set. Please set it.";
    } else {
        const String csv_path = String(rnamake) + String("/cmake/build/motifs_extra_bps.csv");
        if(std::filesystem::exists(csv_path)) {

            lines = base::get_lines_from_file(csv_path);
            lines.erase(lines.begin(),lines.begin()+1);
            for(auto it = lines.begin() ; it != lines.end();) {
                    if(it->empty()) {
                        it = lines.erase(it);
                    } else {
                        ++it;
                    }
            }

        } else {
            LOGW<<"WARNING: The file "<<rnamake <<
                "cmake/build/motifs_extra_bps.csv does not exist. This will limit what can be done in BuildSqliteLibraries::_build_basic_libraries().";
        }

    }

    auto has_extra_bps = std::map<String,int>{};
    for(const auto& l : lines ) {
        auto spl = base::split_str_by_delimiter(l,",");
        if(spl.size() > 1 && std::stof(spl[2]) > 0 ) {
            has_extra_bps[spl[0]] = 1;
        }
    }

    auto libs = std::map<util::MotifType, motif::MotifOPs>{};
    auto mf = motif::MotifFactory{};
    auto folders = StringStringMap{};
    //TODO change this.. a bit hacky but it works for now CJ 09/20
    auto path = std::filesystem::recursive_directory_iterator(parameters_.input_dir,
                                                               std::filesystem::directory_options::skip_permission_denied);

    for(const auto& t : types) {
        auto mlib = motif::MotifOPs{};
        for(auto& dir : _get_valid_dirs(parameters_.input_dir + type_dirs.at(t).string())) {
            try{
                auto m = mf.motif_from_file(dir.second);
                m->mtype(t);
                mlib.push_back(m);

            } catch(std::runtime_error& error ) {
                LOGW<<"WARNING: Unable to generate motif from "<<dir.second<<". Error: "<<error.what()<<". Continuing...";
            }
        }

        const auto type_str = util::type_to_str(t);
        auto successes = std::vector<std::pair<motif::MotifOP,int>>{};

        for(auto k = 0; k < mlib.size(); ++k) {
            const auto& m = mlib[k];
            if(m == nullptr) {
                continue;
            }

           if(t == util::MotifType::TWOWAY) {
               if(has_extra_bps.find(m->name()) != has_extra_bps.end()) {
                   continue;
               }
           }

           if(t == util::MotifType::HAIRPIN) {
               m->block_end_add(-1);
           }

           if(t != util::MotifType::HAIRPIN && m->ends().size() == 1) {
               continue;
           } else if (m->ends().size() == 0 ) {
               continue;
           }

           for(auto ei = 0; ei < m->ends().size(); ++ei) {
               const auto& m_added = mf.can_align_motif_to_end(m, ei);

               if(m_added == nullptr) {
                   continue;
               } else {
                   const auto key = m_added->name() + "-" + m_added->ends()[0]->name();
                   if(bad_keys.find(key) != bad_keys.end()) {
                       continue;
                   }
                   successes.emplace_back(m_added, ei);
               }
           }
        }

        auto data = std::vector<Strings>();

        for(auto i = 0; i < successes.size(); ++i) {
            const auto& s = successes[i];
            const auto& m = s.first;
            const auto& ei = s.second;
            auto m_added = mf.align_motif_to_common_frame(m, ei);

            if (t == util::MotifType::TWOWAY) {
                for(const auto& bp : m_added->basepairs()) {
                    for(auto it = m_added->UNSAFE_basepairs().begin();
                        it != m_added->UNSAFE_basepairs().end(); ) {
                            auto bp = *it;
                            auto in_end(false);
                            for(int i = 0; i < m_added->ends().size(); ++i) {
                                if(*m_added->ends()[i] == *bp) {
                                    in_end = true;
                                    break;
                                }
                            }
                            if(!in_end && bp->bp_type() == "cW-W") {
                                it = m_added->UNSAFE_basepairs().erase(it);
                            } else {
                                ++it;
                            }
                        }
                }

                try {
                    mf._setup_secondary_structure(m_added);
                } catch (std::runtime_error& error) {
                    LOGW<<"Motif generation failed. Error: "<<error.what()<<". Continuing...";
                    continue;
                }
            }

            auto data_entry = Strings{
                m_added->to_str(), m_added->name(), m_added->ends()[0]->name(), m_added->end_ids()[0], std::to_string(i)
            };

            resources::sqlite3_escape(data_entry);
            data.push_back(data_entry);

        }

        auto motif_name = util::type_to_str(t);
        std::transform(motif_name.begin(), motif_name.end(), motif_name.begin(), tolower);
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
        if(!std::filesystem::exists(parameters_.output_dir + pair.second)) {
            LOGW<<"WARNING: The file "<<parameters_.output_dir + pair.second<<" does not exist. Skipping "
            <<pair.first<<" motif state library...";
            continue;
         }
        auto motif_lib = resources::MotifSqliteLibrary{1,parameters_.output_dir + pair.second};
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
    const auto mp = resources::MotifSqliteLibrary::get_libnames();
    for(const auto& libname : libbames) {
        auto motif_lib = resources::MotifSqliteLibrary{1, parameters_.output_dir + "/" + mp.at(libname)};
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
                auto motif_ensemble = motif::MotifEnsemble{cluster.end_id,all_motifs, energies};

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

            auto motif_ensemble = motif::MotifEnsemble{cluster.end_id, clustered_motifs, energies};

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

bool
BuildSqliteLibraries::_db_contains_motif(resources::MotifSqliteLibrary& lib, motif::MotifOP const& motif) {
    if(motif->end_ids().empty()) {
        return false;
    }
    if(lib.contains(motif->name(), motif->end_ids()[0])) {
        return true;
    }  else  {
        return motif->end_ids().size() > 1 && lib.contains(motif->name(), motif->end_ids()[1]);
    }
}

bool
BuildSqliteLibraries::_db_contains_aligned_motif(resources::MotifSqliteLibrary& lib, motif::MotifOP const& motif) const {
    if(motif->end_ids().empty()) {
        return false;
    }

    if(lib.contains(motif->name(),motif->end_ids()[0]) ) {
        auto new_motif = lib.get(motif->name(), motif->end_ids()[0]);
        return new_motif != nullptr &&
               new_motif->ends()[0]->state()->diff(motif->ends()[0]->state()) < parameters_.motif_distance;

    } else if(motif->end_ids().size() > 1 && lib.contains(motif->name(), motif->end_ids()[1])) {
        auto new_motif = lib.get(motif->name(), motif->end_ids()[1]);
        return new_motif != nullptr &&
               new_motif->ends()[0]->state()->diff(motif->ends()[1]->state()) < parameters_.motif_distance;
    } else {
        return false;
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
    auto sql_lib = resources::MotifSqliteLibrary{1,
                         parameters_.output_dir + resources::MotifSqliteLibrary::get_libnames().at("ideal_helices")};
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

    auto mf = motif::MotifFactory{};
    auto clusters = SSandSeqClusters{};
    constexpr auto kB = 1.3806488e-1;
    constexpr auto kBT = kB * 298.15;
    auto motifs = motif::MotifOPs{};
    for(const auto& pair : _get_valid_dirs(hel_dir.string())) {
        auto m = mf.motif_from_file(pair.second);
        if( m != nullptr ) {
            motifs.push_back(std::move(m));
        } else {
            LOGW<<"WARNING: Unable to extract helix motif from "<<pair.second;
        }
    }

    LOGI<<"Beginning cluster construction. Found "<<motifs.size()<<" motifs to use...";

    for(const auto& m: motifs ) {
        const auto spl = base::split_str_by_delimiter(m->name(),".");

        if(spl.size() > 1 && (spl[1] == "IDEAL" || spl[1] == "LE")) {
            continue;
        }

        for(auto i = 0; i<(m->basepairs().size() - 1); ++i) {
            const auto& basepairs = m->basepairs();
            if(basepairs[i]->bp_type() != "cW-W" || basepairs[i+1]->bp_type() != "cW-W") {
                continue;
            }
            auto bps = structure::BasepairOPs{basepairs[i],basepairs[i+1]};

            auto res = std::vector<structure::Residue>();
            for(const auto& bp : bps)  {
                for(const auto& r : bp->residues()) {
                    if(std::find(res.begin(),res.end(),*r) == res.end()) {
                        res.push_back(*r);
                    }
                }
            }

            if(res.size() != 4) {
                continue;
            }

            auto m_bps = mf.motif_from_bps(bps);

            if(m_bps->end_ids().size() != 2 )  {
                continue;
            }

            if(m_bps->end_ids()[0].size() != 11) {
                continue;
            }

            auto matched(0);
            for( auto& c : clusters) {
                if(c.motif_matches(m_bps)) {
                   matched = 1;
                }
            }

            if( not matched) {
                clusters.emplace_back(m_bps->end_ids()[0]);
                clusters.rbegin()->motif_matches(m_bps);
                if (m_bps->end_ids()[0] != m_bps->end_ids()[1]) {
                    clusters.emplace_back(m_bps->end_ids()[1]);
                    clusters.rbegin()->motif_matches(m_bps);
                }
            }

        }
    }

    LOGI<<"Finished cluster construction for "<<motifs.size()<<" files.";

    auto mes_data = std::vector<Strings>{};
    auto all_mes_data = std::vector<String>{};
    auto all_count(1);
    auto motif_data = std::vector<Strings>{};
    auto unique_data = std::vector<Strings>{};
    auto unique = std::unordered_set<String>{};
    auto bp_count(1);

    auto outfile = std::fstream("bp_step_info.csv");
    outfile<<"name,bp_type,pop,d,r\n";

    auto count(1);

    LOGI<<"Beginning cluster calculations. There are "<<clusters.size()<<" clusters...";

    for(auto c_i=0; c_i<clusters.size(); ++c_i) {

        const auto &c = clusters[c_i];
        auto &m = c.motif_and_ends[0].motif;

        if (unique.find(m->end_ids()[0]) != unique.end()) {
            continue;
        }
        if(unique.find(m->end_ids()[1]) != unique.end()) {
            continue;
        }

        unique.insert(m->end_ids()[0]);
        if (unique.find(m->end_ids()[1]) == unique.end())  {
            unique.insert(m->end_ids()[1]);
        }

        ++all_count;

        auto spl = base::split_str_by_delimiter(c.end_id, "_");

        auto aligned_motifs = motif::MotifOPs{};

        for (auto i = 0; i < c.motif_and_ends.size(); ++i) {

            const auto &m_and_e = c.motif_and_ends[i];
            const auto &m = m_and_e.motif;
            const auto &ei = m_and_e.end_index;

            auto m_a = motif::MotifOP{};
            try {
                m_a = mf.can_align_motif_to_end(m, ei);
            } catch (std::runtime_error &error) {
                LOGW << "Error: " << error.what();
                continue;
            }

            if (m_a == nullptr) {
                continue;
            }

            m_a = mf.align_motif_to_common_frame(m_a, ei);
            auto fail(0);

            for (const auto &bp : m_a->basepairs()) {
                // possible area to check things out
                const auto vec1 = bp->res1()->get_atom("C2'")->coords() - bp->res1()->get_atom("O4'")->coords();
                const auto vec2 = bp->res2()->get_atom("C2'")->coords() - bp->res2()->get_atom("O4'")->coords();

                if (vec1.dot(vec2) > 2.f) {
                    fail = 1;
                    break;
                }
            }

            if (fail) {
                continue;
            }
            aligned_motifs.push_back(m_a);
        }

        auto m_clusters = cluster_motifs(aligned_motifs, 0.65f);

        auto clustered_motifs = motif::MotifOPs{};
        auto energies = std::vector<float>{};
        auto ss = std::stringstream{};
        ss<<spl[0][0]<<spl[2][1]<<"="<<spl[0][1]<<spl[2][0];
        String dir_name = ss.str();

        for (auto j = 0; j < m_clusters.size(); ++j) {
            auto& c_motifs = m_clusters[j].motifs;
            auto m = c_motifs[0];
            m->mtype(util::MotifType::HELIX);
            m->name("BP." + std::to_string(c_i) + "." + std::to_string(j));

            auto motif_entry = Strings{
                    m->to_str(), m->name(), m->ends()[0]->name(), c.end_id,
                    std::to_string(count)
            };
            ++count;

            resources::sqlite3_escape(motif_entry);
            motif_data.push_back(motif_entry);

            clustered_motifs.push_back(m);

            const float pop = float(c_motifs.size()) / float(aligned_motifs.size());

            outfile << m->name() << "," << m->end_ids()[0] << "," << pop << ","
                    << m->ends()[1]->d().to_str();
            outfile << "," << m->ends()[1]->r() << "\n";

            energies.push_back(-kBT * std::log(pop)); // is this the same log?
        }

        if(clustered_motifs.empty()) {
            continue;
        }

        auto me = motif::MotifEnsemble{c.end_id, clustered_motifs, energies};
        auto motif = std::make_shared<motif::Motif>(*me.members()[0]->motif);

        auto pop = std::exp(me.members()[0]->energy / -kBT);
        motif->name("BP." + std::to_string(c_i));
        motif->to_pdb("BP." + std::to_string(c_i) + ".pdb");

        outfile << motif->name() << "," << motif->end_ids()[0] << "," << std::to_string(pop) << ","
                << motif->ends()[1]->d().to_str();
        outfile << "," << motif->ends()[1]->r().to_str() << "\n";

        auto mes_entry = Strings{
                me.to_str(),me.id(),std::to_string(bp_count)
        };
        resources::sqlite3_escape(mes_entry);

        mes_data.push_back(mes_entry);

        auto unique_data_entry = Strings{
                motif->to_str(),motif->name(),motif->ends()[0]->name(), me.id(), std::to_string(bp_count)
        };
        ++bp_count;

        resources::sqlite3_escape(unique_data_entry);
        unique_data.push_back(unique_data_entry);

        clustered_motifs.clear();
        energies.clear();

        for (const auto &mem : me.members()) {

            auto m_a = mf.can_align_motif_to_end(mem->motif, 1);
            m_a = mf.align_motif_to_common_frame(m_a, 1);
            clustered_motifs.push_back(m_a);
            energies.push_back(mem->energy);
            pop = std::exp(mem->energy / -kBT);

            if (motif->end_ids()[0] != motif->end_ids()[1]) {
                outfile << m_a->name() << "," << m_a->end_ids()[0] << "," << pop << ",";
                outfile << m_a->ends()[1]->d().to_str() << "," << m_a->ends()[1]->r().to_str()
                        << "\n";
            }

            auto motif_data_entry = Strings{
                    m_a->to_str(),m_a->name(),m_a->ends()[0]->name(),m_a->end_ids()[0],std::to_string(count)
            };
            ++count;
            resources::sqlite3_escape(motif_data_entry);
            motif_data.push_back(motif_data_entry);
        }

        me = motif::MotifEnsemble{clustered_motifs[0]->end_ids()[0],clustered_motifs,energies};

        if(motif->end_ids()[0] != motif->end_ids()[1]) {
            pop = std::exp((*me.members().rbegin())->energy/-kBT);

            auto mes_data_entry = Strings{me.to_str(),me.id(),std::to_string(bp_count)};

            resources::sqlite3_escape(mes_data_entry);

            motif = me.members()[0]->motif; // not sure this is the right way to do it?? CJ 09/20

            motif->name(String{"BP."} + std::to_string(c_i));
            outfile<<motif->name()<<","<<motif->end_ids()[0]<<","<<pop;
            outfile<<motif->ends()[1]->d().to_str()<<","<<motif->ends()[1]->r().to_str()<<std::endl;
        }

        auto unique_entry = Strings{
            motif->to_str(),motif->name(),motif->ends()[0]->name(),me.id(),std::to_string(bp_count)
        };
        ++bp_count;
        resources::sqlite3_escape(unique_entry);
        unique_data.push_back(unique_entry);

    }

    LOGI<<"Finished cluster calculations for "<<clusters.size()<<" clusters ";

    String path = parameters_.output_dir + "/motif_ensemble_libraries_new/bp_steps.db";
    resources::build_sqlite_library(path, mes_data, ensemble_keys_, "id");

    path = parameters_.output_dir + "/motif_libraries_new/bp_steps.db";
    resources::build_sqlite_library(path, motif_data, motif_keys_, "id");

    path = parameters_.output_dir + "/motif_libraries_new/new_bp_steps.db";
    resources::build_sqlite_library(path, unique_data, motif_keys_, "id");
}

void
BuildSqliteLibraries::setup_options()   {

    auto build = app_.add_subcommand("BUILD");

    build->add_flag("--all",parameters_.all,"flag to build all available libraries")
            ->group("Core Inputs");

    build->add_option("--in_dir",parameters_.input_dir,"directory that contains input data folders")
            ->default_val(base::resources_path())
            ->check(CLI::ExistingDirectory)
            ->group("Core Inputs");

    build->add_option("--out_dir",parameters_.output_dir,"directory where the files will be published to")
            ->required()
            ->group("Core Inputs");

    build->add_option("--log_level",parameters_.log_level,"level for global logging")
            ->check(CLI::IsMember(std::set<String>{"debug","error","fatal","info","verbose","warn"}))
            ->default_val("info")
            ->group("Core Inputs");

    for(auto& pair : build_opts_ )  {
        build->add_flag("--" + pair.first,pair.second.selected)
                ->group("Build Options");
    }

    auto validate = app_.add_subcommand("VALIDATE");

    validate->add_option("--orig",parameters_.orig_dir)
            ->default_val(base::resources_path())
            ->check(CLI::ExistingDirectory);

    validate->add_option("--remake",parameters_.remake_dir)
            ->required()
            ->check(CLI::ExistingDirectory);

    validate->add_option("--motif_dist", parameters_.motif_distance)
            ->default_val(0.5f)
            ->description("The distance that basepair end states must be less than to satisfy a motif match")
            ->check(CLI::PositiveNumber);

    validate->add_option("--summary_file", parameters_.summary_file)
            ->description("The .csv file to where the comparison results will be written");

}

void
BuildSqliteLibraries::_validate() {

    LOGI<<"Beginning sqlite database file comparision...";
    orig_files_ = _get_paths(parameters_.orig_dir);
    remake_files_ = _get_paths(parameters_.remake_dir);

    auto mp = resources::MotifSqliteLibrary::get_libnames();

    LOGI<<"Beginning phase I: Checking files that exist in the original db set...";
    for(auto it = mp.begin(); it != mp.end() ; ) {
        const auto db_path = parameters_.orig_dir.string() + it->second ;
        if(std::filesystem::exists(db_path)) {
            LOGI<<"\tThe file "<<db_path<<" exists in the original db set...";
            ++it;
        } else {
            LOGW<<"\tWARNING: The file "<<db_path<<
            " does NOT exist in the original db set. Removing it from further analysis";
            it = mp.erase(it);
        }
    }
    LOGI<<"There are "<<mp.size()<<" files in the original db set...";
    LOGI<<"Checking which files exist in the remake dir "<<parameters_.remake_dir.string()<<"...";

    auto bad(0);
    for(auto it = mp.begin(); it != mp.end(); ) {
        if( std::filesystem::exists(parameters_.remake_dir.string() + it->second))  {
            ++it;
        } else {
            LOGW<<"WARNING: The file "<<parameters_.remake_dir.string() + it->second
                <<" does not exist, so no further analysis will be done on "<<it->second<<" in the original set";
            ++bad;
            it = mp.erase( it );
        }
    }
    LOGI<<"Finished phase I! "<<mp.size()<<" files exist, "<<bad<<" don't exist";
    LOGI<<"Beginning phase II: Making sure that the old motifs are a subset of the new ones...";
    if(parameters_.summary_file != "NONE") {
        writer_.open(parameters_.summary_file,std::ios::out);
        writer_<<"category,motif_dist,matches,misses,aligned,misaligned,missed_names,misaligned_names\n";
    }
    for(const auto& pr : mp) {
        LOGI<<"\tBeginning "<<pr.first<<" comparison...";
        if(parameters_.summary_file != "NONE") {
            writer_<<pr.first<<","<<parameters_.motif_distance<<",";
        }
        auto missed_names = String{};
        auto misaligned_names = String{};
        auto orig_lib = resources::MotifSqliteLibrary(1, parameters_.orig_dir.string() + pr.second);
        orig_lib.load_all();
        auto new_lib = resources::MotifSqliteLibrary(1,parameters_.remake_dir.string() + pr.second);
        new_lib.load_all();
        auto match(0), miss(0), aligned(0), misaligned(0);
        for(const auto& orig_motif : orig_lib)  {
            if(_db_contains_motif(new_lib, orig_motif)) {
                ++match;
            } else {
                missed_names += orig_motif->name() + "|" ;
                ++miss;
            }
            if(_db_contains_aligned_motif(new_lib, orig_motif)) {
                ++aligned;
            } else {
                misaligned_names += orig_motif->name() + "|" ;
                ++misaligned;
            }
        }
        LOGI<<"\tFinished "<<pr.first<<" comparison. Matches: "<<match<<", Misses: "<<miss
        <<". Aligned within distance "<<parameters_.motif_distance<<": "<<aligned<<", Misaligned: "<<misaligned<<".";

        if(!missed_names.empty()) {
            missed_names.pop_back();
        } else {
            missed_names = "NA";
        }

        if(!misaligned_names.empty()) {
            misaligned_names.pop_back();
        } else {
            misaligned_names = "NA";
        }

        if(parameters_.summary_file != "NONE") {
            writer_<<match<<","<<miss<<","<<aligned<<","<<misaligned<<","<<missed_names<<","<<misaligned_names<<"\n";
        }
    }

    LOGI<<"Finished phase II. "<<
    (parameters_.summary_file == "NONE" ?
            "Results not written to a file. Exiting..." :
            "Results written to: " + parameters_.summary_file.string());

}

void
BuildSqliteLibraries::run() {

    base::init_logging(log_level());

    if(app_.got_subcommand("BUILD")) {
        _build();
    } else if (app_.got_subcommand("VALIDATE")) {
       _validate();
    } else {
        LOGF<<"You must enter a subcommand. Options are BUILD or VALIDATE. Exiting...";
        exit(1);
    }

}

int main(int argc, const char** argv) {
    /* NOTES: existing.motifs does not exist... did not ressurrect this. Should ? Therefore, this part does not match
     * */
    std::set_terminate(base::print_backtrace);

    auto env_mgr = base::EnvManager(Strings{"RNAMAKE"});    
    env_mgr.set_envs();
    auto app = BuildSqliteLibraries{};

    app.setup_options();
    CLI11_PARSE(app.app_,argc,argv);
    app.run();
}
