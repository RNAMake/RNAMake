#include <build_sqlite_libraries/build_sqlite_libraries.h>

#define MOTIF_VERSION "/motif_state_librariesV2/"

struct MotifEnsembleMember{
    motif::MotifOP motif;
    float energy;

    MotifEnsembleMember(motif::MotifOP const& motif, float energy) :
                                            motif(motif), energy(energy) {}

    MotifEnsembleMember(MotifEnsembleMember const & other ) :
                                            motif(other.motif), energy(other.energy) {}

    String
    to_str() const {
        return motif->to_str() + "#" + std::to_string(energy); // Not sure how reliable std::to_string is for C++ vs ptyhon... -CJ 08/20
    }

};

using MotifEnsembleMembers = std::vector<MotifEnsembleMember>;

struct MotifEnsemble{

    String id;
    MotifEnsembleMembers members;
    int block_and_end;

    MotifEnsemble() {

    }

    MotifEnsemble(MotifEnsemble const& other ) :
                        id(other.id), members(other.members), block_and_end(other.block_and_end) {}

    void
    setup(String const& id, motif::MotifOPs const& motifs, std::vector<float> const& energies) {

        assert(motifs.size() == energies.size());

        this->id = id;
        this->members.clear();

        const auto num_motifs = motifs.size();

        for(auto ii = 0; ii<num_motifs; ++ii) {
            members.emplace_back(motifs[ii],energies[ii]);
        }

        std::sort(members.begin(), members.end(),
                [] (const auto& m1, const auto& m2) {
            return m1.energy < m2.energy; // I think this is right? -CJ -8/20
        }
        );
    }

    String
    to_str() {
        auto ss = std::stringstream{};

        ss<<id<<"{"<<block_and_end<<"{";
        for(const auto& ms : members) {
            ss<<ms.to_str()<<"{";
        }

        return ss.str();
    }
};

struct MotifandEnd{
    motif::MotifOP motif;
    int index;
    MotifandEnd(motif::MotifOP const& motif, int index) : motif(motif), index(index) {

    }
};

using MotifandEnds = std::vector<MotifandEnd>;

class SSandSeqCluster {
public:
    String id;
    MotifandEnds motif_and_ends;
public:
    SSandSeqCluster(String const& id) : id(id), motif_and_ends(MotifandEnds{}) {

    }

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
    // everything is a string, my dude
    String end_id;
    MotifandEnds motif_and_ends;

    // methods
    //bool //motif_matches_end(m,ei) { //} }; using SSandSeqClusters = std::vector<SSandSeqCluster>; struct MotifCluster{
    MotifCluster(motif::MotifOP const& motif, structure::BasepairStateOP state) : state(std::move(state))
    {
        motifs.push_back(motif);
    }

    structure::BasepairStateOP state; 
    motif::MotifOPs motifs;

};

using MotifClusters = std::vector<MotifCluster>;

MotifClusters
cluster_motifs(resources::MotifSqliteLibrary::iterator && start,resources::MotifSqliteLibrary::iterator const & end, float max_distance=1.5f) {
    
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

void
build_trimmed_ideal_helix_library();

StringStringMap
BuildSqliteLibraries::_get_valid_dirs(String const& base_dir) {
    
    auto paths = StringStringMap{};

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
    
    return paths;
}

void
BuildSqliteLibraries::build_ideal_helices() {
    auto paths = _get_valid_dirs(options_.get_string("dir") + "/" + options_.get_string("motif_type"));
    auto mf = motif::MotifFactory{}; 
    const auto keys = Strings{"data","name","end_name","end_id","id"}; 
    auto data = std::vector<Strings>{};
    auto reversed_data = std::vector<Strings>{};
    auto motif_ct(0); 

    for(auto& [folder, path] : paths) {
        auto tokens = base::split_str_by_delimiter(folder,".");
        
        if(tokens[1] != "LE") {
            continue;
        }

        auto motif = mf.motif_from_file(path);
        
        const auto last_dot = folder.find_last_of('.') + 1;
        const auto num = folder.substr(last_dot);
        motif->name(
            "HELIX.IDEAL" + (num == "0" ? "" : "." + num)
                );

        auto aligned_motif = mf.align_motif_to_common_frame(motif,0); 
        
        auto entry = Strings{};
        entry.reserve(keys.size());
       
        entry.push_back(aligned_motif->to_str());
        entry.push_back(aligned_motif->name());
        entry.push_back(aligned_motif->end_name(0));
        entry.push_back(aligned_motif->end_ids()[0]);
        entry.push_back(std::to_string(motif_ct));

        for(auto& token : entry)  {
            token = base::replace_all(token,"\'","\'\'");
        }
        data.push_back(entry);

        auto reversed_motif = mf.can_align_motif_to_end(aligned_motif,1);
        auto reversed_aligned_motif = mf.align_motif_to_common_frame(reversed_motif,0);
        mf._setup_secondary_structure(reversed_aligned_motif);
        auto reversed_entry = Strings{};
        reversed_entry.reserve(keys.size());
       
        reversed_entry.push_back(reversed_aligned_motif->to_str());
        reversed_entry.push_back(reversed_aligned_motif->name());
        reversed_entry.push_back(reversed_aligned_motif->end_name(0));
        reversed_entry.push_back(reversed_aligned_motif->end_ids()[0]);
        reversed_entry.push_back(std::to_string(motif_ct));
        
        for(auto& token : reversed_entry)  {
            token = base::replace_all(token,"\'","\'\'");
        }

        reversed_data.push_back(reversed_entry);

        ++motif_ct;
    }
    
    const auto path = base::resources_path() +"/motif_librariesV2/ideal_helices.db";
    resources::build_sqlite_library(
                                    path,
                                    data,
                                    keys,
                                    "id"
                                    );

    const auto reversed_path = base::resources_path() +"/motif_librariesV2/ideal_helices_reversed.db";
    resources::build_sqlite_library(
                                    reversed_path,
                                    reversed_data,
                                    keys,
                                    "id"
                                    );
    auto sql_lib = resources::MotifSqliteLibrary{"ideal_helices"};
    auto motif = sql_lib.get("HELIX.IDEAL.3");
    auto outfile = std::ofstream(base::resources_path() + "/motifs/base.motif");
    outfile<<motif->to_str();
    outfile.close();
}

void
BuildSqliteLibraries::build_unique_twoway_library() {
    auto valid_dirs = _get_valid_dirs(base::resources_path() + "/motifsV2/two_ways/");
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
    //auto motif_lib = resources::MotifSqliteLibrary{"twoway"};
    //motif_lib.load_all();
    //auto clusters = cluster_motifs(motif_lib.begin(), motif_lib.end(), 9.0f);

    const auto keys = Strings{"data", "name", "end_name", "end_id", "id"};
    auto data = std::vector<Strings>{};

    const auto mes_keys = Strings{"data", "name", "id"};
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

        auto entry = Strings{};
        entry.reserve(5);

        entry.push_back(lowest->to_str());
        entry.push_back(lowest->name());
        entry.push_back(lowest->ends()[0]->name());
        entry.push_back(lowest->end_ids()[0]);
        entry.push_back(std::to_string(count));

        for(auto& token : entry)  {
            token = base::replace_all(token,"\'","\'\'");
        }

        data.push_back(entry);
    }

    outfile.close();
    const auto path = base::resources_path() + "/motif_librariesV2/unique_twoway.db";

    resources::build_sqlite_library(
                                    path,
                                    data,
                                    keys,
                                    "id");

}

void
BuildSqliteLibraries::build_new_bp_steps() {
    
    auto sql_lib = resources::MotifSqliteLibrary{"bp_steps"};
    sql_lib.load_all();

    const auto keys = Strings{"data","name","end_name","end_id","id"};
    auto motifs = motif::MotifOPs{};
    for(auto& motif : sql_lib)  {
        auto tokens = base::split_str_by_delimiter(motif->name(),".");
        if(tokens.size() == 1) {

            motifs.push_back(motif);
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

        auto entry = Strings{};
        entry.reserve(keys.size());
       
        entry.push_back(motif->to_str());
        entry.push_back(motif->name());
        entry.push_back(motif->end_name(0));
        entry.push_back(motif->end_ids()[0]);
        entry.push_back(std::to_string(motif_ct));
        
        ++motif_ct;

        for(auto& token : entry)  {
            token = base::replace_all(token,"\'","\'\'");
        }

        motif_data.push_back(entry);

        auto aligned_entry = Strings{};
        aligned_entry.reserve(keys.size());
       
        aligned_entry.push_back(motif_aligned->to_str());
        aligned_entry.push_back(motif_aligned->name());
        aligned_entry.push_back(motif_aligned->end_name(0));
        aligned_entry.push_back(motif_aligned->end_ids()[0]);
        aligned_entry.push_back(std::to_string(motif_ct));

        for(auto& token : aligned_entry)  {
            token = base::replace_all(token,"\'","\'\'");
        }
        
        ++motif_ct;

        motif_data.push_back(aligned_entry);
        
        ++ii;

     }


    const auto path = base::resources_path() +"/motif_librariesV2/new_bp_steps.db";
    resources::build_sqlite_library(
                                    path,
                                    motif_data,
                                    keys,
                                    "id"
            );

}

void
BuildSqliteLibraries::build_existing_motif_library() {
    throw std::runtime_error(
        "BuildSqliteLibraries::build_existing_motif_library() has not been implemented yet! Don't call it"
            ) ;
    //auto existing_motifs = Strings{};
    //if(base::file_exists("existing.motifs")) {
    //    existing_motifs = base::get_lines_from_file("existing.motifs");
    //}
    //auto motifs = motif::MotifOPs{};
    //auto mf = motif::MotifFactory{};

    //for(const auto& motif_str : existing_motifs) {
    //    motifs.push_back(
    //        motif::MotifOP{motif_str}
    //            );
    //}
    //for(auto& c : existing_motifs) {
    //    std::cout<<c<<std::endl;
    //}
}
void
BuildSqliteLibraries::build_basic_libraries() {
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

    auto file_lines = base::get_lines_from_file("motifs_extra_bps.csv");
    file_lines.erase(file_lines.begin(),file_lines.begin()+1);

    auto has_extra_bps = std::map<String,int>{};

    for(const auto& line : file_lines ) {
        auto tokens = base::split_str_by_delimiter(line,",");
        if(tokens.size() > 1 && std::stof(tokens[2]) > 0 ) {
            has_extra_bps[tokens[0]] = 1;
        }
    }

    auto folders = StringStringMap{};
    const auto type_dirs = Strings{"../../resources/motifsV2/basepair_steps/",
                                    "../../resources/motifsV2/extras/",
                                    "../../resources/motifsV2/hairpins/",
                                    "../../resources/motifsV2/helices/",
                                    "../../resources/motifsV2/junctions/",
                                    "../../resources/motifsV2/messed_up/",
                                    "../../resources/motifsV2/tertiary_contact_hairpin_hairpin/",
                                    "../../resources/motifsV2/tertiary_contact_junction_hairpin/",
                                    "../../resources/motifsV2/tertiary_contact_junction_junction/",
                                    "../../resources/motifsV2/tertiary_contact_nway_hairpin/",
                                    "../../resources/motifsV2/tertiary_contacts/",
                                    "../../resources/motifsV2/two_ways/"};
    for(auto& dir : type_dirs) {
        auto entries = _get_valid_dirs(dir);
        folders.insert(entries.begin(),entries.end());
    }

    auto mf = motif::MotifFactory{};
    for(const auto& type : motif_types) {
        const auto type_str = util::type_to_str(type);

        //auto mlib = resources::MotifSqliteLibrary(type_str);

        auto successes = std::vector<std::pair<motif::MotifOP,int>>{};

        for(const auto& [dir, file] : folders) {

            if(type_str == "TCONTACT" && dir.find("TC.") == std::string::npos) {
                continue;
            } else if (dir.find(type_str) == std::string::npos) {
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

        const auto keys = Strings{"data","name","end_name","end_id","id"};

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

            auto data_entry = Strings{};
            data_entry.reserve(5);

            data_entry.push_back(m_added->to_str());
            data_entry.push_back(m_added->name()),
            data_entry.push_back(m_added->ends()[0]->name());
            data_entry.push_back(m_added->end_ids()[0]);
            data_entry.push_back(std::to_string(motif_ct++));

            for(auto& token : data_entry)  {
                token = base::replace_all(token,"\'","\'\'");
            }

            data.push_back(std::move(data_entry));

        }

        if(data.empty())  {
            LOGW<<"No data generated for motif type: "<<util::type_to_str(type);
            continue;
        }

        auto motif_name = util::type_to_str(type);
        std::transform(motif_name.begin(), motif_name.end(), motif_name.begin(),tolower);
        const auto path = base::resources_path() + "motif_librariesV2/" + motif_name + ".db";
        resources::build_sqlite_library(
                path,
                data,
                keys,
                "id"
        );
    }

}

void
BuildSqliteLibraries::run()   {
//TODO

//build_ss_and_seq_libraries();


//build_existing_motif_library();


//    builder.build_helix_ensembles()
//#builder.build_flex_helix_library()
//#builder.
//    builder.build_motif_state_libraries()
//    builder.build_motif_ensemble_state_libraries()
//#builder.build_flex_helices()
//

//DONE
            //build_basic_libraries();
            //build_unique_twoway_library();
            //build_new_bp_steps();
            //build_ideal_helices();
            //build_trimmed_ideal_helix_library();
            build_helix_ensembles();
}

void
BuildSqliteLibraries::setup_options()   {
    
//    options.add_argument("-d","--directory")
//            .required();
//   
//    options.add_argument("-m","--motif_type")
//            .required();
//
//    add_option("dir", "", base::OptionType::STRING, true);
//    add_option("motif_type", "", base::OptionType::STRING, true);
}

void
BuildSqliteLibraries::build_trimmed_ideal_helix_library() {

    auto helix_folders = _get_valid_dirs(options_.get_string("dir") + "/helices/" ) ;
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
        auto entry = Strings{};
        entry.reserve(5);

        entry.push_back(m->to_str());
        entry.push_back(m->name());
        entry.push_back(m->end_name(0));
        entry.push_back(m->end_ids()[0]);
        entry.push_back(std::to_string(motif_ct++));

        for(auto& token : entry)  {
            token = base::replace_all(token,"\'","\'\'");
        }

        data.push_back(entry);

    }

    const auto keys = Strings{"data", "name", "end_name", "end_id", "id"};
    const auto path = base::resources_path()  + MOTIF_VERSION + lib_name->first + "_min.db";
    resources::build_sqlite_library(
                                    path,
                                    data,
                                    keys,
                                    "id"
                                    );
}

void
BuildSqliteLibraries::build_helix_ensembles() {

    const auto helix_paths = _get_valid_dirs(base::resources_path() + "/motifsV2/helices/");

    auto clusters = SSandSeqClusters{};
    constexpr auto kB = 1.3806488e-1;
    constexpr auto kBT = kB * 298.15;

    auto mf = motif::MotifFactory{};
    for(const auto& [dir,path] : helix_paths ) {
        const auto tokens = base::split_str_by_delimiter(dir,".");

        if(tokens.size() > 1 && (
                tokens[1] == "IDEAL" || tokens[1] == "LE"
                )) {
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

    const auto mes_keys = Strings{"data","name","id"};
    auto mes_data = std::vector<Strings>{};

    const auto all_mes_keys = Strings{"data","name","id"};
    auto all_mes_data = std::vector<String>{};
    auto all_count(0);

    const auto motif_keys = Strings{"data","name","end_name","end_id","id"};
    auto motif_data = std::vector<Strings>{};

    auto unique_data = std::vector<Strings>{};
    auto unique = std::unordered_set<String>{};
    auto bp_count(1);

    auto outfile = std::fstream("bp_step_info.csv");
    outfile<<"name,bp_type,pop,d\n";

    auto count(1);

    const auto num_clusters = clusters.size();

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

                const auto vec1 = bp->res1()->get_atom("C2\'")->coords() - bp->res1()->get_atom("O4\'")->coords();
                const auto vec2 = bp->res2()->get_atom("C2\'")->coords() - bp->res2()->get_atom("O4\'")->coords();

                if (vec1.dot(vec2) > 2.f) {
                    fail = true;
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
        auto dir_name = String{tokens[0][0]};
        dir_name += tokens[2][1];
        dir_name += "=";
        dir_name += tokens[0][1];
        dir_name + tokens[2][0];

        const auto num_motif_clusters = motif_clusters.size();

        for (auto jj = 0; jj < num_motif_clusters; ++jj) {
            auto cluster_motifs = motif_clusters[jj].motifs;
            auto curr_motif = cluster_motifs[0];
            curr_motif->mtype(util::MotifType::HELIX);
            curr_motif->name("BP." + std::to_string(ii) + "." + std::to_string(jj));

            auto entry = Strings{};
            entry.reserve(5);
            entry.push_back(curr_motif->to_str());
            entry.push_back(curr_motif->name());
            entry.push_back(curr_motif->ends()[0]->name());
            entry.push_back(cluster.id);
            entry.push_back(std::to_string(count++));

            for (auto &token : entry) {
                token = base::replace_all(token, "\'", "\'\'");
            }

            cluster_motifs.push_back(curr_motif);

            const auto pop = float(cluster_motifs.size()) / float(aligned_motifs.size());

            outfile << curr_motif->name() << "," << curr_motif->end_ids()[0] << "," << pop << ","
                    << curr_motif->ends()[1]->d().to_str();
            outfile << "," << curr_motif->ends()[1]->r() << "\n";

            energies.push_back(-kBT * std::log(pop));
        }

        auto motif_ensemble = MotifEnsemble{};

        if(clustered_motifs.empty()) {
            continue;
        }
        std::cout<<"HERE"<<std::endl;
        motif_ensemble.setup(cluster.id, clustered_motifs, energies);
        auto e_motif = motif_ensemble.members[0].motif;
        auto pop = std::exp(motif_ensemble.members[0].energy / -kBT);
        e_motif->name("BP." + std::to_string(ii));
        e_motif->to_pdb("BP." + std::to_string(ii) + ".pdb");

        outfile << e_motif->name() << "," << e_motif->end_ids()[0] << "," << std::to_string(pop) << ","
                << e_motif->ends()[1]->d().to_str();
        outfile << "," << e_motif->ends()[1]->r().to_str() << "\n";

        auto mes_entry = Strings{};
        mes_entry.reserve(3);
        mes_entry.push_back(motif_ensemble.to_str());
        mes_entry.push_back(motif_ensemble.id);
        mes_entry.push_back(std::to_string(bp_count));

        for (auto &token : mes_entry) {
            token = base::replace_all(token, "\'", "\'\'");
        }

        mes_data.push_back(mes_entry);

        auto unique_data_entry = Strings{};
        unique_data_entry.reserve(5);
        unique_data_entry.push_back(e_motif->to_str());
        unique_data_entry.push_back(e_motif->name());
        unique_data_entry.push_back(e_motif->ends()[0]->name());
        unique_data_entry.push_back(motif_ensemble.id);
        unique_data_entry.push_back(std::to_string(bp_count++));

        for (auto &token : unique_data_entry) {
            token = base::replace_all(token, "\'", "\'\'");
        }

        unique_data.push_back(unique_data_entry);


        clustered_motifs.clear();
        energies.clear();

        for (const auto &member : motif_ensemble.members) {

            auto motif_aligned = mf.can_align_motif_to_end(member.motif, 1);
            motif_aligned = mf.can_align_motif_to_end(motif_aligned, 1);
            clustered_motifs.push_back(motif_aligned);
            energies.push_back(member.energy);
            pop = std::exp(member.energy / -kBT);

            if (e_motif->end_ids()[0] != e_motif->end_ids()[1]) {
                outfile << motif_aligned->name() << "," << motif_aligned->end_ids()[0] << "," << pop << ",";
                outfile << motif_aligned->ends()[1]->d().to_str() << "," << motif_aligned->ends()[1]->r().to_str()
                        << "\n";
            }

            auto motif_data_entry = Strings{};
            motif_data_entry.reserve(5);
            motif_data_entry.push_back(motif_aligned->to_str());
            motif_data_entry.push_back(motif_aligned->name());
            motif_data_entry.push_back(motif_aligned->ends()[0]->name());
            motif_data_entry.push_back(motif_aligned->end_ids()[0]);
            motif_data_entry.push_back(std::to_string(count++));

            for (auto &token : motif_data_entry) {
                token = base::replace_all(token, "\'", "\'\'");
            }

            motif_data.push_back(motif_data_entry);
        }

        motif_ensemble = MotifEnsemble{};
        motif_ensemble.setup(clustered_motifs[0]->end_ids()[0],clustered_motifs,energies);

        if(motif->end_ids()[0] != motif->end_ids()[1]) {
            //pop = std::exp()I
        }
    }


}

int main(int argc, const char** argv) {
    
    base::init_logging();
    
    auto env_mgr = base::EnvManager(Strings{"RNAMAKE"});    
    env_mgr.set_envs();

    auto app = BuildSqliteLibraries{};
    app.setup_options();
    app.parse_command_line(argc,argv);
    app.run();

}



