#include <build_sqlite_libraries/build_sqlite_libraries.h>

#define MOTIF_VERSION "/motif_state_librariesV2/"

struct MotifandEnd{

};

using MotifandEnds = std::vector<MotifandEnd>;

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
    const auto bad_keys = Strings{
            "TWOWAY.2GDI.4-X20-X45",
            "TWOWAY.1S72.46-02097-02647",
            "TWOWAY.2GDI.6-Y20-Y45"
    };
    /*j
    auto file_lines = base::get_lines_from_file("motifs_extra_bps.csv");
    file_lines.erase(file_lines.begin(),file_lines.begin()+1);


    auto has_extra_bps = std::map<String,int>{};
    for(const auto& line : file_lines ) {
        auto tokens = base::split_str_by_delimiter(line,",");
        if(tokens.size() > 1 && std::stof(tokens[2]) > 0 ) {
            has_extra_bps[tokens[0]] = 1;
        }
    }
*/
    for(auto&& type : motif_types) {
        std::cout<<util::type_to_str(type)<<std::endl;
    }
}

void
BuildSqliteLibraries::run()   {
//TODO


//build_existing_motif_library();


//    builder.build_helix_ensembles()
//#builder.build_flex_helix_library()
//#builder.build_ss_and_seq_libraries()
//    builder.build_motif_state_libraries()
//    builder.build_motif_ensemble_state_libraries()
//#builder.build_flex_helices()
//

//DONE
            build_basic_libraries();
            //build_unique_twoway_library();
            //build_new_bp_steps();
            //build_ideal_helices();
            //build_trimmed_ideal_helix_library();
    //build_helix_ensembles();
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
    const auto allowed_motifs = std::set<String>{
                                                 "basepair_steps",
                                                 "extras",
                                                 "hairpins",
                                                 "helices",
                                                 "junctions",
                                                 "messed_up",
                                                 "tertiary_contact_hairpin_hairpin",
                                                 "tertiary_contact_junction_hairpin",
                                                 "tertiary_contact_junction_junction",
                                                 "tertiary_contact_nway_hairpin",
                                                 "tertiary_contacts",
                                                 "two_ways"
                                                };
    const auto motif_type = options_.get_string("motif_type");

    if(allowed_motifs.find(motif_type) == allowed_motifs.end()) {
        LOG_ERROR<<"The motif type \""<<motif_type<<"\" is not allowed";
        LOG_ERROR<<"Allowed types include: "; 
        
        for(const auto& type : allowed_motifs) { 
            LOG_ERROR<<"\t"<<type;
        }

        LOG_ERROR<<"Now exiting...";
        exit(1);
    }
    
    auto motif_paths = _get_valid_dirs(options_.get_string("dir") + "/" + options_.get_string("motif_type"));
    
    auto mf = motif::MotifFactory{};
    const auto total = motif_paths.size();
    auto i(0);
    for(auto&& [motif,path] : motif_paths) {
        std::cout<<++i<<"\tof\t"<<total<<std::endl; 
        motif_map_[motif] = mf.motif_from_file(path);
        if(i==5) break;
    }
    //auto clusters =  SSandSeqClusters{};
    for(auto&& [name,motif] : motif_map_) {
        auto tks = base::split_str_by_delimiter(name,".");
        
        if(tks[1] == "IDEAL" || tks[1] == "LE") {
            continue;
        }
        const auto& bps = motif->basepairs();
        const auto bp_num = bps.size() - 1;
        
        for(auto ii = 0; ii<bp_num; ++ii) {
            if(bps[ii]->bp_type() != "cW-W" || bps[ii+1]->bp_type() != "cW-W") {
                continue;
            }
            auto basepair = {bps[ii],bps[ii+1]};
            auto residues = structure::ResidueOPs{};
            auto res_names = std::set<uint64_t>{};

            for(auto& bp : basepair) {
                for(auto& res : bp->residues()) {
                    const auto address((long long)(res.get()));
                    if(res_names.find(address) == res_names.end()) {
                        res_names.insert(address);
                        residues.push_back(res);
                    }
                }
            }

            if(residues.size() != 4) {
                continue;
            }

            auto m_bps = mf.motif_from_bps(bps);

            if( m_bps->end_ids().size() != 2) {
                continue;
            }

            if(m_bps->end_ids()[0].size() != 11) {
                continue; //TODO ask joe if this is correct
            }

            auto matched(false);

        }

    }

    //for(auto& c : lib_names_) std::cout<<c.second<<std::endl; 
    //for(const auto& hel_type : {"avg_helices","flex_helices","ideal_helices","ideal_helices_reversed", "le_helices"} ) {
    //    auto motif_lib = resources::MotifSqliteLibrary(hel_type); 
    //    motif_lib.load_all(); 

    //    try{
    //        std::cout<<motif_lib.get_multi().size()<<std::endl;
    //        //for(auto& motif : ) {


    //        //}

    //    } catch (std::runtime_error& error) {
    //        std::cout<<error.what()<<std::endl;
    //    }

    //}

    //for(const auto& pr : lib_names_) {
    //    std::cout<<pr.first<<"\t"<<pr.second<<"\n";
    //}
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


