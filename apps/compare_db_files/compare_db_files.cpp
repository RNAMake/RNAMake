#include <algorithm>
#include <compare_db_files/compare_db_files.h>

void
CompareDBFiles::run() {
    base::init_logging(base::log_level_from_str(parameters_.log_level));

    LOGI<<"Beginning sqlite database file comparision...";
    dir_1_files_ = _get_paths(parameters_.dir_1);

    auto mp = resources::MotifSqliteLibrary::get_libnames();

    LOGI<<"Beginning phase I: Checking that files exist at all...";
    auto bad(0);
    for(auto it = mp.begin(); it != mp.end(); ) {
       if(std::filesystem::exists(parameters_.dir_1.string() + it->second))  {
           ++it;
       } else {
           LOGW<<"WARNING: The file "<<parameters_.dir_1.string() + it->second
           <<" does not exist, so no further analysis will be done on "<<it->second;
            ++bad;
           it = mp.erase( it );
       }
    }
    LOGI<<"Finished phase I! "<<mp.size()<<" files exist, "<<bad<<" don't exist";
    LOGI<<"Beginning phase II: Making sure that the old motifs are a subset of the new ones...";

    for(auto pr : mp) {
        std::cout<<pr.first<<std::endl;
        auto orig_lib = resources::MotifSqliteLibrary(pr.first);
        orig_lib.load_all();
        auto new_lib = resources::MotifSqliteLibrary(1,parameters_.dir_1.string() + pr.second);
        new_lib.load_all();
        auto match(0), miss(0);
        for(const auto& orig_motif : orig_lib)  {
            if(new_lib.contains(orig_motif->to_str()) ) {
                auto new_motif = new_lib.get("","",orig_motif->end_name(0));
                if((*orig_motif->basepairs().begin())->to_str() == (*new_motif->basepairs().begin())->to_str() &&
                        (*orig_motif->basepairs().rbegin())->to_str() == (*new_motif->basepairs().rbegin())->to_str() )  {
                    ++match;
                } else {
                    ++miss;
                }
            } else {
                ++miss;
            }
                    //orig_motif->to_str()//,orig_motif->name(),orig_motif->end_name(0),orig_motif->end_ids()[0],, std::to_string(motif_ct++)
                    //);
            //std::cout<<orig_motif->to_str()<<std::endl; //return;
        }
        std::cout<<"match " <<match<<"\tmiss "<<miss<<std::endl;
    }
}



int main(int argc, char** argv) {
    auto app = CompareDBFiles();
    app.setup_options();
    CLI11_PARSE(app.app_,argc,argv);
    app.run();
}