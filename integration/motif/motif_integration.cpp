#include <fstream>
#include <cstdlib>
#include <iostream>
#include <dirent.h>
#include <unordered_set>

#include <base/types.h>
#include <base/string.h>
#include <math/xyz_vector.h>
#include <motif/motif_factory.h>
#include <util/x3dna.h>

#define CATCH_CONFIG_MAIN
#include "../../unittests/catch.hpp"

#define DEBUG
// helper method that creates a random point 
class CodeCache{
private:
    std::unordered_set<String> codes;
private:
    String path_;
public:
    CodeCache(String const& path) : path_(path) {
        auto infile = std::ifstream(path); 
        auto code = String{};
        
        while(std::getline(infile,code,'|')) {
            codes.insert(code);     
        }

        infile.close();
    }

public:  
    ~CodeCache() {}

public:
    bool
    contains(String const& path) { 
        return codes.find(path) != codes.end();
    }
public:
    void
    add(String const& code) {
        codes.insert(code);
        
        auto afile = std::fstream(path_,std::ios_base::app);
        afile<<code<<'|';
        afile.close();
    }
};


math::Point
random_point(double lower_limit=-100.,double upper_limit=100.) {
    const auto range = upper_limit - lower_limit; 
    const auto x_val = range*double(rand())/double(RAND_MAX) + lower_limit;
    const auto y_val = range*double(rand())/double(RAND_MAX) + lower_limit;
    const auto z_val = range*double(rand())/double(RAND_MAX) + lower_limit;

    return math::Point{x_val,y_val,z_val};
}

Strings
get_file_names(String const& base_dir, String const& ext ) {
    auto paths = Strings{};

    DIR *target_data_dir;
    struct dirent *current_dir;
    if ((target_data_dir = opendir (base_dir.c_str())) != nullptr) {
        /* print all the files and directories within directory */
        while ((current_dir = readdir (target_data_dir)) != nullptr) {
            std::string directory_entry_name(current_dir->d_name);

            if (directory_entry_name.find(ext) == directory_entry_name.size() - ext.size()) {
                    paths.emplace_back(base_dir + '/' + current_dir->d_name);
            }
        }
        closedir (target_data_dir);
    }
    return paths;

}

TEST_CASE("Integration test for motif code") {
#ifdef JSON_BASEPAIRS
    std::cout<<"Using JSON get_basepairs()"<<std::endl;
#endif
    std::srand(1003);
    auto mf = motif::MotifFactory{}; 
    auto cache = CodeCache(".cache");
    auto failed = CodeCache(".failed");
    auto paths = get_file_names("../../pdb/",".pdb"); 
    const auto total = paths.size();
    auto i(1);
    SECTION("random alignment for all the pdbs first") {

        for(const auto& pdb : paths ) {
            if(cache.contains(pdb) || failed.contains(pdb)) {
                ++i; 
                continue;
            }
            std::cout<<i++<<" of "<<total<<std::endl; 
            try{ 

            std::cout<<pdb<<std::endl;
            auto motif = mf.motif_from_file("../../pdb/4QJH.pdb");
            //auto motif = mf.motif_from_file(pdb);
            break; 
            if(motif->ends().empty()) {
                continue;
            }

            const auto pt = random_point();
            motif->move(pt);
            // check that they all can be loaded in 
            // not guaranteed to have an end... keep that in mind
            auto motif2 = mf.motif_from_file(pdb);
            auto m_aligned = get_aligned_motif(motif2->ends()[0], motif->ends()[0], motif);
            
            auto dist = motif2->ends()[0]->d().distance(m_aligned->ends()[0]->d());
            auto r_dist = motif2->ends()[0]->r().difference(m_aligned->ends()[0]->r());

            REQUIRE(dist < 1.0);
            REQUIRE(r_dist < 0.001);
            
            cache.add(pdb); 
            
            } catch(std::runtime_error& E) {
                std::cout<<E.what()<<std::endl;
                failed.add(pdb); 
                continue;
            }

        }


    }
    //        m->move(math::Point(10, 10, 10));

}
