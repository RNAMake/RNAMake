#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_set>
#include <fstream>

#include <util/x3dna.h>
#include <base/types.h>

#define DEBUG

std::unordered_set<String>
get_cache(String const& path) {
    
    auto cache = std::unordered_set<String>();
    auto infile = std::ifstream(path);
    
    auto name = String{};
    while(std::getline(infile,name,'|')) {
        cache.insert(name);
    }
   
    infile.close();

    return cache;
}

void
save_cache(String const& path, std::unordered_set<String> const & cache) {
    
    auto outfile = std::ofstream(path);

    for(const auto& name : cache) {
        outfile<<name<<'|';
    }
    outfile.close();
}

Strings
get_pdb_files(String const& dir) {
    auto pdb_paths = Strings{};
    for (const auto & entry : std::__fs::filesystem::directory_iterator(dir)) {
        if(entry.path().string().find("pdb") != std::string::npos) {
            pdb_paths.push_back(std::move(entry.path().string()));
        } 
    }
    
    return pdb_paths;
}


void
compare_outputs(Strings const& pdbs) {
    auto i = 0; 
    for(const auto& pdb : pdbs) { 
        std::cout<<"Starting "<<++i<<" of "<<pdbs.size();
        try { 
            auto generator = util::X3dna{}; 
            const auto code = pdb.substr(pdb.find_last_of('/')+1);
            
            auto bps = generator.get_basepairs(pdb);
            auto bps_json = generator.get_basepairs_json(pdb);
             
            auto orig_lines = Strings{}; 
            auto json_lines = Strings{};
            
            for(const auto& bp : bps) orig_lines.push_back(bp.to_string());
            for(const auto& bp : bps_json) json_lines.push_back(bp.to_string());

            auto contents = String{base::join_by_delimiter(orig_lines,"_") + "," + base::join_by_delimiter(json_lines,"_") + "\n"}; 
            
            auto outfile = std::ofstream("all_results.csv",std::ios_base::app); 
            outfile<<contents;
            outfile.close();

        } catch (std::runtime_error& r) {
            std::cout<<"ERROR: "<<r.what()<<std::endl;
            continue; 
        }
            std::cout<<" finished\n";
    }
}

Strings
parse_cl(int argc, char** argv) {
    auto paths = Strings{};
    for(auto ii = 2; ii<argc; ++ii) {
        paths.push_back(argv[ii]);
    } 
    
    return paths; 
}

int main(int argc, char** argv) {

    auto pdbs = get_pdb_files("../../pdb/");
    compare_outputs(pdbs); 
}
