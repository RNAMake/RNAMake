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


String
compare_outputs(Strings const& pdbs) {
    auto contents = String{}; 
    for(const auto& pdb : pdbs) { 
        
        auto generator = util::X3dna{}; 
        auto bps = generator.get_basepairs(pdb);
        auto bps_json = generator.get_basepairs_json(pdb);
        const auto code = pdb.substr(pdb.find_last_of('/')+1);
        contents += (code + "," + util::compare_bps(bps,bps_json));
    }
    return contents;
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
    auto paths = parse_cl(argc,argv);
    auto ct = String{argv[1]};
    
    auto contents = compare_outputs(paths); 
    auto outfile = std::ofstream(ct + ".csv");
    outfile << contents<<std::endl;
    outfile.close();
}
