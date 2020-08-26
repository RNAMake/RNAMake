#include<base/sys_interface.h>



namespace base {

String 
execute_command( const char* cmd ) {
    
    std::array<char, 128> buffer;
    String result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    
    return result;
    
}

String
execute_command( String const & cmd) {
    return execute_command(cmd.c_str());
}

nlohmann::json
execute_command_json( const char* cmd ) {
    return nlohmann::json::parse(execute_command(cmd));
}

nlohmann::json
execute_command_json( String const& cmd )  {
    return nlohmann::json::parse(execute_command(cmd));
}

void
json_cleanup() {
    const auto json_temp_files = Strings{
                            "dssr-Aminors.pdb",
                            "dssr-bulges.pdb",
                            "dssr-iloops.pdb",
                            "dssr-2ndstrs.bpseq",
                            "dssr-2ndstrs.ct",
                            "dssr-2ndstrs.dbn",
                            "dssr-atom2bases.pdb",
                            "dssr-hairpins.pdb",
                            "dssr-helices.pdb",
                            "dssr-junctions.pdb",
                            "dssr-multiplets.pdb",
                            "dssr-pairs.pdb",
                            "dssr-splays.pdb",
                            "dssr-stacks.pdb",
                            "dssr-stems.pdb",
                            "dssr-torsions.txt"};
    
    for(const auto& temp_file : json_temp_files) {
        std::remove(temp_file.c_str());
    }
 
}

}
