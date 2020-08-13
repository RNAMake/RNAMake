#include <util/dssr.h>
#include <base/sys_interface.h>



util::DssrNts
util::get_nts(const nlohmann::json& data) {
    //auto res_map = std::map<String,X3Residue>{}; 
    for(auto& nt : *data.find("nts")) {
        for(auto kv : nt.items()) {
            std::cout<<kv.key()<<"\t"<<kv.value()<<std::endl;
        }

            std::cout<<"----------------------------"<<std::endl;
        //auto residue_it = nt.find("nt_id")->get<String>();
        //auto ii = nt.find("nt_resnum")->get<int>();
        //auto chain = nt.find("chain_name")->get<String>()[0];
        //res_map[residue_it] = X3Residue(ii,chain,' ');
    
    }
return DssrNts{};
}

void
util::get_elements(
        String const & pdb_path,
        util::DssrNts & nts
        ) {
    auto dssr_json = base::execute_command_json("../../resources/x3dna/osx/bin/x3dna-dssr -i=" + pdb_path + " --json --more 2> /dev/null");
    if(!nts.empty()) {std::cout<<"Warning, nts vector is NOT empty, but data will be compeletely overwritten";}

    nts = get_nts(dssr_json);
}
