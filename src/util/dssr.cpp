#include <util/dssr.h>
#include <base/sys_interface.h>

double
util::get_double(const nlohmann::json& json, const String& key) {
    
        auto value = json.find(key);
        if(value != json.end() && !value->is_null()) {
            return value->get<double>(); 
        } else {
            return -1.;
        }
}

int
util::get_int(const nlohmann::json& json, const String& key) {
    
        auto value = json.find(key);
        if(value != json.end() && !value->is_null()) {
            return value->get<int>(); 
        } else {
            return -1;
        }
}

char
util::get_char(const nlohmann::json& json, const String& key) {
    
        auto value = json.find(key);
        if(value != json.end() && !value->is_null()) {
            return value->get<String>()[0]; 
        } else {
            return ' ';
        }
}

String
util::get_string(const nlohmann::json& json, const String& key) {
    
        auto value = json.find(key);
        if(value != json.end() && !value->is_null()) {
            return value->get<String>(); 
        } else {
            return "NA";
        }
}

math::Point
util::get_point(const nlohmann::json& json, const String& key) {
        
        auto value = json.find(key);
        if(value != json.end() &&  
                !value->at(0).is_null() && 
                !value->at(1).is_null() && 
                !value->at(2).is_null()) {

            return math::Point{value->at(0),value->at(1),value->at(2)};
        } else {
            return math::Point{-1.,-1.,-1.};
        }

}

math::Quaternion
util::get_quaternion(const nlohmann::json& json, const String& key ) {

        auto value = json.find(key);
         if(value != json.end() &&
                    !value->at(0).is_null() && 
                        !value->at(1).is_null() && 
                            !value->at(2).is_null() && 
                                !value->at(3).is_null()) {
                return math::Quaternion(value->at(0),value->at(1),value->at(2),value->at(3));
         } else {
                return math::Quaternion(-1.,-1.,-1.,-1.);
         }


}

math::Matrix
util::get_matrix(const nlohmann::json& json) {

        auto x_axis = json.find("x_axis");
        auto y_axis = json.find("y_axis");
        auto z_axis = json.find("z_axis");
        
        if(*x_axis != nullptr && *y_axis != nullptr && *z_axis != nullptr && 
                        !x_axis->at(0).is_null() && 
                        !x_axis->at(1).is_null() && 
                        !x_axis->at(2).is_null() && 
                        !y_axis->at(0).is_null() && 
                        !y_axis->at(1).is_null() && 
                        !y_axis->at(2).is_null() && 
                        !z_axis->at(0).is_null() &&
                        !z_axis->at(1).is_null() &&
                        !z_axis->at(2).is_null()             

                    ) {
            return  math::Matrix(
                                    x_axis->at(0),x_axis->at(1),x_axis->at(2),
                                    y_axis->at(0),y_axis->at(1),y_axis->at(2),
                                    z_axis->at(0),z_axis->at(1),z_axis->at(2)
                                    );

        } else {
            return  math::Matrix(-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.);
        }

}


util::DssrNts
util::get_nts(const nlohmann::json& data) {
    
    auto nts = util::DssrNts{}; 
    
    for(auto& nt : *data.find("nts")) {
        
        auto new_nt = util::DssrNt{}; 
        
        new_nt.C5prime_xyz = get_point(nt,"C5prime_xyz");
        new_nt.Dp = get_double(nt,"Dp");       
    
        auto P_xyz = nt.find("P_xyz");
        new_nt.P_xyz = get_point(nt,"P_xyz"); 
        
        new_nt.alpha = get_double(nt,"alpha"); 
        new_nt.amplitude = get_double(nt,"amplitude"); 
        new_nt.bb_type = get_string(nt,"bb_type"); 
        new_nt.beta = get_double(nt,"beta");
        new_nt.bin = get_string(nt,"bin");  
        new_nt.chain_name = get_char(nt,"chain_name");
        new_nt.chi = get_double(nt,"chi");
        new_nt.cluster = get_char(nt,"cluster");
        new_nt.dbn = get_char(nt,"dbn");
        new_nt.epsilon = get_double(nt,"epsilon");
        new_nt.epsilon_zeta = get_double(nt,"epsilon_zeta");
        new_nt.eta = get_double(nt,"eta");
        new_nt.eta_base = get_double(nt,"eta_base");
        new_nt.eta_prime = get_double(nt,"eta_prime");
        new_nt.form = get_char(nt,"form");
        
        auto frame = nt.find("frame");
        
        if(frame != nt.end() && !frame->is_null()) {
            
            new_nt.frame_origin = get_point(*frame,"origin");
            new_nt.frame_quaternion = get_quaternion(*frame,"quaternion");
            new_nt.frame_rmsd = get_double(*frame,"rmsd"); 
            new_nt.ref_frame = get_matrix(*frame);
                            
        }
        
        new_nt.gamma = get_double(nt,"gamma");        
        new_nt.glyco_bond = get_string(nt,"glyco_bond");
        new_nt.index = get_int(nt,"index");
        new_nt.index_chain = get_int(nt,"index_chain");
        new_nt.nt_code = get_char(nt,"nt_code");
        new_nt.nt_id = get_string(nt,"nt_id");
        new_nt.nt_name = get_char(nt,"nt_name");
        new_nt.nt_resnum = get_int(nt,"nt_resnum");
        new_nt.phase_angle = get_double(nt,"bad_key");  
        new_nt.puckering = get_string(nt,"puckering");
        new_nt.splay_angle = get_double(nt,"splay_angle");
        new_nt.splay_distance = get_double(nt,"splay_distance");
        new_nt.splay_ratio = get_double(nt,"splay_ratio");
        new_nt.ssZp = get_double(nt,"ssZp");
        new_nt.sugar_class = get_string(nt,"sugar_class"); 
        new_nt.suiteness = get_double(nt,"suiteness"); 
        new_nt.summary = get_string(nt,"summary"); 
        new_nt.theta = get_double(nt,"theta");
        new_nt.theta_prime = get_double(nt,"theta_prime");        
        new_nt.v0 = get_double(nt,"v0"); 
        new_nt.v1 = get_double(nt,"v1"); 
        new_nt.v2 = get_double(nt,"v2"); 
        new_nt.v3 = get_double(nt,"v3"); 
        new_nt.v4 = get_double(nt,"v4"); 
        new_nt.zeta = get_double(nt,"zeta"); 

        nts.push_back(std::move(new_nt));    
    }

    return nts;
}

util::DssrPairs
util::get_pairs(const nlohmann::json& data) {

    for(auto& nt_pair : *data.find("pairs")) {
        std::cout<<"----------------------------------------"<<std::endl;
        for(auto& kv : nt_pair.items()) 
            std::cout<<kv.key()<<"\t"<<kv.value()<<std::endl;
    }

    return DssrPairs{};
}

void
util::get_elements(
        String const & pdb_path,
        util::DssrNts & nts,
        util::DssrPairs & pairs
        ) {
    auto dssr_json = base::execute_command_json("../../resources/x3dna/osx/bin/x3dna-dssr -i=" + pdb_path + " --json --more 2> /dev/null");
    if(!nts.empty()) {std::cout<<"Warning, nts vector is NOT empty, but data will be compeletely overwritten";}
    for(auto kv :  dssr_json.items()) {std::cout<<kv.key()<<std::endl;}
    nts = get_nts(dssr_json);
    pairs = get_pairs(dssr_json);
}
