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
            return std::numeric_limits<int>::min();
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

Reals
util::get_reals(const nlohmann::json& json, const String& key) {
        
    auto value = json.find(key);
    if(value != json.end() && !value->is_null()) {
        auto values = Reals{}; 
        for(auto vv : *value) {
            values.push_back(vv);
        }
        
        return values;
    } else {
        return Reals{};
    }
}

Ints
util::get_ints(const nlohmann::json& json, const String& key) {
        
    auto value = json.find(key);
    if(value != json.end() && !value->is_null() ) {
        auto values = Ints{}; 
        for(auto& vv : *value) {
            values.push_back(vv);
        }

        return values;
    } else {
        return Ints{};
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
    auto pairs = util::DssrPairs{}; 

    for(auto& nt_pair : *data.find("pairs")) {
        auto new_pair = util::DssrPair{}; 
        
        new_pair.C1C1_dist = get_double(nt_pair,"C1C1_dist");
        new_pair.C6C8_dist = get_double(nt_pair,"C6C8_dist");     
        new_pair.CNNC_torsion = get_double(nt_pair,"CNNC_torsion"); 
        new_pair.DSSR = get_string(nt_pair,"DSSR");
        new_pair.LW = get_string(nt_pair,"LW");
        new_pair.N1N9_dist = get_double(nt_pair,"N1N9_dist");     
        new_pair.Saenger = get_string(nt_pair,"Saenger");
        new_pair.bp = get_string(nt_pair,"bp");
        new_pair.bp_params = get_reals(nt_pair,"bp_params");
        new_pair.chi1 = get_double(nt_pair,"chi1"); 
        new_pair.chi2 = get_double(nt_pair,"chi2"); 
        new_pair.conf1 = get_string(nt_pair,"conf1"); 
        new_pair.conf2 = get_string(nt_pair,"conf2"); 

        auto frame = nt_pair.find("frame");
        
        if(frame != nt_pair.end() && !frame->is_null()) {
            
            new_pair.frame_origin = get_point(*frame,"origin");
            new_pair.frame_quaternion = get_quaternion(*frame,"quaternion");
            new_pair.frame_rmsd = get_double(*frame,"rmsd"); 
            new_pair.ref_frame = get_matrix(*frame);
                            
        }
        
        new_pair.hbonds_desc = get_string(nt_pair,"hbonds_desc"); 
        new_pair.hbonds_num = get_int(nt_pair,"hbonds_num");
        new_pair.index = get_int(nt_pair,"index");
        new_pair.interBase_angle = get_double(nt_pair,"interBase_angle"); 
        new_pair.lambda1 = get_double(nt_pair,"lambda1"); 
        new_pair.lambda2 = get_double(nt_pair,"lambda2"); 
        new_pair.name = get_string(nt_pair,"name"); 
        new_pair.nt1 = get_string(nt_pair,"nt1"); 
        new_pair.nt2 = get_string(nt_pair,"nt2"); 
        new_pair.planarity = get_double(nt_pair,"planarity");
        new_pair.pucker1 = get_string(nt_pair,"pucker1"); 
        new_pair.pucker2 = get_string(nt_pair,"pucker2"); 
        new_pair.simple_Buckle = get_double(nt_pair,"simple_Buckle"); 
        new_pair.simple_Propeller = get_double(nt_pair,"simple_Propeller"); 
        new_pair.simple_Shear = get_double(nt_pair,"simple_Shear"); 
        new_pair.simple_Stretch = get_double(nt_pair,"simple_Stretch"); 

        pairs.push_back(std::move(new_pair)); 
    }

    return pairs;
}

util::DssrHairpins
util::get_hairpins(const nlohmann::json& data) {
    auto hairpins = util::DssrHairpins{};
    
    for(auto& hairpin : *data.find("hairpins")) {
        
        if(hairpin.is_null()) continue;  // need to add something like this for everywhere
        
        auto new_hairpin = util::DssrHairpin{};
        new_hairpin.index = get_int(hairpin,"index");
        new_hairpin.nts_long = get_string(hairpin,"nts_long");
        new_hairpin.nts_short = get_string(hairpin,"nts_short");
        new_hairpin.type = get_string(hairpin,"type"); 
        new_hairpin.bridging_nts = get_ints(hairpin,"bridging_nts");
        new_hairpin.stem_indices = get_ints(hairpin,"stem_indices");
        new_hairpin.summary = get_string(hairpin,"summary"); 
        new_hairpin.num_nts = get_int(hairpin,"num_nts"); 
        new_hairpin.num_stems = get_int(hairpin,"num_stems"); 


        hairpins.push_back(std::move(new_hairpin)); 
    }
    return hairpins;
}

util::DssrHelices
util::get_helices(const nlohmann::json& data) {
    auto helices = DssrHelices{};
    for(auto& helix : *data.find("helices") ) {
        auto new_helix = util::DssrHelix{};
        new_helix.index = get_int(helix,"index");     
        new_helix.num_stems = get_int(helix,"num_stems");     
        new_helix.strand1 = get_string(helix,"strand1");
        new_helix.strand2 = get_string(helix,"strand2");
        new_helix.bp_type = get_string(helix,"bp_type");
        new_helix.helix_form = get_string(helix,"helix_form");
        new_helix.helical_rise = get_double(helix,"helical_rise");
        new_helix.helical_rise_std = get_double(helix,"helical_rise_std");
        new_helix.helical_radius = get_double(helix,"helical_radius");
        new_helix.helical_radius_std = get_double(helix,"helical_radius_std");
        new_helix.helical_axis = get_point(helix,"helical_axis"); 
        new_helix.point1 = get_point(helix,"point2"); 
        new_helix.point2 = get_point(helix,"point2"); 
        new_helix.num_pairs = get_int(helix,"num_pairs");     
        new_helix.pairs = get_pairs(helix);
 

        helices.push_back(std::move(new_helix));
    }
    return helices;
}

util::DssrStems
util::get_stems(const nlohmann::json& data) {
    auto stems = DssrStems{};
    for(auto& stem : *data.find("stems") ) {
        auto new_stem = util::DssrStem{};
        new_stem.index = get_int(stem,"index");     
        new_stem.strand1 = get_string(stem,"strand1");
        new_stem.strand2 = get_string(stem,"strand2");
        new_stem.bp_type = get_string(stem,"bp_type");
        new_stem.helix_form = get_string(stem,"helix_form");
        new_stem.helical_rise = get_double(stem,"helical_rise");
        new_stem.helical_rise_std = get_double(stem,"helical_rise_std");
        new_stem.helical_radius = get_double(stem,"helical_radius");
        new_stem.helical_radius_std = get_double(stem,"helical_radius_std");
        new_stem.helical_axis = get_point(stem,"helical_axis"); 
        new_stem.point1 = get_point(stem,"point2"); 
        new_stem.point2 = get_point(stem,"point2"); 
        new_stem.num_pairs = get_int(stem,"num_pairs");     
        new_stem.pairs = get_pairs(stem);
 

        stems.push_back(std::move(new_stem));
    }
    return stems;


}

util::DssrILoops
util::get_iloops(const nlohmann::json& data ) {
    auto iloops = DssrILoops{};
    for(auto& iloop : *data.find("iloops")) {
        auto new_loop = DssrILoop{};      
        
        new_loop.index = get_int(iloop,"index");
        new_loop.type = get_string(iloop,"type"); 
        new_loop.num_nts = get_int(iloop,"num_nts");  
        new_loop.num_stems = get_int(iloop,"num_stems");  
        new_loop.nts_short = get_string(iloop,"nts_short"); 
        new_loop.nts_long = get_string(iloop,"nts_long"); 
        new_loop.summary = get_string(iloop,"summary"); 
        new_loop.bridging_nts = get_ints(iloop,"bridging_nts"); 
        new_loop.stem_indices = get_ints(iloop,"stem_indices"); 
        
        iloops.push_back(std::move(new_loop)); 
    }
    return iloops;
}



void
util::get_elements(
        String const & pdb_path,
        util::DssrNts & nts,
        util::DssrPairs & pairs,
        util::DssrHairpins & hairpins,
        util::DssrHelices & helices,
        util::DssrStems & stems,
        util::DssrILoops & iloops
        ) {
    auto dssr_json = base::execute_command_json("../../resources/x3dna/osx/bin/x3dna-dssr -i=" + pdb_path + " --json --more 2> /dev/null");
    if(!nts.empty()) {std::cout<<"Warning, nts vector is NOT empty, but data will be compeletely overwritten";}
    for(auto& kv : dssr_json.items()) {std::cout<<kv.key()<<std::endl;} 
    nts = get_nts(dssr_json);
    pairs = get_pairs(dssr_json);
    hairpins = get_hairpins(dssr_json);
    helices = get_helices(dssr_json);
    stems = get_stems(dssr_json);
    iloops = get_iloops(dssr_json);
}
