//
//  x3dna.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include <set>

//RNAMake Headers
#include <base/sys_interface.h>
#include <base/file_io.h>
#include <base/settings.h>
#include <util/x3dna.h>
#include <util/x3dna_src.h>
#include <util/find_pair.h>

namespace util { 
X3dna::X3dna() :
        rebuild_files_(true),
        generated_dssr_(false),
        generated_ref_frames_(false),
        no_ref_frames_(false) {

    auto os_name = base::get_os_name();
}

void
X3dna::_delete_file(
        String const &file_name) const {
    try { std::remove(file_name.c_str()); }
    catch (String const &e) {}
}

math::Point
X3dna::_convert_string_to_point(
        String const &str) const {
    auto doubles = std::vector<double>();
    auto spl = base::split_str_by_delimiter(str, " ");
    for (auto const &s : spl) {
        if (s.length() > 1) {
            doubles.push_back(std::stod(s));
        }
        if (doubles.size() == 3) { break; }
    }
    return math::Point(doubles[0], doubles[1], doubles[2]);
}

void
X3dna::_parse_ref_frame_file(
        String const &pdb_path,
        X3Basepairs &basepairs) const {
    auto finder = PairFinder(pdb_path);
    finder.find_pair(basepairs);
    for (auto const & b : basepairs) {
        std::cout << b.to_string() << std::endl;
        std::cout  << "Res1: " << b.res1.num << std::endl;
        std::cout << "Res2: " << b.res2.num << std::endl;
    }
}

Strings
X3dna::_split_over_white_space(
        String const &str) const {
    Strings spl = base::split_str_by_delimiter(str, " ");
    Strings non_white_space;
    String temp;
    for (auto &s : spl) {
        temp = base::trim(s);
        if (temp.size() == 0) { continue; }
        non_white_space.push_back(temp);
    }
    return non_white_space;
}


X3dna::X3Basepairs
X3dna::get_basepairs_json(
        String const &pdb_path) const {
    auto dssr_json = base::execute_command_json(bin_path_ + "/x3dna-dssr -i=" + pdb_path + " --json --more 2>.error"); 
        
    //deleting the temp files that we don't want. 
    util::json_cleanup();

    auto nt_it = dssr_json.find("nts"); 
    if (nt_it == dssr_json.end() || nt_it->is_null() || nt_it->empty()) {
        return X3Basepairs{}; 
    }
    // loop through the nucleotides and store them in a map 
    auto res_map = std::map<String,X3Residue>{}; 
    for(auto& nt : *dssr_json.find("nts")) {
        auto residue_it = util::get_string(nt,"nt_id");
        auto ii = util::get_int(nt,"nt_resnum");
        auto chain = util::get_char(nt,"chain_name");
        res_map[residue_it] = X3Residue(ii,chain,' ');
    
    }
    // now the basepairs can be constructed
    auto basepairs = X3Basepairs();
    auto pairings = dssr_json.find("pairs"); 
    if(pairings == dssr_json.end() || pairings->empty()) {
        return basepairs;
    }
    // loop through the pairings 
    for(auto pp : *pairings) {
        auto name = util::get_string(pp,"name");
        auto ii =  util::get_string(pp,"Saenger");
        
        auto nt_1 = util::get_string(pp,"nt1");
        auto nt_2 = util::get_string(pp,"nt2");
        auto type = util::get_string(pp,"DSSR");
        // get the reference frame information 
        auto frame = math::Matrix{}; 
        auto origin = math::Point{};
        auto frame_it = pp.find("frame"); 
        
        if(frame_it != pp.end() && !frame_it->is_null()) {
            origin = get_point(*frame_it,"origin");
            frame = get_matrix(*frame_it);
                            
        }
        // define the basepair type and residues        
        auto bp_type = get_x3dna_by_type(type); 
        auto res1 = res_map.at(nt_1); // TODO maybe add an unkown thing?
        auto res2 = res_map.at(nt_2);
        // TODO add check that they are correct 
        // add the new basepair i 
        auto new_bp = X3Basepair{res1,res2,origin,frame,bp_type} ;
        //std::cout<<nt_1<<"\t"<<nt_2<<"\t"<<type<<std::endl; 
        if(new_bp.valid()) {
            basepairs.push_back(std::move(new_bp));
        }
    }
    

    const auto num_basepairs = util::get_int(dssr_json,"num_pairs");

    if ( num_basepairs != basepairs.size() ) {
        throw X3dnaException("Basepair count mismatch. Expected " + std::to_string(num_basepairs) \
                + " but got " + std::to_string(basepairs.size()));
    }
    
    return basepairs; 
}


X3dna::X3Basepairs
X3dna::get_basepairs(
        String const &pdb_path) const {

    auto basepairs = X3Basepairs();

    _parse_ref_frame_file(pdb_path, basepairs);

    return basepairs;

}


X3dnaBPType
get_x3dna_by_type(String const &name) {
    if      (name == "cm-" ) { return X3dnaBPType::cmU; }
    else if (name == "cM-M") { return X3dnaBPType::cMUM; }
    else if (name == "tW+W") { return X3dnaBPType::tWPW; }
    else if (name == "c.+M") { return X3dnaBPType::cDPM; }
    else if (name == ".W+W") { return X3dnaBPType::DWPW; }
    else if (name == "tW-M") { return X3dnaBPType::tWUM; }
    else if (name == "tm-M") { return X3dnaBPType::tmUM; }
    else if (name == "cW+M") { return X3dnaBPType::cWPM; }
    else if (name == ".W-W") { return X3dnaBPType::DWUW; }
    else if (name == "cM+.") { return X3dnaBPType::cMPD; }
    else if (name == "c.-m") { return X3dnaBPType::cDUm; }
    else if (name == "cM+W") { return X3dnaBPType::cMPW; }
    else if (name == "tM+m") { return X3dnaBPType::tMPm; }
    else if (name == "tM-W") { return X3dnaBPType::tMUW; } 
    else if (name == "cm-m") { return X3dnaBPType::cmUm; }
    else if (name == "cM-W") { return X3dnaBPType::cMUW; }
    else if (name == "cW-W") { return X3dnaBPType::cWUW; }
    else if (name == "c.-M") { return X3dnaBPType::cDUM; }
    else if (name == "cm+M") { return X3dnaBPType::cmPM; }
    else if (name == "cm-M") { return X3dnaBPType::cmUM; }
    else if (name == "....") { return X3dnaBPType::DDDD; }
    else if (name == "cm-W") { return X3dnaBPType::cmUW; }
    else if (name == "tM-m") { return X3dnaBPType::tMUm; }
    else if (name == "c.-W") { return X3dnaBPType::cDUW; }
    else if (name == "cM+m") { return X3dnaBPType::cMPm; }
    else if (name == "cM-m") { return X3dnaBPType::cMUm; }
    else if (name == "c...") { return X3dnaBPType::cDDD; }
    else if (name == "tW+m") { return X3dnaBPType::tWPm; }
    else if (name == "c.+m") { return X3dnaBPType::cDPm; }
    else if (name == "tm+m") { return X3dnaBPType::tmPm; }
    else if (name == "tW+.") { return X3dnaBPType::tWPD; }
    else if (name == "tm+W") { return X3dnaBPType::tmPW; }
    else if (name == "t...") { return X3dnaBPType::tDDD; }
    else if (name == "cW-.") { return X3dnaBPType::cWUD; }
    else if (name == "cW-M") { return X3dnaBPType::cWUM; }
    else if (name == "t.-W") { return X3dnaBPType::tDUW; }
    else if (name == "tM+M") { return X3dnaBPType::tMPM; }
    else if (name == "t.-M") { return X3dnaBPType::tDUM; }
    else if (name == "cM-.") { return X3dnaBPType::cMUD; }
    else if (name == "cW-m") { return X3dnaBPType::cWUm; }
    else if (name == "t.+m") { return X3dnaBPType::tDPm; }
    else if (name == "tM-.") { return X3dnaBPType::tMUD; }
    else if (name == "cm+W") { return X3dnaBPType::cmPW; }
    else if (name == "cM+M") { return X3dnaBPType::cMPM; }
    else if (name == "cm+.") { return X3dnaBPType::cmPD; }
    else if (name == "cm-.") { return X3dnaBPType::cmUD; }
    else if (name == "c.-.") { return X3dnaBPType::cDUD; }
    else if (name == "cW+W") { return X3dnaBPType::cWPW; }
    else if (name == "t.-.") { return X3dnaBPType::tDUD; }
    else if (name == "t.+W") { return X3dnaBPType::tDPW; }
    else if (name == "tm-m") { return X3dnaBPType::tmUm; }
    else if (name == "cW+.") { return X3dnaBPType::cWPD; }
    else if (name == "tm+.") { return X3dnaBPType::tmPD; }
    else if (name == "t.+.") { return X3dnaBPType::tDPD; }
    else if (name == "c.+.") { return X3dnaBPType::cDPD; }
    else if (name == "t.-m") { return X3dnaBPType::tDUm; }
    else if (name == "t.+M") { return X3dnaBPType::tDPM; }
// added by CJ 
    else if (name == "tW-.") { return X3dnaBPType::tWUD; }
    else if (name == "tm-W") { return X3dnaBPType::tmUW; }
    else if (name == "tM-M") { return X3dnaBPType::tMUM; }
    else if (name == "tM+.") { return X3dnaBPType::tMPD; }
    else if (name == "c.+W") { return X3dnaBPType::cDPW; }
    else if (name == "tm+M") { return X3dnaBPType::tmPM; }
    else if (name == "tW-m") { return X3dnaBPType::tWUm; }
    else if (name == "cW+m") { return X3dnaBPType::cWPm; }
    else if (name == "tm-.") { return X3dnaBPType::tmUD; }
    else if (name == "tW+M") { return X3dnaBPType::tWPM; }
    else if (name == ".W+m") { return X3dnaBPType::DWPm; }
    else if (name == "tM+W") { return X3dnaBPType::tMPW; }
    else if (name == "..+m") { return X3dnaBPType::DDPm; }
    else if (name == "tW-W") { return X3dnaBPType::tWUW; }
    else if (name == "cm+m") { return X3dnaBPType::cmPm; }
    else if (name == ".W-m") { return X3dnaBPType::DWUm; }
    else if (name == ".M+m") { return X3dnaBPType::DMPm; }
    else if (name == ".W+M") { return X3dnaBPType::DWPM; }
    else if (name == ".M+M") { return X3dnaBPType::DMPM; }
    else if (name == ".m+W") { return X3dnaBPType::DmPW; }
    else if (name == ".W-M") { return X3dnaBPType::DWUM; }
    else if (name == ".m+m") { return X3dnaBPType::DmPm; }
    else if (name == "..-M") { return X3dnaBPType::DDUM; }
    else if (name == ".M-m") { return X3dnaBPType::DMUm; }
    else if (name == "..-m") { return X3dnaBPType::DDUm; }
    else if (name == ".M+W") { return X3dnaBPType::DMPW; }
    else if (name == ".M+.") { return X3dnaBPType::DMPD; }
    else if (name == ".M-M") { return X3dnaBPType::DMUM; }
    else if (name == ".m-m") { return X3dnaBPType::DmUm; }
    else if (name == ".M-W") { return X3dnaBPType::DMUW; }
    else if (name == ".W-.") { return X3dnaBPType::DWUD; }
    else { throw X3dnaException("cannot get x3dna type with: " + name); }
}

String
get_str_from_x3dna_type(
        X3dnaBPType type) {
    if      (type == X3dnaBPType::cmU) { return "cm-";   }
    else if (type == X3dnaBPType::cMUM) { return "cM-M"; }
    else if (type == X3dnaBPType::tWPW) { return "tW+W"; }
    else if (type == X3dnaBPType::cDPM) { return "c.+M"; }
    else if (type == X3dnaBPType::DWPW) { return ".W+W"; }
    else if (type == X3dnaBPType::tWUM) { return "tW-M"; }
    else if (type == X3dnaBPType::tmUM) { return "tm-M"; }
    else if (type == X3dnaBPType::cWPM) { return "cW+M"; }
    else if (type == X3dnaBPType::DWUW) { return ".W-W"; }
    else if (type == X3dnaBPType::cMPD) { return "cM+."; }
    else if (type == X3dnaBPType::cDUm) { return "c.-m"; }
    else if (type == X3dnaBPType::cMPW) { return "cM+W"; }
    else if (type == X3dnaBPType::tMPm) { return "tM+m"; }
    else if (type == X3dnaBPType::tMUW) { return "tM-W"; }
    else if (type == X3dnaBPType::cmUm) { return "cm-m"; }
    else if (type == X3dnaBPType::cMUW) { return "cM-W"; }
    else if (type == X3dnaBPType::cWUW) { return "cW-W"; }
    else if (type == X3dnaBPType::cDUM) { return "c.-M"; }
    else if (type == X3dnaBPType::cmPM) { return "cm+M"; }
    else if (type == X3dnaBPType::cmUM) { return "cm-M"; }
    else if (type == X3dnaBPType::DDDD) { return "...."; }
    else if (type == X3dnaBPType::cmUW) { return "cm-W"; }
    else if (type == X3dnaBPType::tMUm) { return "tM-m"; }
    else if (type == X3dnaBPType::cDUW) { return "c.-W"; }
    else if (type == X3dnaBPType::cMPm) { return "cM+m"; }
    else if (type == X3dnaBPType::cMUm) { return "cM-m"; }
    else if (type == X3dnaBPType::cDDD) { return "c..."; }
    else if (type == X3dnaBPType::tWPm) { return "tW+m"; }
    else if (type == X3dnaBPType::cDPm) { return "c.+m"; }
    else if (type == X3dnaBPType::tmPm) { return "tm+m"; }
    else if (type == X3dnaBPType::tWPD) { return "tW+."; }
    else if (type == X3dnaBPType::tmPW) { return "tm+W"; }
    else if (type == X3dnaBPType::tDDD) { return "t..."; }
    else if (type == X3dnaBPType::cWUD) { return "cW-."; }
    else if (type == X3dnaBPType::cWUM) { return "cW-M"; }
    else if (type == X3dnaBPType::tDUW) { return "t.-W"; }
    else if (type == X3dnaBPType::tMPM) { return "tM+M"; }
    else if (type == X3dnaBPType::tDUM) { return "t.-M"; }
    else if (type == X3dnaBPType::cMUD) { return "cM-."; }
    else if (type == X3dnaBPType::cWUm) { return "cW-m"; }
    else if (type == X3dnaBPType::tDPm) { return "t.+m"; }
    else if (type == X3dnaBPType::tMUD) { return "tM-."; }
    else if (type == X3dnaBPType::cmPW) { return "cm+W"; }
    else if (type == X3dnaBPType::cMPM) { return "cM+M"; }
    else if (type == X3dnaBPType::cmPD) { return "cm+."; }
    else if (type == X3dnaBPType::cmUD) { return "cm-."; }
    else if (type == X3dnaBPType::cDUD) { return "c.-."; }
    else if (type == X3dnaBPType::cWPW) { return "cW+W"; }
    else if (type == X3dnaBPType::tDUD) { return "t.-."; }
    else if (type == X3dnaBPType::tDPW) { return "t.+W"; }
    else if (type == X3dnaBPType::tmUm) { return "tm-m"; }
    else if (type == X3dnaBPType::cWPD) { return "cW+."; }
    else if (type == X3dnaBPType::tmPD) { return "tm+."; }
    else if (type == X3dnaBPType::tDPD) { return "t.+."; }
    else if (type == X3dnaBPType::cDPD) { return "c.+."; }
    else if (type == X3dnaBPType::tDUm) { return "t.-m"; }
    else if (type == X3dnaBPType::tDPM) { return "t.+M"; }
// added by CJ
    else if (type == X3dnaBPType::tWUD) { return "tW-."; }
    else if (type == X3dnaBPType::tmUW) { return "tm-W"; }
    else if (type == X3dnaBPType::tMUM) { return "tM-M"; }
    else if (type == X3dnaBPType::tMPD) { return "tM+."; }
    else if (type == X3dnaBPType::cDPW) { return "c.+W"; }
    else if (type == X3dnaBPType::tmPM) { return "tm+M"; }
    else if (type == X3dnaBPType::tWUm) { return "tW-m"; }
    else if (type == X3dnaBPType::cWPm) { return "cW+m"; }
    else if (type == X3dnaBPType::tmUD) { return "tm-."; }
    else if (type == X3dnaBPType::tWPM) { return "tW+M"; }
    else if (type == X3dnaBPType::DWPm) { return ".W+m"; }
    else if (type == X3dnaBPType::tMPW) { return "tM+W"; }
    else if (type == X3dnaBPType::DDPm) { return "..+m"; }
    else if (type == X3dnaBPType::tWUW) { return "tW-W"; }
    else if (type == X3dnaBPType::cmPm) { return "cm+m"; }
    else if (type == X3dnaBPType::DWUm) { return ".W-m"; }
    else if (type == X3dnaBPType::DMPm) { return ".M+m"; }
    else if (type == X3dnaBPType::DWPM) { return ".W+M"; }
    else if (type == X3dnaBPType::DMPM) { return ".M+M"; }
    else if (type == X3dnaBPType::DmPW) { return ".m+W"; }
    else if (type == X3dnaBPType::DWUM) { return ".W-M"; }
    else if (type == X3dnaBPType::DmPm) { return ".m+m"; }
    else if (type == X3dnaBPType::DDUM) { return "..-M"; }
    else if (type == X3dnaBPType::DMUm) { return ".M-m"; }
    else if (type == X3dnaBPType::DDUm) { return "..-m"; }
    else if (type == X3dnaBPType::DMPW) { return ".M+W"; }
    else if (type == X3dnaBPType::DMPD) { return ".M+."; }
    else if (type == X3dnaBPType::DMUM) { return ".M-M"; }
    else if (type == X3dnaBPType::DmUm) { return ".m-m"; }
    else if (type == X3dnaBPType::DMUW) { return ".M-W"; }
    else if (type == X3dnaBPType::DWUD) { return ".W-."; }

    else { throw X3dnaException("unknown x3dna bp type");}

}

/*
X3dna::X3Motifs
X3dna::get_motifs(
    String const & pdb_path) {
    
    auto dssr_file_path = _get_dssr_file_path(pdb_path);
    auto sections = _divide_dssr_file_into_sections(dssr_file_path);
    X3Motifs all_motifs;
    
    StringStringMap types;
    types["hairpin"] = "HAIRPIN";
    types["bulges"]   = "TWOWAY";
    types["internal"] = "TWOWAY";
    types["junction"]= "NWAY";
    types["non-loop"]= "SSTRAND";
    
    for(auto const & kv : types) {
        if(sections.find(kv.first) != sections.end()) {
            auto motifs = _parse_dssr_section(sections[kv.first], kv.second);
            for(auto const & m : motifs) { all_motifs.push_back(m); }
        }
    }
    if(sections.find("stems") != sections.end()) {
        auto motifs = _parse_dssr_helix_section(sections["stems"]);
        for(auto const & m : motifs) { all_motifs.push_back(m); }
    }
    
    return all_motifs;
    
}
X3dna::X3Motifs
X3dna::_parse_dssr_section(
    Strings const & section,
    String const & mtype) {
    
    X3Motifs motifs;
    X3Residues seen_res;
    int count = 0;
    for(auto const & l : section) {
        auto spl = _split_over_white_space(l);
        if(spl.size() == 0) { continue; }
        try {
            if(spl[0].length() < 3 || spl[0].substr(0,3) != "nts") { continue; }
        } catch(...) { continue; }
        if(spl.size() < 3) { continue; }
        
        auto res_strs = split_str_by_delimiter(spl[2], ",");
        X3Residues res;
        for(auto const & res_str : res_strs) {
            auto res_obj = _parse_dssr_res_str(res_str);
            res.push_back(res_obj);
        }
        count = 0;
        for(auto const & r : res) {
            for(auto const & r2 : seen_res) {
                if(r == r2) { count += 1; break; }
            }
        }
        if(count == res.size()) {
            continue;
        }
        for(auto const & r : res) { seen_res.push_back(r); }
        motifs.push_back(X3Motif{res, mtype});
        
    }
    
    return motifs;
    
}
X3dna::X3Motifs
X3dna::_parse_dssr_helix_section(
    Strings const & section) {
    
    X3Motifs motifs;
    X3Residues res;
    int i = 0;
    for(auto const & l : section) {
        auto spl = _split_over_white_space(l);
        if(spl.size() == 0) { continue; }
        try {
            i = std::stoi(spl[0]);
        } catch(...) { continue; }
        if(i == 1 && res.size() > 0) {
            motifs.push_back(X3Motif{res, "HELIX"});
            res = X3Residues();
        }
        res.push_back(_parse_dssr_res_str(spl[1]));
        res.push_back(_parse_dssr_res_str(spl[2]));
    }
    
    if(res.size() > 0) {
        motifs.push_back(X3Motif{res, "HELIX"});
    }
    
    return motifs;
}
*/ //struct X3Basepair {
    //    X3Residue res1, res2;
    //    math::Point d;
    //    math::Matrix r;
    //    X3dnaBPType bp_type;
    //};

String
X3dna::X3Basepair::to_string() const {
    auto ss = std::stringstream();
    ss << util::get_str_from_x3dna_type(bp_type) << "|" <<res1.num << "|" << res2.num << "|" << d << "|" << math::matrix_to_str(r);
    return ss.str();
}

void
json_cleanup() {

    auto cleanup_cmd = String{base::x3dna_path()} +  String{"/bin/x3dna-dssr --clean 2>"};
    std::system(cleanup_cmd.c_str());

}

String
compare_bps(X3dna::X3Basepairs& lhs, X3dna::X3Basepairs& rhs) {
    // first, we want to build the maps 
    auto left_map = std::map<triplet<int,int,int>,X3dna::X3Basepair>(); 
    auto right_map = std::map<triplet<int,int,int>,X3dna::X3Basepair>(); 
    
    for(auto& bp :lhs) {
        auto key = triplet<int,int,int>(); 
        key.first = round(100.*bp.d.get_x());
        key.second = round(100.*bp.d.get_y());
        key.third = round(100.*bp.d.get_z());
        left_map[key] = bp;
    }
    for(auto& bp : rhs) {
        auto key = triplet<int,int,int>(); 
        key.first = round(100.*bp.d.get_x());
        key.second = round(100.*bp.d.get_y());
        key.third = round(100.*bp.d.get_z());
        right_map[key] = bp;
    }
    // quick check that there are no basepair repeats 
    if (left_map.size() != lhs.size()) {
        throw std::runtime_error("error, redudant residue numbers in lhs bps");
    }
    
    if (right_map.size() != rhs.size()) {
        throw std::runtime_error("error, redudant residue numbers in rhs bps");
    }
    auto matches(0); 
    auto lhs_it = left_map.begin();
    const auto lhs_end = left_map.end();
    while(lhs_it != lhs_end) {
        // check if right_map has the iterator
        auto rhs_it = right_map.find(lhs_it->first);
        // if so, then check if the two are equivalent. delete if so 
        if(rhs_it != right_map.end()) {
            // need to check that the reference frame, origin and bp_type are identical 
            const auto same_origin_flag = roughly_equal(rhs_it->second.d,lhs_it->second.d,0.1);
            const auto same_frame_flag = roughly_equal(rhs_it->second.r,lhs_it->second.r,0.1) || \
                                         roughly_equal(rhs_it->second.r.get_flip_orientation(),lhs_it->second.r,0.1);
            const auto same_bp_type = rhs_it->second.bp_type == lhs_it->second.bp_type;
            if(same_origin_flag && same_frame_flag && same_bp_type ) {
                lhs_it = left_map.erase(lhs_it);
                right_map.erase(rhs_it);
                ++matches;
            } else {
                // otherwise iterate
                ++lhs_it; 
            }
        } else {
                ++lhs_it;
        }

    }
    auto summary = std::to_string(matches) + "," + std::to_string(left_map.size()) + "," + std::to_string(right_map.size()) + ",";
    // check the sizes of the maps... anything leftotver means they are not unique to that method
    if(!left_map.empty()) {
        auto diffs = Strings{}; 
        for(auto& kv : left_map) { 
            diffs.push_back(kv.second.to_string()); 
        }
        summary += base::join_by_delimiter(diffs,"_");
    }
    summary += ",";
    if (!right_map.empty()) {
        auto diffs = Strings{}; 
        for(auto& kv : right_map) { 
            diffs.push_back(kv.second.to_string()); 
        }
        summary += base::join_by_delimiter(diffs,"_");

    }
    summary += "\n";

    return summary;
}
}
