//
//  x3dna.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

//RNAMake Headers
#include "base/file_io.h"
#include "util/x3dna.h"

X3dna::X3dna() {
    String env = "X3DNA="+x3dna_path();
    s_ = strdup(env.c_str());
    putenv(s_);
    
    bin_path_ = x3dna_path()+"/bin/";
}

void
X3dna::generate_ref_frame(String const & path) {
    
    String fname = filename(path).substr(0,-4);
    fname = fname.substr(0, fname.length()-4);
    if(!file_exists(path)) {
        throw X3dnaException("cannot find pdb for ref_frames.dat\n");
    }
    
    String find_pair_path = bin_path_ + "find_pair ";
    String analyze_path   = bin_path_ + "analyze ";
    String command = find_pair_path + path + " 2> /dev/null stdout | " + analyze_path + "stdin >& /dev/null";
    char * s = strdup(command.c_str());
    int result = std::system(s);

    if(result != 0) {
        throw X3dnaException("could not call find_pair properly, please make sure you have it set up properly\n");
    }
    String files_str = "auxiliary.par,bestpairs.pdb,bp_helical.par,bp_order.dat,bp_step.par,cf_7methods.par,col_chains.scr,col_helices.scr,hel_regions.pdb,hstacking.pdb,poc_haxis.r3d,stacking.pdb";
    Strings files = split_str_by_delimiter(files_str, ",");
    files.push_back(fname+".out");
    for(auto const & f : files) {
        try { std::remove(f.c_str()); }
        catch(String const & e) {}
    }
    delete s;
    
}

void
X3dna::generate_dssr_file(String const & path) {
    
    String fname = filename(path);
    fname = fname.substr(0, fname.length()-4);
    if(!file_exists(path)) {
        throw X3dnaException("cannot find pdb for generate_dssr_file\n");
    }
    
    String dssr_path = bin_path_ + "x3dna-dssr ";
    String command = dssr_path + "-i="+path+" -o="+fname+"_dssr.out --non-pair >& /dev/null";
    char * s = strdup(command.c_str());
    int result = std::system(s);

    if(result != 0) {
        throw X3dnaException("could not call x3dna-dssr properly, please make sure you have it set up properly\n");
    }
    
    String filename_str = "dssr-2ndstrs.ct,dssr-2ndstrs.dbn,dssr-helices.pdb,dssr-pairs.pdb,dssr-stems.pdb,hel_regions.pdb,hstacking.pdb,poc_haxis.r3d,stacking.pdb,dssr-torsions.dat,dssr-Kturns.pdb,dssr-multiplets.pdb,dssr-hairpins.pdb,dssr-Aminors.pdb";
    Strings files = split_str_by_delimiter(filename_str, ",");
    for(auto const & f : files) {
        try { std::remove(f.c_str()); }
        catch(String const & e) {}
    }
    delete s;

}

String
X3dna::_get_ref_frame_path(
    String const & pdb_name,
    bool force_build_files) {
    
    String basedir = base_dir(pdb_name);
    String ref_frames_path = "";
    if(force_build_files) {
        generate_ref_frame(pdb_name);
        ref_frames_path = "ref_frames.dat";
    }
    
    // is pdb
    else if(file_exists(basedir + "/ref_frames.dat")) {
        ref_frames_path = basedir + "/ref_frames.dat";
    }
    // is directory
    else if(file_exists(pdb_name + "/ref_frames.dat")) {
        ref_frames_path = pdb_name + "/ref_frames.dat";
    }
    
    else if(file_exists("ref_frames.dat")) {
        ref_frames_path = "ref_frames.dat";
    }
    else {
        generate_ref_frame(pdb_name);
        ref_frames_path = "ref_frames.dat";
    }
    
    return ref_frames_path;
    
}

String
X3dna::_get_dssr_file_path(
    String const & pdb_path,
    bool force_build_files) {
    String basedir = base_dir(pdb_path);
    String fname = filename(pdb_path);
    fname = fname.substr(0, fname.length()-4);
    String dssr_name = fname + "_dssr.out";
    String dssr_file_path = "";
    
    if(force_build_files) {
        generate_dssr_file(pdb_path);
        dssr_file_path = dssr_name;
    }
    
    else if(file_exists(basedir + "/" + dssr_name)) {
        dssr_file_path = basedir + "/" + dssr_name;
    }
    
    else if(file_exists(pdb_path + "/" + dssr_name )) {
        dssr_file_path = pdb_path + "/" + dssr_name;
    }
    else {
        generate_dssr_file(pdb_path);
        dssr_file_path = dssr_name;
        
    }
    return dssr_file_path;
}

Point
X3dna::_convert_strings_to_point(Strings const & spl) {
    std::vector<double> f;
    for (auto const & s : spl) {
        if(s.length() > 1) {
            f.push_back(std::stof(s));
        }
        if(f.size() == 3) { break; }
    }
    return Point(f);
}

void
X3dna::_parse_ref_frame_file(String const & ref_frames_path) {
    Strings lines = get_lines_from_file(ref_frames_path);
    basepairs_ = X3Basepairs();

    Point d;
    Points rs;
    Matrix r;
    int res_num_1=0, res_num_2=0;
    String res_id_1, res_id_2, res_num_1p, res_num_2p;
    int start_bp = 0;
    
    for(auto const & l : lines) {
        if(l.length() < 3) { continue; }
        if(l.substr(0,3).compare("...") == 0 ) {
            start_bp = 1;
            res_id_1 = l.substr(18,1);
            res_id_2 = l.substr(35,1);
            res_num_1p = ""; res_num_2p = "";
            for(auto const & c : l.substr(20,4)) {
                if(isdigit(c)) { res_num_1p += c; }
            }
            for(auto const & c : l.substr(37,4)) {
                if(isdigit(c)) { res_num_2p += c; }
            }
            res_num_1 = std::stoi(res_num_1p);
            res_num_2 = std::stoi(res_num_2p);
            continue;
        }
        if     (start_bp == 0) { continue; }
        if     (start_bp == 1) {
            Strings spl = split_str_by_delimiter(l, " ");
            d = _convert_strings_to_point(spl);
            start_bp++;
        }
        else if(start_bp > 1)  {
            Strings spl = split_str_by_delimiter(l, " ");
            rs.push_back(_convert_strings_to_point(spl));
            start_bp++;
        }
        
        if(start_bp == 5) {
            start_bp = 0;
            
            r.xx(rs[0][0]); r.xy(rs[0][1]); r.xz(rs[0][2]);
            r.yx(rs[1][0]); r.yy(rs[1][1]); r.yz(rs[1][2]);
            r.zx(rs[2][0]); r.zy(rs[2][1]); r.zz(rs[2][2]);
            auto res1 = X3Residue{res_num_1, res_id_1[0], ' '};
            auto res2 = X3Residue{res_num_2, res_id_2[0], ' '};
            auto bp   = X3Basepair{res1, res2, d, r, X3dnaBPType::cDDD};
            basepairs_.push_back(bp);
            rs.resize(0);
        }
    }
}

std::map<String, Strings>
X3dna::_divide_dssr_file_into_sections(String const & dssr_file_path) {
    Strings lines = get_lines_from_file(dssr_file_path);
    
    std::map<String, Strings> sections;
    Strings section;
    String section_name;
    
    for(auto const & l : lines) {
        if(l.substr(0,4).compare("List") == 0 && section.size() == 0) {
            Strings spl = split_str_by_delimiter(l, " ");
            section_name = spl[3];
        }
        if(l[0] == '*') {
            if(section.size() != 0 && sections.find(section_name) == sections.end()) {
                sections[section_name] = section;
            }
            section.resize(0);
            continue;
        }
        section.push_back(l);
    }
    
    return sections;
    
}

Strings
X3dna::_split_over_white_space(String const & str) {
    Strings spl = split_str_by_delimiter(str, " ");
    Strings non_white_space;
    String temp;
    for(auto & s : spl) {
        //if(!(s.len) {
        temp = trim(s);
        if(temp.size() == 0) { continue; }
        non_white_space.push_back(temp);
        //}
    }
    return non_white_space;
}

X3dna::X3Residue
X3dna::_parse_dssr_res_str(String const & res_str) {
    auto spl = split_str_by_delimiter(res_str, ".");
    auto chain = spl[0][0];
    auto rnum = spl[1].substr(1);
    int num;
    try{
        num = std::stoi(rnum);
    } catch(...) {
        std::cout << res_str << std::endl;
        throw "could not parse " + res_str + " into a residue\n";
    }

    return X3Residue{num, chain, ' '};
    
}

X3dna::X3Basepairs const &
X3dna::get_basepairs(
    String const & pdb_path,
    bool force_build_files) {
    
    auto ref_frames_path = _get_ref_frame_path(pdb_path, force_build_files);
    auto dssr_file_path  = _get_dssr_file_path(pdb_path, force_build_files);
    
    _parse_ref_frame_file(ref_frames_path);
    auto sections = _divide_dssr_file_into_sections(dssr_file_path);
    auto section = sections["base"];
    int found = 0;
    
    for(auto const & l : section) {
        if(l.length() < 3) { continue; }
        auto spl = _split_over_white_space(l);
        auto spl2 = split_str_by_delimiter(spl[1], ".");
        if(spl2.size() == 1) {continue; }
        if(spl.size() < 6) { continue; }
        auto res1 = _parse_dssr_res_str(spl[1]);
        auto res2 = _parse_dssr_res_str(spl[2]);
        auto bp_type = X3dnaBPType::cDDD;
        //TODO look into why this is happening, sometimes will error out if I dont do this check
        if(spl.size() > 7) {
            bp_type = get_x3dna_by_type(spl[7]);
        }
        found = 0;
        for(auto & bp : basepairs_) {
            if(bp.res1 == res1 && bp.res2 == res2) {
                bp.bp_type = bp_type;
                found = 1;
            }
            if(bp.res1 == res2 && bp.res2 == res1) {
                bp.bp_type = bp_type;
                found = 1;
            }
            
            if(found) { break; }
        }
        
        if(!found) {
            basepairs_.push_back(X3Basepair{res1, res2, Point(-1,-1,-1), Matrix(),
                                            X3dnaBPType::cDDD});
        }
        
    }    
    return basepairs_;
    


}

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


X3dna::X3dnaBPType
get_x3dna_by_type(String const & name) {
    if     (name == "cm-")  { return X3dna::X3dnaBPType::cmU; }
    else if(name == "cM-M") { return X3dna::X3dnaBPType::cMUM; }
    else if(name == "tW+W") { return X3dna::X3dnaBPType::tWPW; }
    else if(name == "c.+M") { return X3dna::X3dnaBPType::cDPM; }
    else if(name == ".W+W") { return X3dna::X3dnaBPType::DWPW; }
    else if(name == "tW-M") { return X3dna::X3dnaBPType::tWUM; }
    else if(name == "tm-M") { return X3dna::X3dnaBPType::tmUM; }
    else if(name == "cW+M") { return X3dna::X3dnaBPType::cWPM; }
    else if(name == ".W-W") { return X3dna::X3dnaBPType::DWUW; }
    else if(name == "cM+.") { return X3dna::X3dnaBPType::cMPD; }
    else if(name == "c.-m") { return X3dna::X3dnaBPType::cDUm; }
    else if(name == "cM+W") { return X3dna::X3dnaBPType::cMPW; }
    else if(name == "tM+m") { return X3dna::X3dnaBPType::tMPm; }
    else if(name == "tM-W") { return X3dna::X3dnaBPType::tMUW; }
    else if(name == "cm-m") { return X3dna::X3dnaBPType::cmUm; }
    else if(name == "cM-W") { return X3dna::X3dnaBPType::cMUW; }
    else if(name == "cW-W") { return X3dna::X3dnaBPType::cWUW; }
    else if(name == "c.-M") { return X3dna::X3dnaBPType::cDUM; }
    else if(name == "cm+M") { return X3dna::X3dnaBPType::cmPM; }
    else if(name == "cm-M") { return X3dna::X3dnaBPType::cmUM; }
    else if(name == "....") { return X3dna::X3dnaBPType::DDDD; }
    else if(name == "cm-W") { return X3dna::X3dnaBPType::cmUW; }
    else if(name == "tM-m") { return X3dna::X3dnaBPType::tMUm; }
    else if(name == "c.-W") { return X3dna::X3dnaBPType::cDUW; }
    else if(name == "cM+m") { return X3dna::X3dnaBPType::cMPm; }
    else if(name == "cM-m") { return X3dna::X3dnaBPType::cMUm; }
    else if(name == "c...") { return X3dna::X3dnaBPType::cDDD; }
    else if(name == "tW+m") { return X3dna::X3dnaBPType::tWPm; }
    else if(name == "c.+m") { return X3dna::X3dnaBPType::cDPm; }
    else if(name == "tm+m") { return X3dna::X3dnaBPType::tmPm; }
    else if(name == "tW+.") { return X3dna::X3dnaBPType::tWPD; }
    else if(name == "tm+W") { return X3dna::X3dnaBPType::tmPW; }
    else if(name == "t...") { return X3dna::X3dnaBPType::tDDD; }
    else if(name == "cW-.") { return X3dna::X3dnaBPType::cWUD; }
    else if(name == "cW-M") { return X3dna::X3dnaBPType::cWUM; }
    else if(name == "t.-W") { return X3dna::X3dnaBPType::tDUW; }
    else if(name == "tM+M") { return X3dna::X3dnaBPType::tMPM; }
    else if(name == "t.-M") { return X3dna::X3dnaBPType::tDUM; }
    else if(name == "cM-.") { return X3dna::X3dnaBPType::cMUD; }
    else if(name == "cW-m") { return X3dna::X3dnaBPType::cWUm; }
    else if(name == "t.+m") { return X3dna::X3dnaBPType::tDPm; }
    else if(name == "tM-.") { return X3dna::X3dnaBPType::tMUD; }
    else if(name == "cm+W") { return X3dna::X3dnaBPType::cmPW; }
    else if(name == "cM+M") { return X3dna::X3dnaBPType::cMPM; }
    else if(name == "cm+.") { return X3dna::X3dnaBPType::cmPD; }
    else if(name == "cm-.") { return X3dna::X3dnaBPType::cmUD; }
    else if(name == "c.-.") { return X3dna::X3dnaBPType::cDUD; }
    else if(name == "cW+W") { return X3dna::X3dnaBPType::cWPW; }
    else if(name == "t.-.") { return X3dna::X3dnaBPType::tDUD; }
    else if(name == "t.+W") { return X3dna::X3dnaBPType::tDPW; }
    else if(name == "tm-m") { return X3dna::X3dnaBPType::tmUm; }
    else if(name == "cW+.") { return X3dna::X3dnaBPType::cWPD; }
    else if(name == "tm+.") { return X3dna::X3dnaBPType::tmPD; }
    else if(name == "t.+.") { return X3dna::X3dnaBPType::tDPD; }
    else if(name == "c.+.") { return X3dna::X3dnaBPType::cDPD; }
    else if(name == "t.-m") { return X3dna::X3dnaBPType::tDUm; }
    else if(name == "t.+M") { return X3dna::X3dnaBPType::tDPM; }
    else                    {throw X3dnaException("cannot get x3dna type with: " + name); }
}

































