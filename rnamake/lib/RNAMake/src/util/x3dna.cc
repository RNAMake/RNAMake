//
//  x3dna.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

//RNAMake Headers
#include "util/x3dna.h"
#include "util/file_io.h"

X3dna::X3dna() {
    String env = "X3DNA="+x3dna_path();
    char * s = strdup(env.c_str());
    putenv(s);
    
    bin_path_ = x3dna_path()+"/bin/";
}

void
X3dna::generate_ref_frame(String const & path) {
    
    String full_path = "";
    String fname = filename(path);
    if     (file_exists(path + ".pdb")) {
        full_path = path + ".pdb";
    }
    else if(file_exists(path + "/" + fname + ".pdb")) {
        full_path = path + "/" + fname + ".pdb";
    }
    else {
        throw "cannot find pdb for generate_ref_frames\n";
    }
    
    String find_pair_path = bin_path_ + "find_pair ";
    String analyze_path   = bin_path_ + "analyze ";
    String command = find_pair_path + full_path + " 2> /dev/null stdout | " + analyze_path + "stdin >& /dev/null";
    char * s = strdup(command.c_str());
    int result = std::system(s);

    if(result != 0) { throw "could not call x3dna properly\n"; }
    
    String files_str = "auxiliary.par,bestpairs.pdb,bp_helical.par,bp_order.dat,bp_step.par,cf_7methods.par,col_chains.scr,col_helices.scr,hel_regions.pdb,hstacking.pdb,poc_haxis.r3d,stacking.pdb";
    Strings files = split_str_by_delimiter(files_str, ",");
    files.push_back(fname+".out");
    for(auto const & f : files) {
        try { std::remove(f.c_str()); }
        catch(String const & e) {}
    }
    
}

void
X3dna::generate_dssr_file(String const & path) {
    
    String full_path = "";
    String fname = filename(path);
    if     (file_exists(path + ".pdb")) {
        full_path = path + ".pdb";
    }
    else if(file_exists(path + "/" + fname + ".pdb")) {
        full_path = path + "/" + fname + ".pdb";
    }
    else {
        throw "cannot find pdb for generate_ref_frames\n";
    }
    
    String dssr_path = bin_path_ + "x3dna-dssr ";
    String command = dssr_path + "-i="+full_path+" -o="+fname+"_dssr.out --non-pair >& /dev/null";
    char * s = strdup(command.c_str());
    int result = std::system(s);

    if(result != 0) { throw "could not call x3dna properly\n"; }
    
    String filename_str = "dssr-2ndstrs.ct,dssr-2ndstrs.dbn,dssr-helices.pdb,dssr-pairs.pdb,dssr-stems.pdb,hel_regions.pdb,hstacking.pdb,poc_haxis.r3d,stacking.pdb,dssr-torsions.dat,dssr-Kturns.pdb,dssr-multiplets.pdb,dssr-hairpins.pdb,dssr-Aminors.pdb";
    Strings files = split_str_by_delimiter(filename_str, ",");
    for(auto const & f : files) {
        try { std::remove(f.c_str()); }
        catch(String const & e) {}
    }

}

String
X3dna::_get_ref_frame_path(String const & pdb_name) {
    
    String basedir = base_dir(pdb_name);
    String ref_frames_path = "";
    // is pdb
    if     (file_exists(basedir + "/ref_frames.dat")) {
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
X3dna::_get_dssr_file_path(String const & pdb_path) {
    String basedir = base_dir(pdb_path);
    String pdb_name = filename(pdb_path);
    String dssr_name = pdb_name + "_dssr.out";
    String dssr_file_path = "";
    
    if     (file_exists(basedir + "/" + dssr_name)) {
        dssr_file_path = basedir + "/" + dssr_name;
    }
    
    else if(file_exists(pdb_path + "/" + dssr_name )) {
        dssr_file_path = pdb_path + "/" + dssr_name;
    }
    else {
        generate_dssr_file(pdb_name);
        dssr_file_path = dssr_name;
        
    }
    return dssr_file_path;
}

Point
X3dna::_convert_strings_to_point(Strings const & spl) {
    std::vector<float> f;
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
            X3Residue res1(res_num_1, res_id_1, "");
            X3Residue res2(res_num_2, res_id_2, "");
            X3Basepair bp(res1, res2, r, d);
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

X3Residue
X3dna::_parse_dssr_res_str(String const & res_str) {
    Strings spl = split_str_by_delimiter(res_str, ".");
    String chain = spl[0];
    String rnum = spl[1].substr(1);
    int num;
    try{
        num = std::stoi(rnum);
    } catch(...) {
        std::cout << res_str << std::endl;
        throw "could not parse " + res_str + " into a residue\n";
    }

    return X3Residue(num, chain, "");
    
}

X3Basepairs const &
X3dna::get_basepairs(String const & pdb_path) {
    String ref_frames_path = _get_ref_frame_path(pdb_path);
    String dssr_file_path  = _get_dssr_file_path(pdb_path);
    
    _parse_ref_frame_file(ref_frames_path);
    std::map<String, Strings> sections = _divide_dssr_file_into_sections(dssr_file_path);
    Strings section = sections["base"];
    int found = 0;
    
    for(auto const & l : section) {
        if(l.length() < 3) { continue; }
        Strings spl = _split_over_white_space(l);
        Strings spl2 = split_str_by_delimiter(spl[1], ".");
        if(spl2.size() == 1) {continue; }
        if(spl.size() < 6) { continue; }
        X3Residue res1 = _parse_dssr_res_str(spl[1]);
        X3Residue res2 = _parse_dssr_res_str(spl[2]);
        String bp_type = spl[7];
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
            basepairs_.push_back(X3Basepair(res1, res2, Matrix(), Point(-1,-1,-1)));
        }
        
    }    
    return basepairs_;
    


}









