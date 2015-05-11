//
//  x3dna.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

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
        ref_frames_path = "ref_Frames.dat";
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











