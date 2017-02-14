//
//  x3dna.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__x3dna__
#define __RNAMake__x3dna__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <stdexcept>


//RNAMake Headers
#include "base/string.h"
#include "base/settings.h"
#include "math/xyz_vector.h"
#include "math/xyz_matrix.h"

class X3dnaException : public std::runtime_error {
public:
    X3dnaException(
        String const & message) :
    std::runtime_error(message)
    {}
    
};

class X3dna {
public:
    
    X3dna();
    
    ~X3dna() { delete s_; }

public:
    //- = U
    //+ = P
    //. = D
    enum X3dnaBPType {
        cmU  = 0,  //cm-
        cMUM = 1,  //cM-M
        tWPW = 2,  //tW+W
        cDPM = 3,  //c.+M
        DWPW = 4,  //.W+W
        tWUM = 5,  //tW-M
        tmUM = 6,  //tm-M
        cWPM = 7,  //cW+M
        DWUW = 8,  //.W-W
        cMPD = 9,  //cM+.
        cDUm = 10, //c.-m
        cMPW = 11, //cM+W
        tMPm = 12, //tM+m
        tMUW = 13, //tM-W
        cmUm = 14, //cm-m
        cMUW = 15, //cM-W
        cWUW = 16, //cW-W
        cDUM = 17, //c.-M
        cmPM = 18, //cm+M
        cmUM = 19, //cm-M
        DDDD = 20, //....
        cmUW = 21, //cm-W
        tMUm = 22, //tM-m
        cDUW = 23, //c.-W
        cMPm = 24, //cM+m
        cMUm = 25, //cM-m
        cDDD = 26, //c...
        tWPm = 27, //tW+m
        cDPm = 28, //c.+m
        tmPm = 29, //tm+m
        tWPD = 30, //tW+.
        tmPW = 31, //tm+W
        tDDD = 32, //t...
        cWUD = 33, //cW-.
        cWUM = 34, //cW-M
        tDUW = 35, //t.-W
        tMPM = 36, //tM+M
        tDUM = 37, //t.-M
        cMUD = 38, //cM-.
        cWUm = 39, //cW-m
        tDPm = 40, //t.+m
        tMUD = 41, //tM-.
        cmPW = 42, //cm+W
        cMPM = 43, //cM+M
        cmPD = 44, //cm+.
        cmUD = 45, //cm-.
        cDUD = 46, //c.-.
        cWPW = 47, //cW+W
        tDUD = 48, //t.-.
        tDPW = 49, //t.+W
        tmUm = 50, //tm-m
        cWPD = 51, //cW+.
        tmPD = 52, //tm+.
        tDPD = 53, //t.+.
        cDPD = 54, //c.+.
        tDUm = 55, //t.-m
        tDPM = 56, //t.+M
    };

public:
    struct X3Residue {
        bool
        operator == (X3Residue const & r) const {
            if(num == r.num && chain_id == r.chain_id && i_code == r.i_code) { return 1; }
            else { return 0; }
        }

        int num;
        char chain_id, i_code;
    };

    typedef std::vector<X3Residue>  X3Residues;

    struct X3Basepair {
        X3Residue res1, res2;
        Point d;
        Matrix r;
        X3dnaBPType bp_type;
    };

    typedef std::vector<X3Basepair> X3Basepairs;

    struct X3Motif {
        X3Residues residues;
        String mtype;
    };

    typedef std::vector<X3Motif>    X3Motifs;

public:
    void
    generate_ref_frame(String const &);
    
    void
    generate_dssr_file(String const &);
    
    X3Basepairs const &
    get_basepairs(String const & pdb_path,
                  bool force_build_files = false);
    
    X3Motifs
    get_motifs(String const &);
    
private:
    
    String
    _get_ref_frame_path(String const & pdb_path,
                        bool force_build_files = false);
    
    String
    _get_dssr_file_path(String const & pdb_path,
                        bool force_build_files = false);
    
    Point
    _convert_strings_to_point(Strings const &);
    
    void
    _parse_ref_frame_file(String const &);
    
    std::map<String, Strings>
    _divide_dssr_file_into_sections(String const &);
    
    Strings
    _split_over_white_space(String const &);
    
    X3Residue
    _parse_dssr_res_str(String const &);
    
    X3Motifs
    _parse_dssr_section(
        Strings const &,
        String const &);
    
    X3Motifs
    _parse_dssr_helix_section(
        Strings const &);
    
    
private:
    String bin_path_;
    X3Basepairs basepairs_;
    char * s_;
};

X3dna::X3dnaBPType
get_x3dna_by_type(String const &);

#endif /* defined(__RNAMake__x3dna__) */
