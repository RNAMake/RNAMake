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
#include <regex>


//RNAMake Headers
#include <base/string.h>
#include <math/xyz_vector.h>
#include <math/xyz_matrix.h>

namespace util {

class X3dnaException : public std::runtime_error {
public:
    X3dnaException(
            String const & message) :
            std::runtime_error(message) {}

};

//- = U
//+ = P
//. = D
enum class X3dnaBPType {
    cmU = 0,  //cm-
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

class X3dna {
public:

    X3dna();

    ~X3dna() { delete s_; }

public:
    struct X3Residue {
        inline
        X3Residue(
                int nnum=-1,
                char nchain_id=' ',
                char ni_code=' '):
                num(nnum),
                chain_id(nchain_id),
                i_code(ni_code) {}


        inline
        bool
        operator == (
                X3Residue const & r) const {
            if (num == r.num && chain_id == r.chain_id && i_code == r.i_code) { return 1; }
            else { return 0; }
        }

        int num;
        char chain_id, i_code;
    };

    typedef std::vector<X3Residue> X3Residues;

    struct X3Basepair {
        X3Residue res1, res2;
        math::Point d;
        math::Matrix r;
        X3dnaBPType bp_type;
    };

    typedef std::vector<X3Basepair> X3Basepairs;

    struct X3Motif {
        X3Residues residues;
        String mtype;
    };

    typedef std::vector<X3Motif> X3Motifs;

    struct X3BPInfo {
        inline
        X3BPInfo(
                std::smatch const & match) {
            auto i = -1;
            for(auto const & v: match) {
                i++;
                if     (i == 0) { continue; }
                else if(i == 1) { res1_chain_id = String(v)[0]; }
                else if(i == 2) { res1_num = std::stoi(String(v)); }
                else if(i == 3) { res1_type_name = String(v); }
                else if(i == 4) { res1_name = String(v)[0]; }
                else if(i == 5) { res2_chain_id = String(v)[0];}
                else if(i == 6) { res2_num = std::stoi(String(v)); }
                else if(i == 7) { res2_type_name = String(v); }
                else if(i == 8) { res2_name = String(v)[0]; }
            }
            if(i != 8) {
                throw X3dnaException("could not properly parse X3dna ref frame line:");
            }
        }
        String res1_type_name, res2_type_name;
        char res1_name, res2_name;
        char res1_chain_id, res2_chain_id;
        int res1_num, res2_num;

    };

public:
    X3Basepairs
    get_basepairs(
            String const &) const;
    
    X3Basepairs
    get_basepairs_json(
            String const &) const;
    

    X3Motifs
    get_motifs(
            String const &) const;

public:
    void
    set_rebuild_files(
            bool rebuild_files) const{ rebuild_files_ = rebuild_files; }

private:
    void
    _generate_ref_frame(
            String const &) const;

    void
    generate_dssr_file(
            String const &) const;

private:
    void
    _delete_files(
            Strings const &) const;

    void
    _delete_file(
            String const &) const;

    math::Point
    _convert_string_to_point(
            String const &) const;

    void
    _parse_ref_frame_file(
            String const &,
            X3Basepairs &) const;

    std::map<String, Strings>
    _parse_dssr_file_into_sections(
            String const &) const;

    Strings
    _split_over_white_space(
            String const &) const;

    X3Residue *
    _parse_dssr_res_str(
            String const &) const;

    X3Motifs
    _parse_dssr_section(
            Strings const &,
            String const &) const;

    X3Motifs
    _parse_dssr_helix_section(
            Strings const &) const;


private:
    String bin_path_;
    char *s_;
    Strings ref_frame_files_to_delete_, dssr_files_to_delete_;
    // flags to decide whether to build files or not
    mutable bool rebuild_files_;
    mutable bool generated_ref_frames_, generated_dssr_;
    mutable bool no_ref_frames_;
};

X3dnaBPType
get_x3dna_by_type(String const &);

String
get_str_from_x3dna_type(X3dnaBPType);

void
compare_bps(X3dna::X3Basepairs&, X3dna::X3Basepairs&);

}

#endif /* defined(__RNAMake__x3dna__) */
