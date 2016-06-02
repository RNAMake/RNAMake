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


struct X3Residue {
    X3Residue(
        int nnum,
        String const & nchain_id,
        String const & ni_code):
    num(nnum),
    chain_id(nchain_id),
    i_code(ni_code)
    {}
    
    ~X3Residue() {}
    
    bool
    operator == (X3Residue const & r) const {
        if(num == r.num && chain_id.compare(r.chain_id) == 0 && i_code.compare(r.i_code) == 0) {
            return 1;
        }
        else {
            return 0;
        }
    }
    
public:
    int num;
    String chain_id, i_code;
};

struct X3Basepair {
public:
    X3Basepair(
        X3Residue const & nres1,
        X3Residue const & nres2,
        Matrix const & nr,
        Point const & nd):
    res1(nres1),
    res2(nres2),
    r(nr),
    d(nd),
    bp_type("c...")
    {}
    
    ~X3Basepair() {}

public:
    X3Residue res1, res2;
    Point d;
    Matrix r;
    String bp_type;
    
};

typedef std::vector<X3Residue>  X3Residues;
typedef std::vector<X3Basepair> X3Basepairs;

struct X3Motif {
    X3Motif(
        X3Residues const & nresidues,
        String const & nmtype):
    residues(nresidues),
    mtype(nmtype)
    {}
    
    X3Residues residues;
    String mtype;
};

typedef std::vector<X3Motif>    X3Motifs;

class X3dna {
public:
    
    X3dna();
    
    ~X3dna() {}

public:
    void
    generate_ref_frame(String const &);
    
    void
    generate_dssr_file(String const &);
    
    X3Basepairs const &
    get_basepairs(String const &);
    
    X3Motifs
    get_motifs(String const &);
    
private:
    
    String
    _get_ref_frame_path(String const &);
    
    String
    _get_dssr_file_path(String const &);
    
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
    
};


#endif /* defined(__RNAMake__x3dna__) */
