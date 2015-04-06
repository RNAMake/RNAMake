//
//  vienna_clone.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__vienna_clone__
#define __RNAMake__vienna_clone__

#include <stdio.h>
#include "types.h"
#include "energy_par.h"


typedef std::vector<Ints>   Ints2d;
typedef std::vector<Ints2d> Ints3d;
typedef std::vector<Ints3d> Ints4d;
typedef std::vector<Ints4d> Ints5d;
typedef std::vector<Ints5d> Ints6d;

struct bondT {
    unsigned int i;
    unsigned int j;
};

struct model_detailsT {
    int     dangles;      /**<  \brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
                           \note   Some function do not implement all dangle model but only a subset of
                           (0,1,2,3). Read the documentaion of the particular recurrences or
                           energy evaluation function for information about the provided dangle
                           model.
                           */
    int     special_hp;   /**<  \brief  Include special hairpin contributions for tri, tetra and hexaloops */
    int     noLP;         /**<  \brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
    int     noGU;         /**<  \brief  Do not allow GU pairs */
    int     noGUclosure;  /**<  \brief  Do not allow loops to be closed by GU pair */
    int     logML;        /**<  \brief  Use logarithmic scaling for multi loops */
    int     circ;         /**<  \brief  Assume molecule to be circular */
    int     gquad;        /**<  \brief  Include G-quadruplexes in structure prediction */
    int     canonicalBPonly;  /**<  \brief  remove non-canonical bp's from constraint structures  */
};

struct paramT {
    int id;
    int stack[NBPAIRS+1][NBPAIRS+1];
    int hairpin[31];
    int bulge[MAXLOOP+1];
    int internal_loop[MAXLOOP+1];
    int mismatchExt[NBPAIRS+1][5][5];
    int mismatchI[NBPAIRS+1][5][5];
    int mismatch1nI[NBPAIRS+1][5][5];
    int mismatch23I[NBPAIRS+1][5][5];
    int mismatchH[NBPAIRS+1][5][5];
    int mismatchM[NBPAIRS+1][5][5];
    int dangle5[NBPAIRS+1][5];
    int dangle3[NBPAIRS+1][5];
    int int11[NBPAIRS+1][NBPAIRS+1][5][5];
    int int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
    int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
    int ninio[5];
    double  lxc;
    int     MLbase;
    int     MLintern[NBPAIRS+1];
    int     MLclosing;
    int     TerminalAU;
    int     DuplexInit;
    int     Tetraloop_E[200];
    char    Tetraloops[1401];
    int     Triloop_E[40];
    char    Triloops[241];
    int     Hexaloop_E[40];
    char    Hexaloops[1801];
    int     TripleC;
    int     MultipleCA;
    int     MultipleCB;
    double  temperature;            /**<  \brief  Temperature used for loop contribution scaling */
    
    model_detailsT model_details;   /**<  \brief  Model details to be used in the recursions */
    
};

typedef std::vector<bondT> bondTs;


class ViennaClone {
public:
    
    inline
    ViennaClone() {
        size_ = 1;
        c       = Ints(size_);
        fML     = Ints(size_);
        f5      = Ints(size_);
        f53     = Ints(size_);
        cc      = Ints(size_);
        cc1     = Ints(size_);
        Fmi     = Ints(size_);
        DMLi    = Ints(size_);
        DMLi1   = Ints(size_);
        DMLi2   = Ints(size_);
        DMLi_a  = Ints(size_);
        DMLi_o  = Ints(size_);
        DMLi1_a = Ints(size_);
        DMLi1_o = Ints(size_);
        DMLi2_a = Ints(size_);
        DMLi2_o = Ints(size_);
        ptype   = String();
        ptype.resize(size_);
        base_pair2 = bondTs(size_);
        idx = Ints(size_);
        parameters = paramT();
    }
    
    inline
    void
    get_arrays(int size) {
        if(size > size_) {
            c.resize(size);
            fML.resize(size);
            f5.resize(size);
            cc.resize(size);
            cc1.resize(size);
            Fmi.resize(size);
            DMLi.resize(size);
            DMLi1.resize(size);
            DMLi2.resize(size);
            DMLi_a.resize(size);
            DMLi_o.resize(size);
            DMLi1_a.resize(size);
            DMLi1_o.resize(size);
            DMLi2_a.resize(size);
            DMLi2_o.resize(size);
            ptype.resize(size);
            base_pair2.resize(size);
            idx.resize(size);
            
        }
        
        size_ = size;
    }
    
    void
    init_fold(int length) {
        get_arrays(length);
        get_indx();
    }
    
    inline
    void
    get_indx() {
        idx[0] = 0;
        for (int i = 1; i <= size_; i++) {
            idx[i] = (i*(i-1)) >> 1;
        }
    }
    
    
private:
    Ints c, fML, fM1, f5, f53, cc, cc1, Fmi, DMLi, DMLi1, DMLi2, DMLi_a, DMLi_o, DMLi1_a;
    Ints DMLi1_o, DMLi2_a, DMLi2_o;
    Ints idx;
    String ptype;
    bondTs base_pair2;
    paramT parameters;
    int size_;
    
    
};


#endif /* defined(__RNAMake__vienna_clone__) */
