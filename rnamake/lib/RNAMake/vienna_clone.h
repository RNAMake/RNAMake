//
//  vienna_clone.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__vienna_clone__
#define __RNAMake__vienna_clone__
#define TURN 3

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "energy_par.h"
#include "pair_mat.h"


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
        temp = 37.0;
        ptype   = String();
        ptype.resize(size_);
        base_pair2 = bondTs(size_);
        idx = Ints(size_);
        S = Shorts(size_);
        S1 = Shorts(size_);
        BP = Ints(size_);
        params = paramT();
        setup_part_func();
    }
    
    inline
    void
    setup_part_func() {
        q = Floats(size_);
        
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
            S.resize(size+2);
            S1.resize(size+2);
            BP.resize(size+2);
        }
        
        size_ = size;
    }
    
    void
    init_fold(int length) {
        get_arrays(length);
        get_indx();
        update_fold_params_par();
        make_pair_matrix();
    }
    
    float
    fold(String const &);
    
    inline
    void
    get_indx() {
        idx[0] = 0;
        for (int i = 1; i <= size_; i++) {
            idx[i] = (i*(i-1)) >> 1;
        }
    }
    
    void
    update_fold_params_par() {
        unsigned int i,j,k,l;
        double tempf;
        model_detailsT md = params.model_details;
        params.temperature = temp;
        tempf                 = ((temp+K0)/Tmeasure);

        for (i=0; i<31; i++) {
            params.hairpin[i]  = hairpindH[i] - (hairpindH[i] - hairpin37[i])*tempf;
        }
        
        for (i=0; i<=MAXLOOP; i++) {
            params.bulge[i]          = bulgedH[i] - (bulgedH[i] - bulge37[i]) * tempf;
            params.internal_loop[i]  = internal_loopdH[i] - (internal_loopdH[i] - internal_loop37[i]) * tempf;
        }
        
        params.lxc = lxc37*tempf;
        for (i=0; i<=MAXLOOP; i++) {
            params.bulge[i] = params.bulge[30]+(int)(params.lxc*log((double)(i)/30.));
            params.internal_loop[i] = params.internal_loop[30]+(int)(params.lxc*log((double)(i)/30.));
        }
        
        params.ninio[2] = niniodH - (niniodH - ninio37) * tempf;
        
        params.TripleC = TripleCdH - (TripleCdH - TripleC37) * tempf;
        params.MultipleCA = MultipleCAdH - (MultipleCAdH - MultipleCA37) * tempf;
        params.MultipleCB = MultipleCBdH - (MultipleCBdH - MultipleCB37) * tempf;
        

        
        for (i=0; (i*7)<strlen(Tetraloops); i++) {
            params.Tetraloop_E[i] = TetraloopdH[i] - (TetraloopdH[i]-Tetraloop37[i])*tempf;
        }
        for (i=0; (i*5)<strlen(Triloops); i++) {
            params.Triloop_E[i] =  TriloopdH[i] - (TriloopdH[i]-Triloop37[i])*tempf;
        }
        for (i=0; (i*9)<strlen(Hexaloops); i++) {
            params.Hexaloop_E[i] =  HexaloopdH[i] - (HexaloopdH[i]-Hexaloop37[i])*tempf;
        }
        
        params.TerminalAU = TerminalAUdH - (TerminalAUdH - TerminalAU37) * tempf;
        params.DuplexInit = DuplexInitdH - (DuplexInitdH - DuplexInit37) *tempf;
        params.MLbase = ML_BASEdH - (ML_BASEdH - ML_BASE37) * tempf;

        for (i=0; i<=NBPAIRS; i++)
            params.MLintern[i] = ML_interndH - (ML_interndH - ML_intern37) * tempf;
        
        params.MLclosing = ML_closingdH - (ML_closingdH - ML_closing37) * tempf;
        
        /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
        for (i=0; i<=NBPAIRS; i++) {
            for (j=0; j<=NBPAIRS; j++) {
                params.stack[i][j] = stackdH[i][j] - (stackdH[i][j] - stack37[i][j])*tempf;
            }
        }
        
        /* mismatches */
        for (i=0; i<=NBPAIRS; i++) {
            for (j=0; j<5; j++) {
                for (k=0; k<5; k++) {
                    int mm;
                    params.mismatchI[i][j][k]    = mismatchIdH[i][j][k] - (mismatchIdH[i][j][k] - mismatchI37[i][j][k])*tempf;
                    params.mismatchH[i][j][k]    = mismatchHdH[i][j][k] - (mismatchHdH[i][j][k] - mismatchH37[i][j][k])*tempf;
                    params.mismatch1nI[i][j][k]  = mismatch1nIdH[i][j][k]-(mismatch1nIdH[i][j][k]-mismatch1nI37[i][j][k])*tempf;/* interior nx1 loops */
                    params.mismatch23I[i][j][k]  = mismatch23IdH[i][j][k]-(mismatch23IdH[i][j][k]-mismatch23I37[i][j][k])*tempf;/* interior 2x3 loops */
                    if(md.dangles){
                        mm                      = mismatchMdH[i][j][k] - (mismatchMdH[i][j][k] - mismatchM37[i][j][k])*tempf;
                        params.mismatchM[i][j][k]    = (mm > 0) ? 0 : mm;
                        mm                      = mismatchExtdH[i][j][k] - (mismatchExtdH[i][j][k] - mismatchExt37[i][j][k])*tempf;
                        params.mismatchExt[i][j][k]  = (mm > 0) ? 0 : mm;
                    }
                    else{
                        params.mismatchM[i][j][k] = params.mismatchExt[i][j][k] = 0;
                    }
                }
            }
        }
        
        /* dangles */
        for (i=0; i<=NBPAIRS; i++) {
            for (j=0; j<5; j++) {
                int dd;
                dd = dangle5_dH[i][j] - (dangle5_dH[i][j] - dangle5_37[i][j])*tempf;
                params.dangle5[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
                dd = dangle3_dH[i][j] - (dangle3_dH[i][j] - dangle3_37[i][j])*tempf;
                params.dangle3[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
            }
        }
    
        /* interior 1x1 loops */
        for (i=0; i<=NBPAIRS; i++) {
            for (j=0; j<=NBPAIRS; j++) {
                for (k=0; k<5; k++) {
                    for (l=0; l<5; l++) {
                        params.int11[i][j][k][l] = int11_dH[i][j][k][l] - (int11_dH[i][j][k][l] - int11_37[i][j][k][l])*tempf;
                    }
                }
            }
        }
        
        /* interior 2x1 loops */
        for (i=0; i<=NBPAIRS; i++) {
            for (j=0; j<=NBPAIRS; j++) {
                for (k=0; k<5; k++) {
                    for (l=0; l<5; l++) {
                        int m;
                        for (m=0; m<5; m++)
                            params.int21[i][j][k][l][m] = int21_dH[i][j][k][l][m] - (int21_dH[i][j][k][l][m] - int21_37[i][j][k][l][m])*tempf;
                    }
                }
            }
        }
        
        /* interior 2x2 loops */
        for (i=0; i<=NBPAIRS; i++) {
            for (j=0; j<=NBPAIRS; j++) {
                for (k=0; k<5; k++) {
                    for (l=0; l<5; l++) {
                        int m,n;
                        for (m=0; m<5; m++) {
                            for (n=0; n<5; n++) {
                                params.int22[i][j][k][l][m][n] = int22_dH[i][j][k][l][m][n] - (int22_dH[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*tempf;
                            }
                        }
                    }
                }
            }
        }
        
        strncpy(params.Tetraloops, Tetraloops, 281);
        strncpy(params.Triloops, Triloops, 241);
        strncpy(params.Hexaloops, Hexaloops, 361);
        
        params.id++;
        
    }
    
private:
    
    void
    make_ptypes(
        Shorts const &,
        String const &);
    
    
private:
    //variables from fold
    Ints c, fML, fM1, f5, f53, cc, cc1, Fmi, DMLi, DMLi1, DMLi2, DMLi_a, DMLi_o, DMLi1_a;
    Ints DMLi1_o, DMLi2_a, DMLi2_o;
    Ints idx, BP;
    String ptype, structure;
    bondTs base_pair2;
    paramT params;
    Shorts S, S1;
    double temp;
    int size_;
    //varibles from part_func
    Floats q, qb, qm, qm1, qm2, probs, q1k, qln, qq, qq1, qqm, qqm1, prm_l, prm_l1;
    Floats expMLbase, scale, Gj, Gj1;
    Ints my_iindx, iindx, jindx;
    
    
};


#endif /* defined(__RNAMake__vienna_clone__) */
