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
#define MAXSECTORS        500
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MIN3(A, B, C)   (MIN2(  (MIN2((A),(B))) ,(C)))

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

struct sect {
    int  i;
    int  j;
    int ml;
};


typedef std::vector<bondT> bondTs;
typedef std::vector<sect> sects;


class ViennaClone {
public:
    
    inline
    ViennaClone() {
        size_ = 10;
        c       = Ints((size_*(size_+1)/2+2));
        fML     = Ints((size_*(size_+1)/2+2));
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
        ptype   = Chars((size_*(size_+1)/2+2));
        base_pair2 = bondTs(4*(1+size_/2));
        indx = Ints(size_);
        S = Shorts(size_);
        S1 = Shorts(size_);
        BP = Ints(size_);
        params = paramT();
        setup_part_func();
        uniq_ML = 0;
        sector = sects(MAXSECTORS);
        backtrack_type = 'F';
    }
    
    inline
    void
    setup_part_func() {
        q         = Floats(size_);
        qb        = Floats(size_);
        qm        = Floats(size_);
        probs     = Floats(size_);
        q1k       = Floats(size_+1);
        qln       = Floats(size_+2);
        qq        = Floats(size_+2);
        qq1       = Floats(size_+2);
        qqm       = Floats(size_+2);
        qqm1      = Floats(size_+2);
        prm_l     = Floats(size_+2);
        prm_l1    = Floats(size_+2);
        prml      = Floats(size_+2);
        expMLbase = Floats(size_+1);
        scale     = Floats(size_+1);
        Gj        = Floats(size_+2);
        Gj1       = Floats(size_+2);
        
        my_iindx  = Ints(size_);
        iindx     = Ints(size_);
        jindx     = Ints(size_);
        
        get_iindx(my_iindx);
        get_iindx(iindx);
        get_indx(jindx);
    }
    
    inline
    void
    get_arrays(int size) {
        if(size > size_) {
            c.resize((size*(size+1)/2+2));
            fML.resize((size*(size+1)/2+2));
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
            ptype.resize((size*(size+1)/2+2));
            base_pair2.resize(4*(1+size/2));
            indx.resize(size);
            S.resize(size+2);
            S1.resize(size+2);
            BP.resize(size+2);
            get_indx(indx);
            get_arrays_part_func(size);
        }
        
        size_ = size;
    }
    
    inline
    void
    get_arrays_part_func(int size) {
        q.resize(size);
        qb.resize(size);
        qm.resize(size);
        probs.resize(size);
        q1k.resize(size+1);
        qln.resize(size+2);
        qq.resize(size+2);
        qq1.resize(size+2);
        qqm.resize(size+2);
        qqm1.resize(size+2);
        prm_l.resize(size+2);
        prm_l1.resize(size+2);
        prml.resize(size+2);
        expMLbase.resize(size+1);
        scale.resize(size+1);
        
        my_iindx.resize(size);
        iindx.resize(size);
        jindx.resize(size);
        
        get_iindx(my_iindx);
        get_iindx(iindx);
        get_indx(jindx);
    }
    
    void
    init_fold(int length) {
        params.model_details.dangles = 2;
        params.model_details.special_hp = 1;
        get_arrays(length);
        get_indx(indx);
        update_fold_params_par();
        make_pair_matrix();
    }
    
    float
    fold(String const &);
    
    inline
    void
    get_indx(Ints & c_idx) {
        c_idx[0] = 0;
        for (int i = 1; i <= size_; i++) {
            c_idx[i] = (i*(i-1)) >> 1;
        }
    }
    
    inline
    void
    get_iindx(Ints & c_idx) {
        unsigned int i;
        for (i = 1; i <= size_; i++) {
            c_idx[i] = (i*(i-1)) >> 1;        /* i(i-1)/2 */
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
        
        for (i=0; i<=MIN2(30,MAXLOOP); i++) {
            params.bulge[i]          = bulgedH[i] - (bulgedH[i] - bulge37[i]) * tempf;
            params.internal_loop[i]  = internal_loopdH[i] - (internal_loopdH[i] - internal_loop37[i]) * tempf;
        }
        
        params.lxc = lxc37*tempf;
        for (; i<=MAXLOOP; i++) {
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
    
    int
    fill_arrays(
        String const &);
    
    void
    backtrack(
        String const &,
        int);
    
    void
    parenthesis_structure (int);
    
private:
    
    void
    make_ptypes(
        Shorts const &,
        String const &);
    
    
private: //Energy calculations
    
    inline
    int
    E_Hairpin(
        int size,
        int type,
        int si1,
        int sj1,
        const char *string) {
        
        int energy = (size <= 30) ? params.hairpin[size] : params.hairpin[30]+(int)(params.lxc*log((size)/30.));
        if (params.model_details.special_hp){
            if (size == 4) { /* check for tetraloop bonus */
                char tl[7]={0}, *ts;
                strncpy(tl, string, 6);
                if ((ts=strstr(params.Tetraloops, tl)))
                    return (params.Tetraloop_E[(ts - params.Tetraloops)/7]);
            }
            else if (size == 6) {
                char tl[9]={0}, *ts;
                strncpy(tl, string, 8);
                if ((ts=strstr(params.Hexaloops, tl)))
                    return (energy = params.Hexaloop_E[(ts - params.Hexaloops)/9]);
            }
            else if (size == 3) {
                char tl[6]={0,0,0,0,0,0}, *ts;
                strncpy(tl, string, 5);
                if ((ts=strstr(params.Triloops, tl))) {
                    return (params.Triloop_E[(ts - params.Triloops)/6]);
                }
                return (energy + (type>2 ? params.TerminalAU : 0));
            }
        }
        energy += params.mismatchH[type][si1][sj1];

        return energy;
        
    }
    
    inline
    int
    E_IntLoop(
        int n1,
        int n2,
        int type,
        int type_2,
        int si1,
        int sj1,
        int sp1,
        int sq1) {
        /* compute energy of degree 2 loop (stack bulge or interior) */
        int nl, ns, energy;
        energy = INF;
        
        if (n1>n2) { nl=n1; ns=n2;}
        else {nl=n2; ns=n1;}
        
        if (nl == 0) {
            return params.stack[type][type_2];  /* stack */
        }
            
        if (ns==0) {                      /* bulge */
            energy = (nl<=MAXLOOP) ? params.bulge[nl]:
            (params.bulge[30]+(int)(params.lxc*log(nl/30.)));
            if (nl==1) energy += params.stack[type][type_2];
            else {
                if (type>2) energy += params.TerminalAU;
                if (type_2>2) energy += params.TerminalAU;
            }

            return energy;
        }
        else {                            /* interior loop */
            if (ns==1) {
                if (nl==1)                    /* 1x1 loop */
                    return params.int11[type][type_2][si1][sj1];
                if (nl==2) {                  /* 2x1 loop */
                    if (n1==1)
                        energy = params.int21[type][type_2][si1][sq1][sj1];
                    else
                        energy = params.int21[type_2][type][sq1][si1][sp1];
                    return energy;
                }
                else {  /* 1xn loop */
                    energy = (nl+1<=MAXLOOP)?(params.internal_loop[nl+1]) : (params.internal_loop[30]+(int)(params.lxc*log((nl+1)/30.)));
                    energy += MIN2(MAX_NINIO, (nl-ns)*params.ninio[2]);
                    energy += params.mismatch1nI[type][si1][sj1] + params.mismatch1nI[type_2][sq1][sp1];
                    return energy;
                }
            }
            else if (ns==2) {
                if(nl==2)      {              /* 2x2 loop */
                    return params.int22[type][type_2][si1][sp1][sq1][sj1];}
                else if (nl==3){              /* 2x3 loop */
                    energy = params.internal_loop[5]+params.ninio[2];
                    energy += params.mismatch23I[type][si1][sj1] + params.mismatch23I[type_2][sq1][sp1];
                    return energy;
                }
                
            }
            { /* generic interior loop (no else here!)*/
                energy = (n1+n2<=MAXLOOP)?(params.internal_loop[n1+n2]) : (params.internal_loop[30]+(int)(params.lxc*log((n1+n2)/30.)));
                
                energy += MIN2(MAX_NINIO, (nl-ns)*params.ninio[2]);
                
                energy += params.mismatchI[type][si1][sj1] + params.mismatchI[type_2][sq1][sp1];
            }
        }
        return energy;
    }
    
    inline
    int
    E_MLstem(
        int type,
        int si1,
        int sj1) {
        
        int energy = 0;
        if(si1 >= 0 && sj1 >= 0){
            energy += params.mismatchM[type][si1][sj1];
        }
        else if (si1 >= 0){
            energy += params.dangle5[type][si1];
        }
        else if (sj1 >= 0){
            energy += params.dangle3[type][sj1];
        }
        
        if(type > 2)
            energy += params.TerminalAU;
        
        energy += params.MLintern[type];
        
        return energy;
    }

    inline
    int
    E_ExtLoop(
        int type,
        int si1,
        int sj1) {
        
        int energy = 0;
        if(si1 >= 0 && sj1 >= 0){
            energy += params.mismatchExt[type][si1][sj1];
        }
        else if (si1 >= 0){
            energy += params.dangle5[type][si1];
        }
        else if (sj1 >= 0){
            energy += params.dangle3[type][sj1];
        }
        
        if(type > 2)
            energy += params.TerminalAU;
        
        return energy;
        
    }
    

    
    
private:
    //variables from fold
    Ints c, fML, fM1, f5, f53, cc, cc1, Fmi, DMLi, DMLi1, DMLi2, DMLi_a, DMLi_o, DMLi1_a;
    Ints DMLi1_o, DMLi2_a, DMLi2_o;
    Ints indx, BP;
    String structure;
    Chars ptype;
    bondTs base_pair2;
    paramT params;
    Shorts S, S1;
    double temp;
    int uniq_ML;
    int size_;
    sects sector;
    char backtrack_type;
    //varibles from part_func
    Floats q, qb, qm, probs, q1k, qln, qq, qq1, qqm, qqm1, prml, prm_l, prm_l1;
    Floats expMLbase, scale, Gj, Gj1;
    Ints my_iindx, iindx, jindx;
    
    
};


#endif /* defined(__RNAMake__vienna_clone__) */
