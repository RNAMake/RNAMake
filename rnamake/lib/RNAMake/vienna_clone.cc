//
//  vienna_clone.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "vienna_clone.h"
#include "pair_mat.h"

#define LOCALITY          0.      /* locality parameter for base-pairs */

float
ViennaClone::fold(
    String const & string) {
    
    //update length of data arrays if size is greater then previous
    get_arrays((int)string.length());
    encode_sequence(string, S, 0);
    encode_sequence(string, S1, 1);
    make_ptypes(S, structure);
    
    fill_arrays(string);
    
    return 1;
}

void
ViennaClone::make_ptypes(
    Shorts const & S,
    String const & structure) {
    
    int n,i,j,k,l, noLP;
    noLP = params.model_details.noLP;
    n = S[0];
    
    n=S[0];
    for (k=1; k<n-TURN; k++) {
        for (l=1; l<=2; l++) {
            int type,ntype=0,otype=0;
            i=k; j = i+TURN+l; if (j>n) continue;
            type = BP_pair[S[i]][S[j]];
            while ((i>=1)&&(j<=n)) {
                if ((i>1)&&(j<n)) ntype = BP_pair[S[i-1]][S[j+1]];
                ptype[indx[j]+i] = (char) type;
                otype =  type;
                type  = ntype;
                i--; j++;
            }
        }
    }
}

int
ViennaClone::fill_arrays(
    String const & string)  {
    
    const char * s = string.c_str();
    
    int   i, j, k, length, energy, en, mm5, mm3;
    int   decomp, new_fML;
    int   no_close, type, type_2, tt, max_separation;
    int   bonus=0;
    int   dangle_model, noGUclosure, with_gquads;
    int   noLonelyPairs=0;
    int   circular=0;

    dangle_model  = params.model_details.dangles;
    noGUclosure   = params.model_details.noGUclosure;
    length = string.length();
    
    max_separation = (int) ((1.-LOCALITY)*(double)(length-2)); /* not in use */
    
    for (j=1; j<=length; j++) {
        Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
    }
    
    for (j = 1; j<=length; j++) {
        for (i=(j>TURN?(j-TURN):1); i<j; i++) {
            c[indx[j]+i] = fML[indx[j]+i] = INF;
            if (uniq_ML) fM1[indx[j]+i] = INF;
        }
    }
    
    if (length <= TURN) return 0;
    
    for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
        for (j = i+TURN+1; j <= length; j++) {

            int p, q, ij, jj, ee;
            int minq, maxq, l1, up, c0, c1, c2, c3;
            int MLenergy;

            ij = indx[j]+i;
            bonus = 0;
            type = ptype[ij];
            energy = INF;
            no_close = (((type==3)||(type==4))&&noGUclosure&&(bonus==0));
            
            if (j-i-1 > max_separation) type = 0;  /* forces locality degree */

            if (type) {   /* we have a pair */
                int new_c=0, stackEnergy=INF;
                /* hairpin ----------------------------------------------*/
                new_c = (no_close) ? FORBIDDEN : E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], s+i-1);
                
                /*--------------------------------------------------------
                 check for elementary structures involving more than one
                 closing pair.
                 --------------------------------------------------------*/
                
                for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
                    minq = j-i+p-MAXLOOP-2;
                    if (minq<p+1+TURN) minq = p+1+TURN;
                    for (q = minq; q < j; q++) {
                        type_2 = ptype[indx[q]+p];
                        
                        if (type_2==0) continue;
                        type_2 = rtype[type_2];
                        
                        if (noGUclosure) {
                            if (no_close||(type_2==3)||(type_2==4)) {
                                if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */
                            }
                        }
                        
                        energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                                           S1[i+1], S1[j-1], S1[p-1], S1[q+1]);
                        
                        ee = energy+c[indx[q]+p];
                        std::cout << ij << " " << ee << " " << energy << std::endl;
                        new_c = MIN2(new_c, ee);
                        if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */
                        
                    } /* end q-loop */
                } /* end p-loop */
                /* multi-loop decomposition ------------------------*/

                if (!no_close) {
                    decomp = DMLi1[j-1];
                    tt = rtype[type];
                    switch(dangle_model){
                            /* no dangles */
                        case 0:   decomp += E_MLstem(tt, -1, -1);
                            break;
                            
                            /* double dangles */
                        case 2:   decomp += E_MLstem(tt, S1[j-1], S1[i+1]);
                            break;
                            
                            /* normal dangles, aka dangles = 1 || 3 */
                        default:  decomp += E_MLstem(tt, -1, -1);
                            decomp = MIN2(decomp, DMLi2[j-1] + E_MLstem(tt, -1, S1[i+1]) + params.MLbase);
                            decomp = MIN2(decomp, DMLi2[j-2] + E_MLstem(tt, S1[j-1], S1[i+1]) + 2*params.MLbase);
                            decomp = MIN2(decomp, DMLi1[j-2] + E_MLstem(tt, S1[j-1], -1) + params.MLbase);
                            break;
                    }
                    MLenergy = decomp + params.MLclosing;
                    new_c = MIN2(new_c, MLenergy);
                }
                
                new_c = MIN2(new_c, cc1[j-1]+stackEnergy);
                cc[j] = new_c + bonus;
                //std::cout << ij << " " << cc[j] << " " << cc1[j-1] << " " << stackEnergy << std::endl;
                if (noLonelyPairs)
                    c[ij] = cc1[j-1]+stackEnergy+bonus;
                else
                    c[ij] = cc[j];
                
            } /* end >> if (pair) << */
            
            else c[ij] = INF;
            
            /* done with c[i,j], now compute fML[i,j] and fM1[i,j] */
            
            /* (i,j) + MLstem ? */
            new_fML = INF;
            if(type){
                new_fML = c[ij];
                switch(dangle_model){
                    case 2:
                        new_fML += E_MLstem(type, (i==1) ? S1[length] : S1[i-1], S1[j+1]);
                        //std::cout << type << " " << ((i==1) ? S1[length] : S1[i-1]) << " " << S1[j+1] << " " << i << " " << j << " " << ij << " " << new_fML <<  " " << c[ij] << std::endl;
                        break;
                    default:
                        new_fML += E_MLstem(type, -1, -1);
                        break;
                }
            }
            
            if (uniq_ML){
                fM1[ij] = MIN2(fM1[indx[j-1]+i] + params.MLbase, new_fML);
            }
            
            /* free ends ? -----------------------------------------*/
            /*  we must not just extend 3'/5' end by unpaired nucleotides if
             *   dangle_model == 1, this could lead to d5+d3 contributions were
             *   mismatch must be taken!
             */
            switch(dangle_model){
                    /* no dangles */
                case 0:
                    new_fML = MIN2(new_fML, fML[ij+1]+params.MLbase);
                    new_fML = MIN2(fML[indx[j-1]+i]+params.MLbase, new_fML);
                    break;
                    
                    /* double dangles */
                case 2:
                    //std::cout << ij << " " << fML[ij+1]+params.MLbase << " " << fML[indx[j-1]+i]+params.MLbase << std::endl;
                    new_fML = MIN2(new_fML, fML[ij+1]+params.MLbase);
                    new_fML = MIN2(fML[indx[j-1]+i]+params.MLbase, new_fML);
                    break;
                    
                    /* normal dangles, aka dangle_model = 1 || 3 */
                default:
                    mm5 = ((i>1) || circular) ? S1[i] : -1;
                    mm3 = ((j<length) || circular) ? S1[j] : -1;
                    new_fML = MIN2(new_fML, fML[ij+1] + params.MLbase);
                    new_fML = MIN2(new_fML, fML[indx[j-1]+i] + params.MLbase);
                    tt = ptype[ij+1];
                    if(tt) new_fML = MIN2(new_fML, c[ij+1] + E_MLstem(tt, mm5, -1) + params.MLbase);
                    tt = ptype[indx[j-1]+i];
                    if(tt) new_fML = MIN2(new_fML, c[indx[j-1]+i] + E_MLstem(tt, -1, mm3) + params.MLbase);
                    tt = ptype[indx[j-1]+i+1];
                    if(tt) new_fML = MIN2(new_fML, c[indx[j-1]+i+1] + E_MLstem(tt, mm5, mm3) + 2*params.MLbase);
                    break;
            }
            
            
            /* modular decomposition -------------------------------*/
            for (decomp = INF, k = i + 1 + TURN; k <= j - 2 - TURN; k++)
                decomp = MIN2(decomp, Fmi[k]+fML[indx[j]+k+1]);
            DMLi[j] = decomp;               /* store for use in ML decompositon */
            new_fML = MIN2(new_fML,decomp);
            
            /* coaxial stacking */
            if (dangle_model==3) {
                /* additional ML decomposition as two coaxially stacked helices */
                for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++) {
                    type = ptype[indx[k]+i]; type = rtype[type];
                    type_2 = ptype[indx[j]+k+1]; type_2 = rtype[type_2];
                    if (type && type_2)
                        decomp = MIN2(decomp,
                                      c[indx[k]+i]+c[indx[j]+k+1]+params.stack[type][type_2]);
                }
                
                decomp += 2*params.MLintern[1];        /* no TermAU penalty if coax stack */
                new_fML = MIN2(new_fML, decomp);
            }
            fML[ij] = Fmi[j] = new_fML;     /* substring energy */
            
            Ints FF;
            FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
            FF = cc1; cc1=cc; cc=FF;
            for (int k=1; k<=length; k++) { cc[k]=Fmi[k]=DMLi[k]=INF; }

        }
        
    }
    
    
    return 0;
}












