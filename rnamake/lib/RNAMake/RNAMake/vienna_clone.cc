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
    
    actual_size_ = string.length();
    
    //update length of data arrays if size is greater then previous
    //get_arrays((int)string.length());
    encode_sequence(string, S, 0);
    encode_sequence(string, S1, 1);
    make_ptypes(S, structure);
    
    //sector = sects(MAXSECTORS);
    //base_pair2 = bondTs(4*(1+size_/2));
    
    int energy = fill_arrays(string);
    
    backtrack(string, 0);
    parenthesis_structure((int)string.length());
    
    return (float) energy/100.;
}

void
ViennaClone::dotplot(
    String const & string) {
    
    
}

void
ViennaClone::get_boltzmann_factors(
    float temperature,
    float betaScale,
    model_detailsT const & md,
    float pf_scale) {
    
}


void
ViennaClone::backtrack(
    String const & string,
    int s) {
    
    int   i, j, ij, k, l1, mm5, mm3, length, energy, en, nnew;
    int   no_close, type, type_2, tt, minq, maxq, c0, c1, c2, c3;
    int   bonus;
    int   b=0;
    int   dangle_model = params.model_details.dangles;
    int   no_closingGU   = params.model_details.noGUclosure;
    int   noLonelyPairs=0;

    const char * str = string.c_str();
    
    length = (int)string.length();
    
    if(s == 0) {
        sector[++s].i = 1;
        sector[s].j = length;
        sector[s].ml = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')? 2: 0);
    }
    
    while (s>0) {
        int ml, fij, fi, cij, traced, i1, j1, p, q, jj=0, gq=0;
        int canonical = 1;     /* (i,j) closes a canonical structure */
        i  = sector[s].i;
        j  = sector[s].j;
        ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                                occur in the fML- (1) or in the f-array (0) */
        if (ml==2) {
            base_pair2[++b].i = i;
            base_pair2[b].j   = j;
            goto repeat1;
        }
        
        else if(ml==7) { /* indicates that i,j are enclosing a gquadruplex */
            /* actually, do something here */
        }
        
        if (j < i+TURN+1) continue; /* no more pairs in this interval */
        
        fij = (ml == 1)? fML[indx[j]+i] : f5[j];
        fi  = (ml == 1)?(fML[indx[j-1]+i]+params.MLbase): f5[j-1];
        
        if (fij == fi) {  /* 3' end is unpaired */
            sector[++s].i = i;
            sector[s].j   = j-1;
            sector[s].ml  = ml;
            continue;
        }
        
        if (ml == 0) { /* backtrack in f5 */
            switch(dangle_model){
                case 0:   /* j is paired. Find pairing partner */
                    for(k=j-TURN-1,traced=0; k>=1; k--){
                        
                        type = ptype[indx[j]+k];
                        if(type)
                            if(fij == E_ExtLoop(type, -1, -1) + c[indx[j]+k] + f5[k-1]){
                                traced=j; jj = k-1;
                                break;
                            }
                    }
                    break;
                    
                case 2:   mm3 = (j<length) ? S1[j+1] : -1;
                    for(k=j-TURN-1,traced=0; k>=1; k--){
                        
                        type = ptype[indx[j]+k];
                        if(type)
                            if(fij == E_ExtLoop(type, (k>1) ? S1[k-1] : -1, mm3) + c[indx[j]+k] + f5[k-1]){
                                traced=j; jj = k-1;
                                break;
                            }
                    }
                    break;
                    
                default:  for(traced = 0, k=j-TURN-1; k>1; k--){
                 
                    type = ptype[indx[j] + k];
                    if(type){
                        en = c[indx[j] + k];
                        if(fij == f5[k-1] + en + E_ExtLoop(type, -1, -1)){
                            traced = j;
                            jj = k-1;
                            break;
                        }
                        if(fij == f5[k-2] + en + E_ExtLoop(type, S1[k-1], -1)){
                            traced = j;
                            jj = k-2;
                            break;
                        }
                    }
                    type = ptype[indx[j-1] + k];
                    if(type){
                        en = c[indx[j-1] + k];
                        if(fij == f5[k-1] + en + E_ExtLoop(type, -1, S1[j])){
                            traced = j-1;
                            jj = k-1;
                            break;
                        }
                        if(fij == f5[k-2] + en + E_ExtLoop(type, S1[k-1], S1[j])){
                            traced = j-1;
                            jj = k-2;
                            break;
                        }
                    }
                }
                    if(!traced){
                        
                        type = ptype[indx[j]+1];
                        if(type){
                            if(fij == c[indx[j]+1] + E_ExtLoop(type, -1, -1)){
                                traced = j;
                                jj = 0;
                                break;
                            }
                        }
                        type = ptype[indx[j-1]+1];
                        if(type){
                            if(fij == c[indx[j-1]+1] + E_ExtLoop(type, -1, S1[j])){
                                traced = j-1;
                                jj = 0;
                                break;
                            }
                        }
                    }
                    break;
            }
            
            if (!traced){
                std::cout << "backtrack failed in f5" << std::endl;
                exit(1);
            }
            /* push back the remaining f5 portion */
            sector[++s].i = 1;
            sector[s].j   = jj;
            sector[s].ml  = ml;
            
            /* trace back the base pair found */
            i=k; j=traced;
            
            base_pair2[++b].i = i;
            base_pair2[b].j   = j;
            goto repeat1;
        }
        else { /* trace back in fML array */
            if (fML[indx[j]+i+1]+params.MLbase == fij) { /* 5' end is unpaired */
                sector[++s].i = i+1;
                sector[s].j   = j;
                sector[s].ml  = ml;
                continue;
            }
            
            ij  = indx[j]+i;
            
            tt  = ptype[ij];
            en  = c[ij];
            switch(dangle_model){
                case 0:   if(fij == en + E_MLstem(tt, -1, -1)){
                    base_pair2[++b].i = i;
                    base_pair2[b].j   = j;
                    goto repeat1;
                }
                    break;
                    
                case 2:   if(fij == en + E_MLstem(tt, S1[i-1], S1[j+1])){
                    base_pair2[++b].i = i;
                    base_pair2[b].j   = j;
                    goto repeat1;
                }
                    break;
                    
                default:  if(fij == en + E_MLstem(tt, -1, -1)){
                    base_pair2[++b].i = i;
                    base_pair2[b].j   = j;
                    goto repeat1;
                }
                    tt = ptype[ij+1];
                    if(fij == c[ij+1] + E_MLstem(tt, S1[i], -1) + params.MLbase){
                        base_pair2[++b].i = ++i;
                        base_pair2[b].j   = j;
                        goto repeat1;
                    }
                    tt = ptype[indx[j-1]+i];
                    if(fij == c[indx[j-1]+i] + E_MLstem(tt, -1, S1[j]) + params.MLbase){
                        base_pair2[++b].i = i;
                        base_pair2[b].j   = --j;
                        goto repeat1;
                    }
                    tt = ptype[indx[j-1]+i+1];
                    if(fij == c[indx[j-1]+i+1] + E_MLstem(tt, S1[i], S1[j]) + 2*params.MLbase){
                        base_pair2[++b].i = ++i;
                        base_pair2[b].j   = --j;
                        goto repeat1;
                    }
                    break;
            }
            
            for(k = i + 1 + TURN; k <= j - 2 - TURN; k++)
                if(fij == (fML[indx[k]+i]+fML[indx[j]+k+1]))
                    break;
            
            if ((dangle_model==3)&&(k > j - 2 - TURN)) { /* must be coax stack */
                ml = 2;
                for (k = i+1+TURN; k <= j - 2 - TURN; k++) {
                    type    = rtype[ptype[indx[k]+i]];
                    type_2  = rtype[ptype[indx[j]+k+1]];
                    if (type && type_2)
                        if (fij == c[indx[k]+i]+c[indx[j]+k+1]+params.stack[type][type_2]+
                            2*params.MLintern[1])
                            break;
                }
            }
            sector[++s].i = i;
            sector[s].j   = k;
            sector[s].ml  = ml;
            sector[++s].i = k+1;
            sector[s].j   = j;
            sector[s].ml  = ml;
            
            if (k>j-2-TURN) {
                std::cout << "backtrack failed in fML" << std::endl;
                exit(1);
            }
            continue;
        }
        
    repeat1:
        
        /*----- begin of "repeat:" -----*/
        ij = indx[j]+i;
        if (canonical)  cij = c[ij];
        
        type = ptype[ij];
        
        bonus = 0;
        if (noLonelyPairs)
            if (cij == c[ij]){
                /* (i.j) closes canonical structures, thus
                 (i+1.j-1) must be a pair                */
                type_2 = ptype[indx[j-1]+i+1]; type_2 = rtype[type_2];
                cij -= params.stack[type][type_2] + bonus;
                base_pair2[++b].i = i+1;
                base_pair2[b].j   = j-1;
                i++; j--;
                canonical=0;
                goto repeat1;
            }
        canonical = 1;
        
        
        no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));
        if (no_close) {
            if (cij == FORBIDDEN) continue;
        } else
            if (cij == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1],str+i-1)+bonus)
                continue;
        
        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
            minq = j-i+p-MAXLOOP-2;
            if (minq<p+1+TURN) minq = p+1+TURN;
            for (q = j-1; q >= minq; q--) {
                
                type_2 = ptype[indx[q]+p];
                if (type_2==0) continue;
                type_2 = rtype[type_2];
                if (no_closingGU)
                    if (no_close||(type_2==3)||(type_2==4))
                        if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */
                
                /* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
                energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                                   S1[i+1], S1[j-1], S1[p-1], S1[q+1]);
                
                nnew = energy+c[indx[q]+p]+bonus;
                traced = (cij == nnew);
                if (traced) {
                    base_pair2[++b].i = p;
                    base_pair2[b].j   = q;
                    i = p, j = q;
                    goto repeat1;
                }
            }
        }
        
        /* end of repeat: --------------------------------------------------*/
        
        /* (i.j) must close a multi-loop */
        tt = rtype[type];
        i1 = i+1; j1 = j-1;
        
        sector[s+1].ml  = sector[s+2].ml = 1;
        
        switch(dangle_model){
            case 0:   en = cij - E_MLstem(tt, -1, -1) - params.MLclosing - bonus;
                for(k = i+2+TURN; k < j-2-TURN; k++){
                    if(en == fML[indx[k]+i+1] + fML[indx[j-1]+k+1])
                        break;
                }
                break;
                
            case 2:   en = cij - E_MLstem(tt, S1[j-1], S1[i+1]) - params.MLclosing - bonus;
                for(k = i+2+TURN; k < j-2-TURN; k++){
                    if(en == fML[indx[k]+i+1] + fML[indx[j-1]+k+1])
                        break;
                }
                break;
                
            default:  for(k = i+2+TURN; k < j-2-TURN; k++){
                en = cij - params.MLclosing - bonus;
                if(en == fML[indx[k]+i+1] + fML[indx[j-1]+k+1] + E_MLstem(tt, -1, -1)){
                    break;
                }
                else if(en == fML[indx[k]+i+2] + fML[indx[j-1]+k+1] + E_MLstem(tt, -1, S1[i+1]) + params.MLbase){
                    i1 = i+2;
                    break;
                }
                else if(en == fML[indx[k]+i+1] + fML[indx[j-2]+k+1] + E_MLstem(tt, S1[j-1], -1) + params.MLbase){
                    j1 = j-2;
                    break;
                }
                else if(en == fML[indx[k]+i+2] + fML[indx[j-2]+k+1] + E_MLstem(tt, S1[j-1], S1[i+1]) + 2*params.MLbase){
                    i1 = i+2;
                    j1 = j-2;
                    break;
                }
                /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
                /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
                if(dangle_model == 3){
                    type_2 = rtype[ptype[indx[k]+i+1]];
                    if (type_2) {
                        en = c[indx[k]+i+1]+params.stack[type][type_2]+fML[indx[j-1]+k+1];
                        if (cij == en+2*params.MLintern[1]+params.MLclosing) {
                            ml = 2;
                            sector[s+1].ml  = 2;
                            traced = 1;
                            break;
                        }
                    }
                    type_2 = rtype[ptype[indx[j-1]+k+1]];
                    if (type_2) {
                        en = c[indx[j-1]+k+1]+params.stack[type][type_2]+fML[indx[k]+i+1];
                        if (cij == en+2*params.MLintern[1]+params.MLclosing) {
                            sector[s+2].ml = 2;
                            traced = 1;
                            break;
                        }
                    }
                }
            }
                break;
        }
        
        if (k<=j-3-TURN) { /* found the decomposition */
            sector[++s].i = i1;
            sector[s].j   = k;
            sector[++s].i = k+1;
            sector[s].j   = j1;
        } else {
            std::cout << "backtracking failed in repeat" << std::endl;
            exit(1);
        }
        
        continue; /* this is a workarround to not accidentally proceed in the following block */
        
    repeat_gquad:
        /*
         now we do some fancy stuff to backtrace the stacksize and linker lengths
         of the g-quadruplex that should reside within position i,j
         */
        {
            int l[3], L, a;
            L = -1;
            
            if(L != -1){
                /* fill the G's of the quadruplex into base_pair2 */
                for(a=0;a<L;a++){
                    base_pair2[++b].i = i+a;
                    base_pair2[b].j   = i+a;
                    base_pair2[++b].i = i+L+l[0]+a;
                    base_pair2[b].j   = i+L+l[0]+a;
                    base_pair2[++b].i = i+L+l[0]+L+l[1]+a;
                    base_pair2[b].j   = i+L+l[0]+L+l[1]+a;
                    base_pair2[++b].i = i+L+l[0]+L+l[1]+L+l[2]+a;
                    base_pair2[b].j   = i+L+l[0]+L+l[1]+L+l[2]+a;
                }
                goto repeat_gquad_exit;
            }
            
            std::cout << "backtracking failed in repeat_gquad" << std::endl;
            exit(1);
            
        }
    repeat_gquad_exit:
        asm("nop");
        
    } /* end of infinite while loop */
    
    base_pair2[0].i = b;    /* save the total number of base pairs */
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
    int   dangle_model, noGUclosure;
    int   noLonelyPairs=0;
    int   circular=0;

    dangle_model  = params.model_details.dangles;
    noGUclosure   = params.model_details.noGUclosure;
    length = (int)string.length();
    
    max_separation = (int) ((1.-LOCALITY)*(double)(length-2)); /* not in use */
    
    for (j=1; j<=length; j++) {
        Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
        cc[j]=cc1[j]=0;
        //f5[j]=f53[j]=0;
    }
    
    for (j=1; j<=length; j++) {
        for (i=(j>TURN?(j-TURN):1); i<j; i++) {
            c[indx[j]+i] = fML[indx[j]+i] = Fmi[indx[j]+i] =INF;
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
        
        }
        
        Ints FF;
        FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
        FF = cc1; cc1=cc; cc=FF;
        for (int k=1; k<=length; k++) { cc[k]=Fmi[k]=DMLi[k]=INF; }
        
    }
    
    /* calculate energies of 5' and 3' fragments */
    f5[TURN+1]= 0;
    
    switch(dangle_model){
            /* dont use dangling end and mismatch contributions at all */
        case 0:   for(j=TURN+2; j<=length; j++){
            f5[j] = f5[j-1];
            for (i=j-TURN-1; i>1; i--){
  
                type = ptype[indx[j]+i];
                if(!type) continue;
                en = c[indx[j]+i];
                f5[j] = MIN2(f5[j], f5[i-1] + en + E_ExtLoop(type, -1, -1));
            }
            
            type=ptype[indx[j]+1];
            if(!type) continue;
            en = c[indx[j]+1];
            f5[j] = MIN2(f5[j], en + E_ExtLoop(type, -1, -1));
        }
            break;
            
            /* always use dangles on both sides */
        case 2:   for(j=TURN+2; j<length; j++){
            f5[j] = f5[j-1];
            for (i=j-TURN-1; i>1; i--){
                
                type = ptype[indx[j]+i];
                if(!type) continue;
                en = c[indx[j]+i];
                f5[j] = MIN2(f5[j], f5[i-1] + en + E_ExtLoop(type, S1[i-1], S1[j+1]));
            }
            
            type=ptype[indx[j]+1];
            if(!type) continue;
            en = c[indx[j]+1];
            f5[j] = MIN2(f5[j], en + E_ExtLoop(type, -1, S1[j+1]));
        }
            f5[length] = f5[length-1];
            for (i=length-TURN-1; i>1; i--){
                
                type = ptype[indx[length]+i];
                if(!type) continue;
                en = c[indx[length]+i];
                f5[length] = MIN2(f5[length], f5[i-1] + en + E_ExtLoop(type, S1[i-1], -1));
            }

            type=ptype[indx[length]+1];
            if(!type) break;
            en = c[indx[length]+1];
            f5[length] = MIN2(f5[length], en + E_ExtLoop(type, -1, -1));
            
            
            break;
            
            /* normal dangles, aka dangle_model = 1 || 3 */
        default:  for(j=TURN+2; j<=length; j++){
            f5[j] = f5[j-1];
            for (i=j-TURN-1; i>1; i--){
          
                type = ptype[indx[j]+i];
                if(type){
                    en = c[indx[j]+i];
                    f5[j] = MIN2(f5[j], f5[i-1] + en + E_ExtLoop(type, -1, -1));
                    f5[j] = MIN2(f5[j], f5[i-2] + en + E_ExtLoop(type, S1[i-1], -1));
                }
                type = ptype[indx[j-1]+i];
                if(type){
                    en = c[indx[j-1]+i];
                    f5[j] = MIN2(f5[j], f5[i-1] + en + E_ExtLoop(type, -1, S1[j]));
                    f5[j] = MIN2(f5[j], f5[i-2] + en + E_ExtLoop(type, S1[i-1], S1[j]));
                }
            }
            
            type = ptype[indx[j]+1];
            if(type) f5[j] = MIN2(f5[j], c[indx[j]+1] + E_ExtLoop(type, -1, -1));
            type = ptype[indx[j-1]+1];
            if(type) f5[j] = MIN2(f5[j], c[indx[j-1]+1] + E_ExtLoop(type, -1, S1[j]));
        }
    }
    
    return f5[length];

}


void
ViennaClone::parenthesis_structure(
    int length) {

    int n, k;
    
    for (n = 0; n < length; n++) {
        structure[n] = '.';
    }
    structure[length] = '\0';
    
    for (k = 1; k <= base_pair2[0].i; k++){
        
        //std::cout << "BP i" << base_pair2[k].i << " " << base_pair2[k].j << std::endl;
        if(base_pair2[k].i == base_pair2[k].j){ /* Gquad bonds are marked as bp[i].i == bp[i].j */
            structure[base_pair2[k].i-1] = '+';
        } else { /* the following ones are regular base pairs */
            structure[base_pair2[k].i-1] = '(';
            structure[base_pair2[k].j-1] = ')';
        }
    }
}











