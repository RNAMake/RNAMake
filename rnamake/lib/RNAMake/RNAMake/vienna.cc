//
//  vienna_clone.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <float.h>
#include "vienna.h"
#include "pair_mat.h"

#define LOCALITY          0.      /* locality parameter for base-pairs */

#define SCALE 10
/**
 *** dangling ends should never be destabilizing, i.e. expdangle>=1<BR>
 *** specific heat needs smooth function (2nd derivative)<BR>
 *** we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,
 *** in the interval b<x<sqrt(3)/2
 */
#define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):\
SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))


float
Vienna::fold(
    String const & string) {
    
    actual_size_ = (int)string.length();
    
    for(int i = 0; i < S1.size(); i++) {
        S1[i] = 0;
    }

    encode_sequence(string, S, 0);
    encode_sequence(string, S1, 1);
    make_ptypes(S, structure);
    
    int energy = fill_arrays(string);
    
    backtrack(string, 0);
    parenthesis_structure((int)string.length());
    free_energy_ = (float) energy/100.;

    return free_energy_;
}

plists const &
Vienna::bp_probabilities(
    String const & sequence) {
    
    get_iindx(my_iindx, (int)sequence.size());
    get_iindx(iindx, (int)sequence.size());
    
    float min_en = fold(sequence);
    int n = (int)sequence.length();

    double temperature = 37;
    double betaScale   = 1.;
    double sfact       = 1.07;
    double kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */
    double bppmThreshold = 1e-5;

    
    pf.pf_scale = exp(-(sfact*min_en)/kT/n);
    scale_pf_params(n);
    
    encode_sequence(sequence, S, 0);
    encode_sequence(sequence, S1, 1);
    make_ptypes_2(S, structure);
    pf_linear(sequence);
    
    //reset basepair probabilities
    for(int i =0; i < n*n; i++) {
        pl[i].p = 0; pl[i].i = 0; pl[i].j = 0;
    }
    
    double Q = q[my_iindx[1]-n];
    double free_energy = (-log(Q)-n*log(pf.pf_scale))*pf.kT/1000.0;
    pf_create_bppm(sequence);
    assign_plist_from_pr(n, bppmThreshold);
    
    return pl;
    
}

void
Vienna::assign_plist_from_pr(
    int length,
    double cut_off) {
    
    int i, j, n, count;
    count = 0;
    n     = 2;
    
    get_iindx(diindx,length);
    
    for (i=1; i<length; i++) {
        for (j=i+1; j<=length; j++) {
            if(probs[diindx[i] - j] < cut_off) { continue; }
            pl[count].i = i;
            pl[count].j = j;
            pl[count].p  = sqrt(probs[diindx[i] - j]); //added sqrt to match dot.ps
            pl[count].type = 0;
            count++;
        }
    }
}


void
Vienna::pf_create_bppm(
    String const & sequence) {
    
    int n, i,j,k,l, ij, kl, ii, i1, ll, type, type_2, tt, u1, ov=0, qo=0.;
    double  temp, Qmax=0, prm_MLb;
    double  prmt,prmt1;
    Floats tmp, G;
    double  tmp2;
    double  expMLclosing = pf.expMLclosing;
    double      max_real;
    double  expMLstem = 0;
    max_real = DBL_MAX;
    const char * s = sequence.c_str();
    int circular = 0;
    int with_gquad = 0;
    int VRNA_GQUAD_MIN_BOX_SIZE =0;
    
    if((S.size() != 0) && (S1.size() != 0)){
        n = S[0];
        Qmax=0;
        
        for (k=0; k<=n*n; k++) {
            probs[k] = 0;
        }
        
        for (k=1; k<=n; k++) {
            q1k[k] = q[my_iindx[1] - k];
            qln[k] = q[my_iindx[k] -n];
            prml[k] = 0;
            prm_l[k] = prm_l1[k] = 0;
        }
        q1k[0] = 1.0;
        qln[n+1] = 1.0;
        
        /*  pr = q; */     /* recycling */
        
        
        /* 1. exterior pair i,j and initialization of pr array */
        if(circular){
            for (i=1; i<=n; i++) {
                for (j=i; j<=MIN2(i+TURN,n); j++)
                    probs[my_iindx[i]-j] = 0;
                for (j=i+TURN+1; j<=n; j++) {
                    ij = my_iindx[i]-j;
                    type = ptype[ij];
                    if (type&&(qb[ij]>0.)) {
                        probs[ij] = 1./qo;
                        int rt = rtype[type];
                        
                        /* 1.1. Exterior Hairpin Contribution */
                        int u = i + n - j -1;
                        /* get the loop sequence */
                        char loopseq[10];
                        if (u<7){
                            strcpy(loopseq , s+j-1);
                            strncat(loopseq, s, i);
                        }
                        tmp2 = exp_E_Hairpin(u, rt, S1[j+1], S1[i-1], loopseq) * scale[u];
                        
                        /* 1.2. Exterior Interior Loop Contribution                    */
                        /* 1.2.1. i,j  delimtis the "left" part of the interior loop    */
                        /* (j,i) is "outer pair"                                                */
                        for(k=1; k < i-TURN-1; k++){
                            int ln1, lstart;
                            ln1 = k + n - j - 1;
                            if(ln1>MAXLOOP) break;
                            lstart = ln1+i-1-MAXLOOP;
                            if(lstart<k+TURN+1) lstart = k + TURN + 1;
                            for(l=lstart; l < i; l++){
                                int ln2, type_2;
                                type_2 = ptype[my_iindx[k]-l];
                                if (type_2==0) continue;
                                ln2 = i - l - 1;
                                if(ln1+ln2>MAXLOOP) continue;
                                tmp2 += qb[my_iindx[k] - l]
                                * exp_E_IntLoop(ln1,
                                                ln2,
                                                rt,
                                                rtype[type_2],
                                                S1[j+1],
                                                S1[i-1],
                                                S1[k-1],
                                                S1[l+1])
                                * scale[ln1 + ln2];
                            }
                        }
                        /* 1.2.2. i,j  delimtis the "right" part of the interior loop  */
                        for(k=j+1; k < n-TURN; k++){
                            int ln1, lstart;
                            ln1 = k - j - 1;
                            if((ln1 + i - 1)>MAXLOOP) break;
                            lstart = ln1+i-1+n-MAXLOOP;
                            if(lstart<k+TURN+1) lstart = k + TURN + 1;
                            for(l=lstart; l <= n; l++){
                                int ln2, type_2;
                                type_2 = ptype[my_iindx[k]-l];
                                if (type_2==0) continue;
                                ln2 = i - 1 + n - l;
                                if(ln1+ln2>MAXLOOP) continue;
                                tmp2 += qb[my_iindx[k] - l]
                                * exp_E_IntLoop(ln2,
                                                ln1,
                                                rtype[type_2],
                                                rt,
                                                S1[l+1],
                                                S1[k-1],
                                                S1[i-1],
                                                S1[j+1])
                                * scale[ln1 + ln2];
                            }
                        }
                        /* 1.3 Exterior multiloop decomposition */
                        /* 1.3.1 Middle part                    */
                        if((i>TURN+2) && (j<n-TURN-1))
                            tmp2 += qm[my_iindx[1]-i+1]
                            * qm[my_iindx[j+1]-n]
                            * expMLclosing
                            * exp_E_MLstem(type, S1[i-1], S1[j+1]);
                        
                        /* 1.3.2 Left part                      */
                        for(k=TURN+2; k < i-TURN-2; k++)
                            tmp2 += qm[my_iindx[1]-k]
                            * qm1[jindx[i-1]+k+1]
                            * expMLbase[n-j]
                            * expMLclosing
                            * exp_E_MLstem(type, S1[i-1], S1[j+1]);
                        
                        /* 1.3.3 Right part                      */
                        for(k=j+TURN+2; k < n-TURN-1;k++)
                            tmp2 += qm[my_iindx[j+1]-k]
                            * qm1[jindx[n]+k+1]
                            * expMLbase[i-1]
                            * expMLclosing
                            * exp_E_MLstem(type, S1[i-1], S1[j+1]);
                        
                        /* all exterior loop decompositions for pair i,j done  */
                        probs[ij] *= tmp2;
                        
                    }
                    else probs[ij] = 0;
                }
            }
        } /* end if(circular)  */
        else {
            for (i=1; i<=n; i++) {
                for (j=i; j<=MIN2(i+TURN,n); j++) probs[my_iindx[i]-j] = 0;
                for (j=i+TURN+1; j<=n; j++) {
                    ij = my_iindx[i]-j;
                    type = ptype[ij];
                    if (type&&(qb[ij]>0.)) {
                        probs[ij] = q1k[i-1]*qln[j+1]/q1k[n];
                        probs[ij] *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1);
                    }
                    else
                        probs[ij] = 0.;
                }
            }
        } /* end if(!circular)  */
        
        for (l=n; l>TURN+1; l--) {
            
            /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
            for (k=1; k<l-TURN; k++) {
                kl = my_iindx[k]-l;
                type_2 = ptype[kl];
                if (type_2==0) continue;
                type_2 = rtype[type_2];
                if (qb[kl]==0.) continue;
                
                tmp2 = 0.;
                for (i=MAX2(1,k-MAXLOOP-1); i<=k-1; i++)
                    for (j=l+1; j<=MIN2(l+ MAXLOOP -k+i+2,n); j++) {
                        ij = my_iindx[i] - j;
                        type = ptype[ij];
                        if (type && (probs[ij]>0.)) {
                            /* add *scale[u1+u2+2] */
                            tmp2 +=  probs[ij]
                            * (scale[k-i+j-l]
                               * exp_E_IntLoop(k - i - 1,
                                               j - l - 1,
                                               type,
                                               type_2,
                                               S1[i + 1],
                                               S1[j - 1],
                                               S1[k - 1],
                                               S1[l + 1]));
                        }
                    }
                probs[kl] += tmp2;
            }
            
            if(with_gquad){
                /* 2.5. bonding k,l as gquad enclosed by i,j */
                Floats expintern;
                double qe;
                
                if(l < n - 3){
                    for(k = 2; k <= l - VRNA_GQUAD_MIN_BOX_SIZE; k++){
                        kl = my_iindx[k]-l;
                        if (G[kl]==0.) continue;
                        tmp2 = 0.;
                        i = k - 1;
                        for(j = MIN2(l + MAXLOOP + 1, n); j > l + 3; j--){
                            ij = my_iindx[i] - j;
                            type = ptype[ij];
                            if(!type) continue;
                            qe = (type > 2) ? pf.expTermAU : 1.;
                            tmp2 += probs[ij] * qe * expintern[j-l-1] * pf.expmismatchI[type][S1[i+1]][S1[j-1]] * scale[2];
                        }
                        probs[kl] += tmp2 * G[kl];
                    }
                }
                
                if (l < n - 1){
                    for (k=3; k<=l-VRNA_GQUAD_MIN_BOX_SIZE; k++) {
                        kl = my_iindx[k]-l;
                        if (G[kl]==0.) continue;
                        tmp2 = 0.;
                        for (i=MAX2(1,k-MAXLOOP-1); i<=k-2; i++){
                            u1 = k - i - 1;
                            for (j=l+2; j<=MIN2(l + MAXLOOP - u1 + 1,n); j++) {
                                ij = my_iindx[i] - j;
                                type = ptype[ij];
                                if(!type) continue;
                                qe = (type > 2) ? pf.expTermAU : 1.;
                                tmp2 += probs[ij] * qe * expintern[u1+j-l-1] * pf.expmismatchI[type][S1[i+1]][S1[j-1]] * scale[2];
                            }
                        }
                        probs[kl] += tmp2 * G[kl];
                    }
                }
                
                if(l < n){
                    for(k = 4; k <= l - VRNA_GQUAD_MIN_BOX_SIZE; k++){
                        kl = my_iindx[k]-l;
                        if (G[kl]==0.) continue;
                        tmp2 = 0.;
                        j = l + 1;
                        for (i=MAX2(1,k-MAXLOOP-1); i < k - 3; i++){
                            ij = my_iindx[i] - j;
                            type = ptype[ij];
                            if(!type) continue;
                            qe = (type > 2) ? pf.expTermAU : 1.;
                            tmp2 += probs[ij] * qe * expintern[k - i - 1] * pf.expmismatchI[type][S1[i+1]][S1[j-1]] * scale[2];
                        }
                        probs[kl] += tmp2 * G[kl];
                    }
                }
            }
            
            /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
            prm_MLb = 0.;
            if (l<n) for (k=2; k<l-TURN; k++) {
                i = k-1;
                prmt = prmt1 = 0.0;
                
                ii = my_iindx[i];     /* ii-j=[i,j]     */
                ll = my_iindx[l+1];   /* ll-j=[l+1,j-1] */
                tt = ptype[ii-(l+1)]; tt=rtype[tt];
                /* (i, l+1) closes the ML with substem (k,l) */
                if(tt)
                    prmt1 = probs[ii-(l+1)] * expMLclosing * exp_E_MLstem(tt, S1[l], S1[i+1]);
                
                /* (i,j) with j>l+1 closes the ML with substem (k,l) */
                for (j=l+2; j<=n; j++) {
                    tt = ptype[ii-j]; tt = rtype[tt];
                    if(tt)
                        prmt += probs[ii-j] * exp_E_MLstem(tt, S1[j-1], S1[i+1]) * qm[ll-(j-1)];
                }
                kl = my_iindx[k]-l;
                tt = ptype[kl];
                prmt *= expMLclosing;
                prml[ i] = prmt;
                prm_l[i] = prm_l1[i]*expMLbase[1]+prmt1;
                
                prm_MLb = prm_MLb*expMLbase[1] + prml[i];
                /* same as:    prm_MLb = 0;
                 for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
                
                prml[i] = prml[ i] + prm_l[i];
                
                if(with_gquad){
                    if ((!tt) && (G[kl] == 0.)) continue;
                } else {
                    if (qb[kl] == 0.) continue;
                }
                
                temp = prm_MLb;
                
                for (i=1;i<=k-2; i++)
                    temp += prml[i]*qm[my_iindx[i+1] - (k-1)];
                
                if(with_gquad){
                    if(tt)
                        temp    *= exp_E_MLstem(tt, (k>1) ? S1[k-1] : -1, (l<n) ? S1[l+1] : -1) * scale[2];
                    else
                        temp    *= G[kl] * expMLstem * scale[2];
                } else {
                    temp    *= exp_E_MLstem(tt, (k>1) ? S1[k-1] : -1, (l<n) ? S1[l+1] : -1) * scale[2];
                }
                
                probs[kl]  += temp;
                
                if (probs[kl]>Qmax) {
                    Qmax = probs[kl];
                    if (Qmax>max_real/10.)
                        fprintf(stderr, "P close to overflow: %d %d %g %g\n",
                                i, j, probs[kl], qb[kl]);
                }
                if (probs[kl]>=max_real) {
                    ov++;
                    probs[kl]=FLT_MAX;
                }
                
            } /* end for (k=..) */
            tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;
            
        }  /* end for (l=..)   */
        
        for (i=1; i<=n; i++)
            for (j=i+TURN+1; j<=n; j++) {
                ij = my_iindx[i]-j;
                
                if(with_gquad){
                    if (qb[ij] > 0.)
                        probs[ij] *= qb[ij];
                    if (G[ij] > 0.){
                        probs[ij] += q1k[i-1] * G[ij] * qln[j+1]/q1k[n];
                    }
                } else {
                    if (qb[ij] > 0.)
                        probs[ij] *= qb[ij];
                }
            }
        
        if (pf_structure.size() == 0)
            bppm_to_structure(structure, probs, n);
        if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
                          "you might try a smaller pf_scale than %g\n",
                          ov, pf.pf_scale);
    } /* end if((S != NULL) && (S1 != NULL))  */
    else
        printf("bppm calculations have to be done after calling forward recursion\n");
    return;

}


char
Vienna::bppm_symbol(
    const float *x) {
    if( x[0] > 0.667 )  return '.';
    if( x[1] > 0.667 )  return '(';
    if( x[2] > 0.667 )  return ')';
    if( (x[1]+x[2]) > x[0] ) {
        if( (x[1]/(x[1]+x[2])) > 0.667) return '{';
        if( (x[2]/(x[1]+x[2])) > 0.667) return '}';
        else return '|';
    }
    if( x[0] > (x[1]+x[2]) ) return ',';
    return ':';
}


void
Vienna::bppm_to_structure(
    String & structure,
    Floats & p,
    unsigned int length) {
    
    int    i, j;
    get_iindx(diindx,length);
    float  P[3];   /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */
    
    for( j=1; j<=length; j++ ) {
        P[0] = 1.0;
        P[1] = P[2] = 0.0;
        for( i=1; i<j; i++) {
            P[2] += p[diindx[i]-j];    /* j is paired downstream */
            P[0] -= p[diindx[i]-j];    /* j is unpaired */
        }
        for( i=j+1; i<=length; i++ ) {
            P[1] += p[diindx[j]-i];    /* j is paired upstream */
            P[0] -= p[diindx[j]-i];    /* j is unpaired */
        }
        structure[j-1] = bppm_symbol(P);
    }
    structure[length] = '\0';

}


void
Vienna::get_boltzmann_factors(
    float temp,
    float betaScale,
    model_detailsT const & md,
    float pf_scale) {
    
    unsigned  int i, j, k, l;
    double        kT, TT;
    double        GT;
    
    pf.model_details = md;
    pf.temperature   = temp;
    pf.alpha         = betaScale;
    pf.kT = kT       = betaScale*(temp+K0)*GASCONST;   /* kT in cal/mol  */
    pf.pf_scale      = pf_scale;
    TT                = (temp+K0)/(Tmeasure);
    
    for (i=0; i<31; i++){
        GT  = hairpindH[i] - (hairpindH[i] - hairpin37[i])*TT;
        pf.exphairpin[i] = exp( -GT*10./kT);
    }
    
    for (i=0; i<=MIN2(30, MAXLOOP); i++) {
        GT =  bulgedH[i]- (bulgedH[i] - bulge37[i])*TT;
        pf.expbulge[i] = exp( -GT*10./kT);
        GT =  internal_loopdH[i] - (internal_loopdH[i] - internal_loop37[i])*TT;
        pf.expinternal[i] = exp( -GT*10./kT);
    }
    
    /* special case of size 2 interior loops (single mismatch) */
    pf.expinternal[2] = exp ( -80*10./kT);
    pf.lxc = lxc37*TT;
    
    GT =  DuplexInitdH - (DuplexInitdH - DuplexInit37)*TT;
    pf.expDuplexInit = exp( -GT*10./kT);
    
    for (i=31; i<=MAXLOOP; i++) {
        GT = bulge37[30]*TT + (pf.lxc*log( i/30.));
        pf.expbulge[i] = exp( -GT*10./kT);
        GT = internal_loop37[30]*TT + (pf.lxc*log( i/30.));
        pf.expinternal[i] = exp( -GT*10./kT);
    }
    
    GT = niniodH - (niniodH - ninio37)*TT;
    for (j=0; j<=MAXLOOP; j++)
        pf.expninio[2][j]=exp(-MIN2(MAX_NINIO,j*GT)*10./kT);
    
    for (i=0; (i*7)<strlen(Tetraloops); i++) {
        GT = TetraloopdH[i] - (TetraloopdH[i]-Tetraloop37[i])*TT;
        pf.exptetra[i] = exp( -GT*10./kT);
    }
    for (i=0; (i*5)<strlen(Triloops); i++) {
        GT = TriloopdH[i] - (TriloopdH[i]-Triloop37[i])*TT;
        pf.exptri[i] = exp( -GT*10./kT);
    }
    for (i=0; (i*9)<strlen(Hexaloops); i++) {
        GT = HexaloopdH[i] - (HexaloopdH[i]-Hexaloop37[i])*TT;
        pf.exphex[i] = exp( -GT*10./kT);
    }
    GT =  ML_closingdH - (ML_closingdH - ML_closing37)*TT;
    pf.expMLclosing = exp( -GT*10./kT);
    
    for (i=0; i<=NBPAIRS; i++) {
        GT =  ML_interndH - (ML_interndH - ML_intern37)*TT;
        /* if (i>2) GT += TerminalAU; */
        pf.expMLintern[i] = exp( -GT*10./kT);
    }

    GT = TerminalAUdH - (TerminalAUdH - TerminalAU37)*TT;
    pf.expTermAU = exp(-GT*10./kT);
    
    GT = ML_BASEdH - (ML_BASEdH - ML_BASE37)*TT;
    
    pf.expMLbase=exp(-10.*GT/kT);
    
    /* if dangles==0 just set their energy to 0,
     don't let dangle energies become > 0 (at large temps),
     but make sure go smoothly to 0                        */
    for (i=0; i<=NBPAIRS; i++) {
        for (j=0; j<=4; j++) {
            if (md.dangles) {
                GT = dangle5_dH[i][j] - (dangle5_dH[i][j] - dangle5_37[i][j])*TT;
                pf.expdangle5[i][j] = exp(SMOOTH(-GT)*10./kT);
                GT = dangle3_dH[i][j] - (dangle3_dH[i][j] - dangle3_37[i][j])*TT;
                pf.expdangle3[i][j] =  exp(SMOOTH(-GT)*10./kT);
            } else
                pf.expdangle3[i][j] = pf.expdangle5[i][j] = 1;
        }
    }
    
    /* stacking energies */
    for (i=0; i<=NBPAIRS; i++) {
        for (j=0; j<=NBPAIRS; j++) {
            GT =  stackdH[i][j] - (stackdH[i][j] - stack37[i][j])*TT;
            pf.expstack[i][j] = exp( -GT*10./kT);
        }
    }
        
    /* mismatch energies */
    for (i=0; i<=NBPAIRS; i++) {
        for (j=0; j<5; j++) {
            for (k=0; k<5; k++) {
                GT =  mismatchIdH[i][j][k] - ( mismatchIdH[i][j][k] - mismatchI37[i][j][k])*TT;
                pf.expmismatchI[i][j][k] = exp(-GT*10.0/kT);
                GT = mismatch1nIdH[i][j][k] - (mismatch1nIdH[i][j][k] - mismatch1nI37[i][j][k])*TT;
                pf.expmismatch1nI[i][j][k] = exp(-GT*10.0/kT);
                GT = mismatchHdH[i][j][k] - (mismatchHdH[i][j][k] - mismatchH37[i][j][k])*TT;
                pf.expmismatchH[i][j][k] = exp(-GT*10.0/kT);
                if (md.dangles) {
                    GT = mismatchMdH[i][j][k] - (mismatchMdH[i][j][k] - mismatchM37[i][j][k])*TT;
                    pf.expmismatchM[i][j][k] = exp(SMOOTH(-GT)*10.0/kT);
                    GT = mismatchExtdH[i][j][k] - (mismatchExtdH[i][j][k] - mismatchExt37[i][j][k])*TT;
                    pf.expmismatchExt[i][j][k] = exp(SMOOTH(-GT)*10.0/kT);
                }
                else{
                    pf.expmismatchM[i][j][k] = pf.expmismatchExt[i][j][k] = 1.;
                }
                GT = mismatch23IdH[i][j][k] - (mismatch23IdH[i][j][k] - mismatch23I37[i][j][k])*TT;
                pf.expmismatch23I[i][j][k] = exp(-GT*10.0/kT);
            }
        }
    }

    /* interior lops of length 2 */
    for (i=0; i<=NBPAIRS; i++) {
        for (j=0; j<=NBPAIRS; j++) {
            for (k=0; k<5; k++) {
                for (l=0; l<5; l++) {
                    GT = int11_dH[i][j][k][l] -
                    (int11_dH[i][j][k][l] - int11_37[i][j][k][l])*TT;
                    pf.expint11[i][j][k][l] = exp(-GT*10./kT);
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
                    for (m=0; m<5; m++) {
                        GT = int21_dH[i][j][k][l][m] -
                        (int21_dH[i][j][k][l][m] - int21_37[i][j][k][l][m])*TT;
                        pf.expint21[i][j][k][l][m] = exp(-GT*10./kT);
                    }
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
                            GT = int22_dH[i][j][k][l][m][n] -
                            (int22_dH[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*TT;
                            pf.expint22[i][j][k][l][m][n] = exp(-GT*10./kT);
                        }
                    }
                }
            }
        }
    }

    strncpy(pf.Tetraloops, Tetraloops, 281);
    strncpy(pf.Triloops, Triloops, 241);
    strncpy(pf.Hexaloops, Hexaloops, 361);


}


void
Vienna::scale_pf_params(
    int length) {
    
    unsigned int  i;
    double        scaling_factor;
    
    scaling_factor = pf.pf_scale;
    
    /* scaling factors (to avoid overflows) */
    if (scaling_factor == -1) { /* mean energy for random sequences: 184.3*length cal */
        scaling_factor = exp(-(-185+(pf.temperature-37.)*7.27)/pf.kT);
        if (scaling_factor<1) scaling_factor=1;
        pf.pf_scale = scaling_factor;
    }
    scale[0] = 1.;
    scale[1] = 1./scaling_factor;
    expMLbase[0] = 1;
    expMLbase[1] = pf.expMLbase/scaling_factor;
    for (i=2; i<=length; i++) {
        scale[i] = scale[i/2]*scale[i-(i/2)];
        expMLbase[i] = pow(pf.expMLbase, (double)i) * scale[i];
    }

    
}

void
Vienna::pf_linear(
    String const & sequence) {
    
    const char * s = sequence.c_str();
    int n, i,j,k,l, ij, u,u1,d,ii, type, type_2, tt, minl, maxl;
    int noGUclosure;
    int circular = 0;
    float expMLstem = 0.;
    
    float temp, Qmax=0;
    float qbt1;
    Floats tmp;
    
    float  expMLclosing = pf.expMLclosing;
    float  max_real = FLT_MAX;
    n = (int)sequence.length();
    
    noGUclosure = pf.model_details.noGUclosure;

    for (d=0; d<=TURN; d++) {
        for (i=1; i<=n-d; i++) {
            j=i+d;
            ij = my_iindx[i]-j;
            q[ij]=1.0*scale[d+1];
            qb[ij]=qm[ij]=0.0;
        }
    }
    
    for (i=1; i<=n; i++) {
        qq[i]=qq1[i]=qqm[i]=qqm1[i]=0;
    }
    
    for (j=TURN+2;j<=n; j++) {
        for (i=j-TURN-1; i>=1; i--) {
            /* construction of partition function of segment i,j*/
            /*firstly that given i binds j : qb(i,j) */
            u = j-i-1; ij = my_iindx[i]-j;
            type = ptype[ij];
            if (type!=0) {
                /*hairpin contribution*/
                if (((type==3)||(type==4))&&noGUclosure) qbt1 = 0;
                else {
                    qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], s+i-1) * scale[u+2];
                }
                /* interior loops with interior pair k,l */
                for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
                    u1 = k-i-1;
                    for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
                        type_2 = ptype[my_iindx[k]-l];
                        if (type_2) {
                            type_2 = rtype[type_2];
                            qbt1 += qb[my_iindx[k]-l] * (scale[u1+j-l+1] *
                                                         exp_E_IntLoop(u1, j-l-1, type, type_2,
                                                                       S1[i+1], S1[j-1], S1[k-1], S1[l+1]));
                        }
                    }
                }
                /*multiple stem loop contribution*/
                ii = my_iindx[i+1]; /* ii-k=[i+1,k-1] */
                temp = 0.0;
                for (k=i+2; k<=j-1; k++) temp += qm[ii-(k-1)]*qqm1[k];
                tt = rtype[type];
                qbt1 += temp * expMLclosing * exp_E_MLstem(tt, S1[j-1], S1[i+1]) * scale[2];
                qb[ij] = qbt1;
            }
            /* end if (type!=0) */
            else
                qb[ij] = 0.0;
            
            /* construction of qqm matrix containing final stem
             contributions to multiple loop partition function
             from segment i,j */
            qqm[i] = qqm1[i]*expMLbase[1];
            if (type) {
                qbt1 = qb[ij] * exp_E_MLstem(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1);
                qqm[i] += qbt1;
            }
            
                    
            /*construction of qm matrix containing multiple loop
             partition function contributions from segment i,j */
            temp = 0.0;
            ii = my_iindx[i];  /* ii-k=[i,k-1] */
            for (k=j; k>i; k--) temp += (qm[ii-(k-1)] + expMLbase[k-i])*qqm[k];
            qm[ij] = (temp + qqm[i]);
            
            /*auxiliary matrix qq for cubic order q calculation below */
            qbt1=0.0;
            if (type){
                qbt1 += qb[ij];
                qbt1 *= exp_E_ExtLoop(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1);
            }
      
            qq[i] = qq1[i]*scale[1] + qbt1;
            
            /*construction of partition function for segment i,j */
            temp = 1.0*scale[1+j-i] + qq[i];
            for (k=i; k<=j-1; k++) temp += q[ii-k]*qq[k+1];
            q[ij] = temp;
            if (temp>Qmax) {
                Qmax = temp;
                if (Qmax>max_real/10.) {
                    fprintf(stderr, "Q close to overflow: %d %d %g\n", i,j,temp);
                    exit(1);
                }
            }
            if (temp>=max_real) {
                fprintf(stderr, "overflow in pf_fold while calculating q[%d,%d]\n"
                        "use larger pf_scale", i,j);
                exit(1);
            }
        }
        tmp = qq1;  qq1 =qq;  qq =tmp;
        tmp = qqm1; qqm1=qqm; qqm=tmp;
    }    
}


void
Vienna::backtrack(
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
Vienna::make_ptypes(
    Shorts const & S,
    String const & structure) {
    
    int n,i,j,k,l, noLP;
    noLP = params.model_details.noLP;
    
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

void
Vienna::make_ptypes_2(
    Shorts const & S,
    String const & structure) {
    
    int n,i,j,k,l, noLP;
    noLP = params.model_details.noLP;
    n = S[0];
    
    for (k=1; k<n-TURN; k++) {
        for (l=1; l<=2; l++) {
            int type,ntype=0,otype=0;
            i=k; j = i+TURN+l; if (j>n) continue;
            type = BP_pair[S[i]][S[j]];
            while ((i>=1)&&(j<=n)) {
                if ((i>1)&&(j<n)) ntype = BP_pair[S[i-1]][S[j+1]];
                if (noLP && (!otype) && (!ntype))
                    type = 0; /* i.j can only form isolated pairs */
                qb[my_iindx[i]-j] = 0.;
                ptype[my_iindx[i]-j] = (char) type;
                otype =  type;
                type  = ntype;
                i--; j++;
            }
        }
    }
}


int
Vienna::fill_arrays(
    String const & string)  {
    
    const char * s = string.c_str();

    int   i, j, k, length, energy, en, mm5, mm3;
    int   decomp, new_fML;
    int   no_close, type, type_2, tt, max_separation;
    int   bonus=0;
    int   dangle_model, noGUclosure;
    int   noLonelyPairs=0;
    int   circular=0;
    get_indx(indx);

    dangle_model  = params.model_details.dangles;
    noGUclosure   = params.model_details.noGUclosure;
    length = (int)string.length();
    
    for(j=1;j < S1.size(); j++) {
        if(S1[j] < 0) { S1[j] = 0; }
    }
    
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
    
    s = NULL;
    
    return f5[length];

}


void
Vienna::parenthesis_structure(
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











