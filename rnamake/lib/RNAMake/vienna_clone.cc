//
//  vienna_clone.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "vienna_clone.h"


float
ViennaClone::fold(
    String const & string) {
    
    //update length of data arrays if size is greater then previous
    get_arrays((int)string.length());
    encode_sequence(string, S, 0);
    encode_sequence(string, S1, 1);
    make_ptypes(S, structure);
    
    

    
    return 1;
}

void
ViennaClone::make_ptypes(
    Shorts const & S,
    String const & structure) {
    
    int n,i,j,k,l, noLP;
    noLP = params.model_details.noLP;
    n = S[0];
    
    for (k=1; k<n-TURN; k++) {
        for (l=1; l<=2; l++) {
            int type,ntype=0,otype=0;
            i=k; j = i+TURN+l; if (j>n) continue;
            type = pair[S[i]][S[j]];
            while ((i>=1)&&(j<=n)) {
                if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
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
