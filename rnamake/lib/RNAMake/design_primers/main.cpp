//
//  main.cpp
//  design_primers
//
//  Created by Joseph Yesselman on 3/27/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <math.h>
#include "types.h"
typedef std::vector<float> Array;
typedef std::vector<Array> Array2d;
typedef std::vector<Array2d> Array3d;


struct NN_parameters {
public:
    
    NN_parameters(
        float nT = 310.15) {
        T = nT;
        setup_NN();
        setup_closing_penalty();
        setup_mismatch();
        delH_init = 0.2;
        delS_init = -5.7;
    }
    
public:
    

    void
    setup_NN() {
        float data1[4][4] =
        {
            {-7.6, -8.4,  -7.8, -7.2},
            {-8.5, -8.0,  -9.8, -7.8},
            {-8.2, -9.8,  -8.0, -8.4},
            {-7.2, -8.2,  -8.5, -7.6}
        };
        delH_NN = Array2d(4,std::vector<float>(4));
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                delH_NN[i][j] = data1[i][j];
            }
        }
        
        float data2[4][4] =
        {
            {-21.3, -22.4, -21.0, -20.4},
            {-22.7, -19.9, -24.4, -21.0},
            {-22.2, -24.4, -19.9, -22.4},
            {-21.3, -22.2, -22.7, -21.3}
        };
        delS_NN = Array2d(4,std::vector<float>(4));
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                delS_NN[i][j] = data2[i][j];
            }
        }
        
        delG_NN =  Array2d(4,std::vector<float>(4));
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                delG_NN[i][j] = delH_NN[i][j] - (T * delS_NN[i][j]) / 1000;
            }
        }
        
        
        
        

    }
    
    void
    setup_closing_penalty() {
        delH_AT_closing_penalty = Array(4);
        delS_AT_closing_penalty = Array(4);
        delG_AT_closing_penalty = Array(4);
        
        float data1[4] = { 2.2, 0, 0, 2.2 };
        float data2[4] = { 6.9, 0, 0, 6.9 };

        int i = 0;
        for(i = 0; i < 4; i++) { delH_AT_closing_penalty[i] = data1[i]; }
        for(i = 0; i < 4; i++) { delS_AT_closing_penalty[i] = data2[i]; }

        for(i = 0; i < 4; i++) {
            delG_AT_closing_penalty[i] = delH_AT_closing_penalty[i] - (T * delS_AT_closing_penalty[i]) / 1000;
        }
        
    }
    
    void
    setup_mismatch() {
        delH_mismatch = Array3d(4, Array2d(4, Array(4)));
        delS_mismatch = Array3d(4, Array2d(4, Array(4)));
        delG_mismatch = Array3d(4, Array2d(4, Array(4)));

        //Delta H
        float data1[4][4] =
        {
            {1.2, 2.3, -0.6, -7.6},
            {5.3, 0.0, -8.4,  0.7},
            {-0.7, -7.8, -3.1, 1.0},
            {-7.2, -1.2, -2.5, -2.7}
        };
        float data2[4][4] =
        {
            {-0.9, 1.9, -0.7, -8.5},
            { 0.6, -1.5, -8.0, -0.8},
            {-4.0, -10.6, -4.9, -4.1},
            {-7.8, -1.5, -2.8, -5.0}
        };
        float data3[4][4] =
        {
            {-2.9, 5.2, -0.6, -8.2},
            { -0.7, 3.6, -9.8, 2.3},
            {0.5, -8.0, -6.0, 3.3},
            {-8.4, 5.2, -4.4, -2.2}
        };
        float data4[4][4] =
        {
            {4.7, 3.4, 0.7, -7.2},
            {7.6, 6.1, -8.2, 1.2},
            {3.0, -8.5, 1.6, -0.1},
            { -7.6, 1.0, -1.3, 0.2}
        };
        
        //Delta S
        float data5[4][4] =
        {
            {1.7, 4.6, -2.3, -21.3},
            {14.6, -4.4, -22.4, 0.2},
            { -2.3, -21.0, -9.5, 0.9},
            {-20.4, -6.2, -8.3, -10.8}
        };
        float data6[4][4] =
        {
            { -4.2, 3.7, -2.3, -22.7},
            { -0.6, -7.2, -19.9, -4.5},
            {  -13.2, -27.2, -15.3, -11.7},
            {-21.0, -6.1, -8.0, -15.8}
        };
        float data7[4][4] =
        {
            { -9.8, 14.2, -1.0, -22.2},
            { -3.8, 8.9, -24.4, 5.4},
            {  3.2, -19.9, -15.8, 10.4},
            {-22.4, 13.5, -12.3, -8.4}
        };
        float data8[4][4] =
        {
            {  12.9, 8.0, 0.7, -21.3},
            { 20.2, 16.4, -22.2, 0.7},
            {  7.4, -22.7, 3.6, -1.7},
            {7.4, -22.7, 3.6, -1.7}
        };
        
        //Delta G
        float data9[4][4] =
        {
            { 0.61,   0.88,  0.14, -1.0},
            { 0.77,   1.33,  -1.44,  0.64},
            { 0.02,  -1.28, -0.13,  0.71},
            { -0.88, 0.73,   0.07,  0.69}
        };
        float data10[4][4] =
        {
            { 0.43, 0.75, 0.03, -1.45},
            {  0.79, 0.70, -1.84, 0.62},
            { 0.11, -2.17, -0.11, -0.47},
            { -1.28, 0.40, -0.32, -0.13}
        };
        float data11[4][4] =
        {
            { 0.17, 0.81, -0.25, -1.30},
            { 0.47, 0.79, -2.24, 0.62},
            { -0.52, -1.84, -1.11, 0.08},
            { -1.44, 0.98, -0.59, 0.45}
        };
        float data12[4][4] =
        {
            { 0.69, 0.92, 0.42, -0.58},
            {  1.33, 1.05, -1.30, 0.97},
            { 0.74, -1.45, 0.44, 0.43},
            { -1.00, 0.75, 0.34, 0.68}
        };
        
        
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                delH_mismatch[i][j][0] = data1[i][j];
                delH_mismatch[i][j][1] = data2[i][j];
                delH_mismatch[i][j][2] = data3[i][j];
                delH_mismatch[i][j][3] = data4[i][j];
                delS_mismatch[i][j][0] = data5[i][j];
                delS_mismatch[i][j][1] = data6[i][j];
                delS_mismatch[i][j][2] = data7[i][j];
                delS_mismatch[i][j][3] = data8[i][j];
                delG_mismatch[i][j][0] = data9[i][j];
                delG_mismatch[i][j][1] = data10[i][j];
                delG_mismatch[i][j][2] = data11[i][j];
                delG_mismatch[i][j][3] = data12[i][j];
            }
        }

        
    }
    
    
public:
    Array delH_AT_closing_penalty, delS_AT_closing_penalty, delG_AT_closing_penalty;
    Array2d delH_NN, delS_NN, delG_NN;
    Array3d delH_mismatch, delS_mismatch, delG_mismatch;
    float T, delH_init, delS_init;
    
};


Ints
convert_sequence(
    String sequence) {
    
    Ints int_sequence(sequence.size());
    int i =0;
    for(auto const & s : sequence) {
        if(s == 'A') { int_sequence[i] = 0; } //1
        if(s == 'C') { int_sequence[i] = 1; } //2
        if(s == 'G') { int_sequence[i] = 2; } //3
        if(s == 'T') { int_sequence[i] = 3; } //4
        i++;
        
    }
    
    return int_sequence;
    
}

float
ionic_strength_correction(
    float Tm,
    float monovalent_concentration,
    float divalent_concentration,
    float f_GC,
    float N_BP) {
    
    float R = sqrtf(divalent_concentration) / monovalent_concentration;
    float x;
    float Tm_corrected = Tm;
    if( R < 0.22) {
        x = logf( monovalent_concentration );
        Tm_corrected = 1 / ( (1/Tm) + (4.29*f_GC - 3.95)*1e-5 * x + 9.40e-6*x*x);
    }
    else {
        float a = 3.92e-5;
        float b = -9.11e-6;
        float c = 6.26e-5;
        float d = 1.42e-5;
        float e = -4.82e-4;
        float f = 5.25e-4;
        float g = 8.31e-5;
        
        if( R < 6.0 ) {
            float y = monovalent_concentration;
            a = 3.92e-5 * ( 0.843 - 0.352*sqrt(y)*log(y));
            d = 1.42e-5 * ( 1.279 - 4.03e-3 * logf(y) - 8.03e-3 * (logf( y )*logf(y)));
            g = 8.31e-5 * ( 0.486 - 0.258 * logf(y) + 5.25e-3 * (logf(y)*logf(y)*logf(y)));
        }
    
        Tm_corrected = 1 / ( (1/Tm) + a + b *x  + f_GC*(c + d*x) + (1/(2*(N_BP-1)))*(e + f*x+g*x*x));
        
    }
    
    return Tm_corrected;
    
}


Array2d
precalculate_Tm(
    String sequence) {
 
    float DNA_concentration = 0.2e-6;
    float monovalent_concentration = 0.1;
    float divalent_concentration = 0.0015;
    int N_BP = (int)sequence.size();
    Ints ns = convert_sequence(sequence);
    NN_parameters NN;
    float delS_DNAconc = 1.987 * log( DNA_concentration/2 );
    float delS_init =  NN.delS_init + delS_DNAconc;

    Array2d delH_matrix (N_BP, Array(N_BP));
    Array2d delS_matrix (N_BP, Array(N_BP));
    Array2d f_GC (N_BP, Array(N_BP));
    Array2d len_BP (N_BP, Array(N_BP));
    Array2d Tm (N_BP, Array(N_BP));


    for(int i = 0; i < N_BP; i++) {
        if(ns[i] == 1 || ns[i] == 2) {
            f_GC[i][i] = 1;
        }
    }
    for(int i = 0; i < N_BP; i++) {
        for(int j = 0; j < N_BP; j++) {
            delH_matrix[i][j] = NN.delH_init;
            delS_matrix[i][j] = delS_init;
            len_BP[i][j] = 1;
        }
    }
    
    for(int i = 0; i < N_BP; i++) {
        for(int j = i+1; j < N_BP; j++) {
            delH_matrix[i][j] = delH_matrix[i][j-1] + NN.delH_NN [ ns[j-1] ][ ns[j] ];
            delS_matrix[i][j] = delS_matrix[i][j-1] + NN.delS_NN [ ns[j-1] ][ ns[j] ];
            len_BP[i][j] = len_BP[i][j-1] + 1;
            f_GC[i][j] = f_GC[i][j-1];
            if(ns[i] == 1 || ns[i] == 2) {
                f_GC[i][j] = f_GC[i][j] + 1;
            }
        }
    }
    
    
    //ADD_TERMINAL_PENALTY
    for(int i = 0; i < N_BP; i++) {
        for(int j = i+1; j < N_BP; j++) {
            delH_matrix[i][j] = delH_matrix[i][j] + NN.delH_AT_closing_penalty [ ns[i] ];
            delH_matrix[i][j] = delH_matrix[i][j] + NN.delH_AT_closing_penalty [ ns[j] ];
            
            delS_matrix[i][j] = delS_matrix[i][j] + NN.delS_AT_closing_penalty [ ns[i] ];
            delS_matrix[i][j] = delS_matrix[i][j] + NN.delS_AT_closing_penalty [ ns[j] ];

        }
    }
    
    for(int i = 0; i < N_BP; i++) {
        for(int j = 0; j < N_BP; j++) {
            Tm[i][j] = 1000 * (delH_matrix[i][j] / delS_matrix[i][j]);
            f_GC[i][j] = f_GC[i][j] / len_BP[i][j];
        }
    }
    
    //IONIC_STRENGTH_CORRECTION
    for(int i = 0; i < N_BP; i++) {
        for(int j = i; j < N_BP; j++) {
            Tm[i][j] = ionic_strength_correction(Tm[i][j], monovalent_concentration, divalent_concentration, f_GC[i][j], len_BP[i][j]);
            Tm[i][j] = Tm[i][j] - 273.15;
        }
    }
    
    for(int i = 0; i < N_BP; i++) {
        std::cout << Tm[0][i] << " ";
    }
    std::cout << std::endl;
    
    
    return Tm;
    
    

}


Strings
design_primers(
    String seq) {
    
    Strings primers;
    float min_Tm = 60;
    int MAX_LENGTH = 60;
    int MIN_LENGTH = 15;
    int NUM_PRIMERS = 0;
    int misprime_mode = 0;
    int is_server = 0;
    
    int N_BP = (int)seq.size();
    
    Array2d Tm_precalculated = precalculate_Tm(seq);
    
    return primers;
    
}






int main(int argc, const char * argv[]) {
    // insert code here...
    String seq = "TTTTTAAAAAAAATTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    design_primers(seq);
    
    return 0;
}
