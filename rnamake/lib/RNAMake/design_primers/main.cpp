//
//  main.cpp
//  design_primers
//
//  Created by Joseph Yesselman on 3/27/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "types.h"
typedef std::vector<float> Array;
typedef std::vector<Array> Array2d;
typedef std::vector<Array2d> Array3d;


Ints
convert_sequence(
    String const & sequence) {
    
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


bool
matlab_string_compare(
    String const & s1,
    String const & s2) {
    
    int max = (int)(s1.length() > s2.length() ? s2.length() : s1.length());
    Ints seq1 = convert_sequence(s1);
    Ints seq2 = convert_sequence(s2);
    for(int i = 0 ; i < max; i++) {
        if(seq1[i] == seq2[i]) { continue; }
        if(seq1[i] < seq2[i]) { return true; }
        if(seq1[i] > seq2[i]) { return false; }
    }
    
    if(seq1.size() < seq2.size()) { return true; }
    return false;
    //else { return true; }
    
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
    
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return matlab_string_compare(v[i1],v[i2]);});
    
    return idx;
}

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

struct MisprimingResults {
public:
    MisprimingResults() {}
    
public:
    Array num_match_forward, num_match_reverse, best_match_forward, best_match_reverse;
    Array misprime_score_forward, misprime_score_reverse;
    
};


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
            a = 3.92e-5 * ( 0.843 - 0.352*sqrtf(y)*logf(y));
            d = 1.42e-5 * ( 1.279 - 4.03e-3 * logf(y) - 8.03e-3 * (logf( y )*logf(y)));
            g = 8.31e-5 * ( 0.486 - 0.258 * logf(y) + 5.25e-3 * (logf(y)*logf(y)*logf(y)));
        }
    
        x = logf( divalent_concentration );
        Tm_corrected = 1 / ( (1/Tm) + a + b *x  + f_GC*(c + d*x) + (1/(2*(N_BP-1)))*(e + f*x+g*x*x));
        
    }
    
    return Tm_corrected;
    
}


Array2d
precalculate_Tm(
    String const & sequence) {
 
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
    Array2d f_GC (N_BP, Array(N_BP, 0));
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
            if(ns[j] == 1 || ns[j] == 2) {
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
            /*if(i == 1 && j == 18) {
                std::cout <<
            }*/
            Tm[i][j] = ionic_strength_correction(Tm[i][j], monovalent_concentration, divalent_concentration, f_GC[i][j], len_BP[i][j]);
            Tm[i][j] = Tm[i][j] - 273.15;

        }
    }
    return Tm;
    
    

}


String
reverse_complement(
    String const & seq) {
    
    String reverse;
    char c, c2;
    for(int k = 0; k < seq.length(); k++) {
        c = seq [ seq.length() - k - 1 ];
        if     ( c == 'A') { c2 = 'T'; }
        else if( c == 'C') { c2 = 'G'; }
        else if( c == 'G') { c2 = 'C'; }
        else if( c == 'U') { c2 = 'A'; }
        else               { c2 = 'A'; }
        reverse += c2;
    }
    return reverse;
    
}

char
reverse_complement(
    char const & c) {
    
    char c2 = 'T';
    if     ( c == 'A') { c2 = 'T'; }
    else if( c == 'C') { c2 = 'G'; }
    else if( c == 'G') { c2 = 'C'; }
    else if( c == 'U') { c2 = 'A'; }
    else               { c2 = 'A'; }
    return c2;
}


void
output_num_match(
    String const & sequence,
    Array const & num_match_forward,
    Array const & num_match_reverse) {
    
    float COL_SIZE = 150;
    float N_BP = (int)sequence.size();
    int pos;
    int f, r;
    
    std::stringstream ss1, ss2;
    
    for(int i = 0; i < floorf(N_BP/COL_SIZE)+1; i++) {
        String allow_forward_line = "";
        String allow_reverse_line = "";
        String sequence_line = "";
        
        for(int n = 0; n < COL_SIZE; n++) {
            pos = i * COL_SIZE+n;
            if(pos < N_BP) {
                sequence_line += sequence[pos];
                f = num_match_forward[pos] > 9 ? 9 : num_match_forward[pos];
                r = num_match_reverse[pos] > 9 ? 9 : num_match_reverse[pos];
                ss1 << f;
                ss2 << r;
                
                //allowed_forward_line += num
            }
        }
        allow_forward_line = ss1.str();
        allow_reverse_line = ss2.str();
        printf("%s\n%s\n%s\n\n",allow_forward_line.c_str(),sequence_line.c_str(), allow_reverse_line.c_str());
        ss1.str("");
        ss2.str("");
    }
    
}

MisprimingResults
check_for_mispriming_FAST(
    String const & sequence) {
    int N = (int)sequence.size();
    String sub;
    int m = 20;
    Strings subset_strings(2*N);
    int end_pos;
    int start_pos;
    for(int i = 0; i < N; i++) {
        start_pos = 0;
        end_pos = i+1;
        if(i > m) {
            start_pos = i - m;
            end_pos = m+1;
        }
        sub = sequence.substr(start_pos,end_pos);
        std::reverse(sub.begin(), sub.end());
        subset_strings[i] = sub;

    }
    String sequence_rc = reverse_complement(sequence);
    int j = 2*N-1;
    for(int i = 0; i < N; i++) {
        start_pos = 0;
        end_pos = i+1;
        if(i > m) {
            start_pos = i - m;
            end_pos = m+1;
        }
        sub = sequence_rc.substr(start_pos,end_pos);
        std::reverse(sub.begin(), sub.end());
        subset_strings[j] = sub;
        j--;
    }
    
    std::vector<size_t> sortidx = sort_indexes(subset_strings);
    Strings strings_sorted(2*N);
    j = 0;
    for(auto const & i : sortidx) {
        strings_sorted[j] = subset_strings[i];
        j++;
    }

    Ints match_to_next(2*N);
    Array misprime_score_to_next(2*N);
    String string1, string2;
    int count;
    float misprime_score;
    for(int i = 0 ; i < 2*N; i++) {
        count = 0; misprime_score = 0;
        string1 = strings_sorted[i];
        string2 = strings_sorted[i+1];
        while ( count < string1.length() && count < string2.length() ) {
            if(string1[count] != string2[count]) { break; }
            
            if( string1[count] == 'G' || string1[count] == 'C') {
                misprime_score = misprime_score + 1.25;
            }
            else {
                misprime_score = misprime_score + 1;
            }
            count = count + 1;
        }
        match_to_next[i] = count;
        misprime_score_to_next[i] = misprime_score;
    }
    
    Array match_max(2*N);
    Array best_match(2*N);
    Array misprime_score_max(2*N);
    
    match_max[0] = match_to_next[0];
    best_match[0] = 1;
    misprime_score_max[0] = misprime_score_to_next[0];
    
    match_max[2*N-1] = match_to_next[2*N-1];
    best_match[2*N-1] = 2*N-2;
    misprime_score_max[2*N-1] = misprime_score_to_next[2*N-1];
    
    for(int i = 1; i < 2*N; i++) {
        if( match_to_next[i-1] > match_to_next[i]) {
            best_match[i] = i-1;
            match_max[i] = match_to_next[i-1];
            misprime_score_max[i] = misprime_score_to_next[i-1];
        }
        else {
            best_match[i] = i+1;
            match_max[i] = match_to_next[i];
            misprime_score_max[i] = misprime_score_to_next[i];
        }

    }
    
    Array num_match_forward(N);
    Array num_match_reverse(N);
    Array misprime_score_forward(N);
    Array misprime_score_reverse(N);
    Array best_match_forward(N);
    Array best_match_reverse(N);
    float f;
    
    for(int i = 0; i < 2*N; i++) {
        if( sortidx[i] <= N) {
            num_match_forward[sortidx[i]]= match_max[i];
            misprime_score_forward[sortidx[i]] = misprime_score_max[i];
            f = floorf( (float)(sortidx [ best_match[i] ]) / (float)N);
            best_match_forward[sortidx[i]] = (sortidx [ best_match[i]] - 1) - f*N + 1;
        }
        else {
            //std::cout << sortidx[i]-N << " " << match_max[i] << std::endl;
            num_match_reverse[sortidx[i]-N]= match_max[i];
            misprime_score_reverse[sortidx[i]-N] = misprime_score_max[i];
            f = floorf( (float)(sortidx[ best_match[i]-N] - 1) / (float)N);
            best_match_reverse[sortidx[i]-N] = (sortidx [ best_match[i] ] - 1) - f*N + 1;
        }
    }
    //output_num_match(sequence, num_match_forward, num_match_reverse);
    
    MisprimingResults results;
    results.num_match_forward = num_match_forward;
    results.num_match_reverse = num_match_reverse;
    results.misprime_score_forward = misprime_score_forward;
    results.misprime_score_reverse = misprime_score_reverse;
    results.best_match_forward = best_match_forward;
    results.best_match_reverse = best_match_reverse;
    
    return results;
    
    

}


float
get_max(
    Array const & arr1,
    Array const & arr2) {
    
    float max = -1000;
    for (auto const & n : arr1) {
        if(n > max) { max = n; }
    }
    for (auto const & n : arr2) {
        if(n > max) { max = n; }
    }
    
    return max;
}

float
get_min(
    Array3d const & scores,
    int n) {
    
    float min = 1000000;
    for(int i = 0; i < scores.size(); i++) {
        for(int j = 0; j < scores[i].size(); j++) {
            if(scores[i][j][n] < min) {
                min = scores[i][j][n];
            }
        }
    }
    
    return min;
}

inline
float max(
    float a,
    float b) {
    
    return a > b ? a : b;
}

inline
float min(
    float a,
    float b) {
    
    return a > b ? b : a;
}

Strings
output_primers(
    Array2d const & primers_all,
    String const & sequence,
    String const & tag ) {
    
    
    int OUTPUT_STAGGER = 1;
    int num_primers = (int)primers_all[0].size();
    String blank_line;
    int COLWIDTH = 140;
    for(int k = 0; k < max(sequence.length(), COLWIDTH); k++) {
        blank_line += " ";
    }
    String seq_line_prev = blank_line;
    Strings bp_lines(num_primers);
    Strings seq_lines(num_primers);
    Strings primer_sequences(num_primers);
    String seq_line, overlap_seq;
    String bp_line;
    int seq_start, seq_end, direction, last_bp_pos;
    int a;
    
    
    std::stringstream ss;
    String text_out;
    Array primers(primers_all.size());
    for(int j = 0; j < num_primers; j++) {
        for(int i = 0; i < primers_all.size(); i++) {
            primers[i] = primers_all[i][j];
        }
        seq_line = String(blank_line);
        seq_start = primers[0];
        seq_end = primers[1];
        direction = primers[2];
        
        if ( direction == 1) {
            for(int k = seq_start-1; k < seq_end; k++) {
                seq_line[k] = sequence[k];
            }
            primer_sequences.push_back(sequence.substr(seq_start,seq_start+seq_end));
            if(seq_end + 1 <= sequence.length()) {
                seq_line[seq_end+1] = '-';
            }
            if(seq_end + 2 <= sequence.length()) {
                seq_line[seq_end+2] = '>';
            }
            ss << j;
            text_out = ss.str();
            ss.str("");
            a = 0;
            if(seq_end + 2 + text_out.size() <= sequence.length()) {
                for(int k = seq_end + 2; k < seq_end + 2 + text_out.size(); k++) {
                    seq_line[k] = text_out[a];
                    a++;
                }
            }
        }
        
        else {
            for(int k = seq_start-1; k < seq_end; k++) {
                seq_line[k] = reverse_complement(sequence[k]);
            }

            primer_sequences.push_back(reverse_complement(sequence.substr(seq_start,seq_start+seq_end)));
            if(seq_end + 1 <= sequence.length()) {
                seq_line[seq_start-1] = '-';
            }
            if(seq_end + 2 <= sequence.length()) {
                seq_line[seq_start-2] = '<';
            }
            ss << j;
            text_out = ss.str();
            ss.str("");
            a = 0;
            if(seq_end + 2 + text_out.size() <= sequence.length()) {
                for(int k = seq_end - 2; k > seq_end -2 - text_out.size(); k--) {
                    seq_line[k] = text_out[a];
                    a++;
                }
            }
        }
        
        bp_line = String(blank_line);
        overlap_seq = "";
        last_bp_pos = 1;
        
        for(int k = 0; k < sequence.length(); k++) {
            if(seq_line_prev[k] == 'A' || seq_line_prev[k] == 'C' ||
               seq_line_prev[k] == 'G' || seq_line_prev[k] == 'T' ||
               seq_line[k] == 'A' || seq_line[k] == 'C' ||
               seq_line[k] == 'G' || seq_line[k] == 'T') {
                bp_line[k] = ' ';
            }
            else {
                bp_line[k] = '|';
                last_bp_pos = k;
                overlap_seq += sequence[k];
            }
        }
        
        // implement Tm calc
        if(last_bp_pos > 1) {
            
        }
        
        bp_lines.push_back(bp_line);
        seq_lines.push_back(seq_line);
        seq_line_prev = seq_line;
        
    }
    
    int start_pos, end_pos;
    String out_line;
    
    for(int n = 0; n < floorf(((float)(sequence.length()-1))/(float)COLWIDTH)+1; n++) {
        start_pos = COLWIDTH*n + 1;
        end_pos = min( COLWIDTH*n + COLWIDTH, sequence.length());
        out_line = sequence.substr(start_pos,start_pos+end_pos);
        printf("%s\n", out_line.c_str());
        for(int k = 0; k < seq_lines.size(); k++) {
            bp_line = bp_lines[k];
            seq_line = seq_lines[k];
            printf("%s\n%s\n",bp_line.c_str(), seq_line.c_str());
        }
    }
    
    return primer_sequences;
    
}


Strings
design_primers(
    String seq) {
    
    Strings primer_strs;
    float min_Tm = 60;
    int MAX_LENGTH = 60;
    int MIN_LENGTH = 15;
    int NUM_PRIMERS = 0;
    int misprime_mode = 0;
    int is_server = 0;
    
    int N_BP = (int)seq.size();
    
    Array2d Tm_precalculated = precalculate_Tm(seq);
    Array allow_forward_endpoint ( N_BP, 1);
    Array allow_reverse_endpoint ( N_BP, 1);
    MisprimingResults r = check_for_mispriming_FAST(seq);
    
    if(misprime_mode != 0) { std::cout << "misprime_mode not implemented" << std::endl; exit(0); }
    float misprime_score_weight = 10.0;
    int NEW_SCORE = 1;
    //printf("Doing dynamics programming calculation ...");
    int num_primer_sets = (NUM_PRIMERS/2);
    int num_primer_sets_max = 7; // N_BP / MIN_LENGTH
    float MAX_SCORE = N_BP*2 + 1;
    if (NEW_SCORE) {
        MAX_SCORE = MAX_SCORE + misprime_score_weight * get_max(r.misprime_score_forward, r.misprime_score_reverse) * 2 * num_primer_sets_max;
    }

    Array3d scores_start (N_BP, Array2d(N_BP, Array(num_primer_sets_max, MAX_SCORE)));
    Array3d scores_stop (N_BP, Array2d(N_BP, Array(num_primer_sets_max, MAX_SCORE)));
    Array3d scores_final (N_BP, Array2d(N_BP, Array(num_primer_sets_max, MAX_SCORE)));
    
    Array3d choice_start_p (N_BP, Array2d(N_BP, Array(num_primer_sets_max, 0)));
    Array3d choice_start_q (N_BP, Array2d(N_BP, Array(num_primer_sets_max, 0)));
    Array3d choice_start_i (N_BP, Array2d(N_BP, Array(num_primer_sets_max, 0)));
    Array3d choice_start_j (N_BP, Array2d(N_BP, Array(num_primer_sets_max, 0)));

    float q_min, q_max;
    
    for(int p = MIN_LENGTH; p < MAX_LENGTH; p++) {
        if(!allow_forward_endpoint[p]) { continue; }
        
        q_min = max(0, p - MAX_LENGTH);
        q_max = p;
        
        for (int q = q_min; q < q_max; q++) {
            if(!allow_reverse_endpoint[q]) { continue; }
            
            if( Tm_precalculated[q][p] > min_Tm)  {
                scores_stop[p][q][0] = (q-1) + 2 * (p - q + 1);
                if(NEW_SCORE) {
                    scores_stop[p][q][0] = scores_stop[p][q][0] + misprime_score_weight * ( r.misprime_score_forward[p] +  r.misprime_score_reverse[q] );
                }
            }
            
        }
        
    }
    
    float best_min_score = MAX_SCORE;
    int n = 0, best_n = 0;
    int i, j, min_j, max_j, min_i, max_i, min_p, max_p, max_q, min_q;
    int last_primer_length;
    float min_score;
    float potential_score;
    while (n <= num_primer_sets_max) {
        for(int p = 0; p < N_BP; p++) {
            if(!allow_forward_endpoint[p]) { continue; }
        
            q_min = max(0, p - MAX_LENGTH);
            q_max = p;
        
            for (int q = q_min; q < q_max; q++) {
                if(!allow_reverse_endpoint[q]) { continue; }
                if(scores_stop[p][q][n] < MAX_SCORE) {
                    i = N_BP + 1;
                    j = N_BP;
                    last_primer_length = j - q + 1;
                    if (last_primer_length <= MAX_LENGTH && last_primer_length >= MIN_LENGTH) {
                        scores_final[p][q][n] = scores_stop[p][q][n] + (i - p - 1);
                        if(NEW_SCORE) {
                            scores_final[p][q][n] = scores_final[p][q][n] + misprime_score_weight * ( r.misprime_score_forward[p] +  r.misprime_score_reverse[q] );
                        }
                    }
                }
            }
        }
        
        min_score = get_min(scores_final, n);
        
        if((min_score < best_min_score || n == 0) && min_score > 0.1) {
            best_min_score = min_score;
            best_n = n;
        }
        
        if ( n >= num_primer_sets_max) { break; }
        if ( num_primer_sets > 0 && n == num_primer_sets) { break; }
        
        n += 1;
        
        for(int p = 0; p < N_BP; p++) {
            
            q_min = max(0, p - MAX_LENGTH);
            q_max = p;
            
            for(int q = q_min; q < q_max; q++) {
                
                if( scores_stop[p][q][n-1] < MAX_SCORE) {
                    
                    min_j = max( (p+1), q + MIN_LENGTH - 1);
                    max_j = min(  N_BP, q + MAX_LENGTH - 1);
                    for(j = min_j; j < max_j; j++) {
                        
                        min_i = max( (p+1), j - MAX_LENGTH + 1);
                        max_i = j;
                        
                        for(i = min_i; i < max_i; i++) {
                            if(Tm_precalculated[i][j] > min_Tm) {
                                potential_score = scores_stop[p][q][n-1] + (i - p - 1) + 2 * (j - i + 1 );
                                if(potential_score < scores_start[i][j][n-1]) {
                                    scores_start[i][j][n-1] = potential_score;
                                    choice_start_p[i][j][n-1] = p;
                                    choice_start_q[i][j][n-1] = q;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        for(j = 0; j < N_BP; j++) {
            min_i = max(0, j - MAX_LENGTH + 1);
            max_i = j;
            for(i = min_i; i < max_i; i++) {
                if(scores_start[i][j][n-1] < MAX_SCORE) {
                    min_p = max( (j+1), i + MIN_LENGTH - 1);
                    max_p = min( N_BP , i + MAX_LENGTH - 1);
                    
                    for(int p = min_p; p < max_p; p++) {
                        if(!allow_forward_endpoint[p]) { continue; }
                        
                        min_q = max ((j+1), p - MAX_LENGTH + 1);
                        max_q = p;
                        
                        for(int q = min_q; q < max_q; q++) {
                            if(!allow_reverse_endpoint[q]) { continue; }
                            if( Tm_precalculated[q][p] > min_Tm) {
                                potential_score = scores_start[i][j][n-1] + (q - j - 1) + 2 * (p - q + 1 );
                                if (NEW_SCORE) {
                                    potential_score = potential_score + misprime_score_weight * ( r.misprime_score_forward[p] +  r.misprime_score_reverse[q] );
                                }
                                if( potential_score < scores_stop[p][q][n]) {
                                    scores_stop[p][q][n] = potential_score;
                                    choice_start_i[p][q][n] = i;
                                    choice_start_j[p][q][n] = j;
                                }
                            }
                        }
                    }
                }
            }
        }
        
    }
    
    if( num_primer_sets > 0) {
        n = num_primer_sets;
    }
    else {
        n = best_n;
    }
    int p = 0, q = 0;
    min_score = best_min_score;
    
    if( min_score == MAX_SCORE) {
        std::cout << "No solution found" << std::endl;
        return primer_strs;
    }
    
    
    for(i = 0; i < scores_final.size(); i++) {
        for(j = 0; j < scores_final[i].size(); j++) {
            if (fabs(scores_final[i][j][n] - best_min_score) < 0.01) {
                p = i;
                q = j;
            }
        }
    }
    
    Array2d primers (3, Array(2*(n+1)));
    primers[0].back() = q;
    primers[1].back() = N_BP;
    primers[2].back() = -1;
    
    for(int m = n; m > 1; m--) {
        i = choice_start_i[p][q][m];
        j = choice_start_j[p][q][m];
        primers[0][2*m-1] = i;
        primers[1][2*m-1] = p;
        primers[2][2*m-1] = 1;

        p = choice_start_p[i][j][m-1];
        q = choice_start_q[i][j][m-1];
        primers[0][2*m-2] = q;
        primers[1][2*m-2] = j;
        primers[2][2*m-2] = -1;

        
    }
    
    primers[0][0] = 1;
    primers[1][0] = p;
    primers[2][0] = 1;
    String tag = "primer";
    
    primer_strs = output_primers( primers, seq, tag );

    
    return primer_strs;
    
}






int main(int argc, const char * argv[]) {
    // insert code here...
    String seq = "GTAAAACAGGGGTGGCTGCACTGGGCAATCCTCTTGCCGAACAGTGATGCGCCCCACGGAACAAGATAACACGTAGCAGCAGCTTGTATGGTATATAGTC";
    design_primers(seq);
    
    return 0;
}

















