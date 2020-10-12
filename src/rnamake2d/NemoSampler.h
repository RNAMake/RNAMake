#ifndef __RNAMAKE_NEMO_SAMPLER_H__
#define __RNAMAKE_NEMO_SAMPLER_H__

#include <iostream>

#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <climits>
#include <unistd.h>
#include <rnamake2d/Design.h>

extern "C" {
#include <RNAstruct.h>
#include <fold_vars.h>
#include <fold.h>
#include <params.h>
#include <part_func.h>
#include <utils.h>
#include <convert_epars.h>
#include <read_epars.h>
#include <MEA.h>
}

namespace rnamake2d {
    class NemoSampler {
        int curr_design;
        int    verbosity = 0;
        bool   standard_nmcs = false;
        int    iter = 2500;
        char*  target = nullptr;
        char*  start = nullptr;
        char*  seed = nullptr;
        short* pt = nullptr;       // Pair Table
        short* mt = nullptr;       // Mismatch Table
        int*   lt = nullptr;       // Loop Table
        int*   jct = nullptr;      // Junction Table
        int*   smap = nullptr;     // Strength Map
        int    npairs;
        int*   shuffle = nullptr;
        bool   found = false;

        bool   force_boost = false;

    public:
        NemoSampler() : curr_design(-1) {
            init_rand();
            time_t now = time( nullptr );
            now = 1996;
            srand48( now ^ 0x5A5A5A5A );

        }

        void
        initialize_design(Design& design) {

            if(target != nullptr) {
                free(target);
            }

            if(start != nullptr) {
                free(start);
            }

            target = (char*)malloc(sizeof(char)*(1+design.target.size()));
            start = (char*)malloc(sizeof(char)*(1+design.target.size()));

            memcpy(target,design.target.c_str(),design.target.size()+1 );
            memset(start,'N',design.target.size());
            start[design.target.size()] = '\0'; // gotta add that null terminator

            curr_design  = design.id;
            int len = strlen( start );
            pt = make_pair_table( target );
            mt = make_mismatch_table( pt );
            lt = scan_loops( pt );
            jct = scan_junctions( pt, lt );
            smap = make_strength_map( pt, mt );

            char* p;
            for( npairs = 0, p = target; p[npairs]; p[npairs]=='('? (void)npairs++ : (void)p++ );

            shuffle = (int*) calloc( len, sizeof(int) );
            int k;
            for( k = 0; k < len; k++ ) shuffle[k] = k;

            // finally, assign it to the desgin
            //int len = strlen( start );
            char* copy = strdup( start );
            char* position = strdup( start );
            char* secstr = strdup( target );
            secstr[0] = '\0';

            int closest_bpd = 999999;
            char* closest_seq = strdup( start );
            char* closest_struct = strdup( target );
            double closest_fe;

            int stuck = 0;
            char* last_copy = strdup( start );

            double final, e;
            int bpd;
            int n_iter = iter;

            if( seed ) strcpy( copy, seed );

            design.sequence = String{start};
        }

        void
        mutate(Design& design) const  {
            if(design.id != curr_design) {
                std::cerr<<"incorrect design id... exiting";
                exit(1);
            }
        }

    private:
        short* make_mismatch_table( short* pt )
        {
            short* mt = (short*)calloc( 1+pt[0], sizeof(short) );
            int j;

            for( j = 1; j <= pt[0]; j++ ) {
                if( pt[j] ) {
                    if( j < pt[0] && pt[j] > 1 && pt[j+1] && pt[pt[j]-1]
                        && pt[j+1] != pt[j]-1 ) {
                        mt[j+1] = pt[j] - 1;
                        mt[pt[j] - 1] = j+1;
                    }
                    continue;
                }
                if( j < pt[0] && pt[j+1] && pt[j+1]+1 <= pt[0] ) {
                    mt[j] = pt[j+1] + 1;
                    mt[pt[j+1] + 1] = j;
                } else if( j > 1 && pt[j-1] && pt[j-1]-1 >= 1 ) {
                    mt[j] = pt[j-1] - 1;
                    mt[pt[j-1] - 1] = j;
                }
            }

            return mt;
        }

        int* scan_loops( short* pt )
        {
            int len = pt[0];
            int* map = (int*)calloc( 1+len, sizeof(int) );
            int mark = 1;
            int j;

            for( j = 1; j < len; j++ ) {
                if( pt[j] <= 1 ) continue;
                if( map[j] ) continue;
                if( pt[pt[j]-1] != j+1 ) {
                    map[j] = mark;
                    int k = j;
                    do {
                        do {
                            if( ++k > len ) k = 1;
                        } while( pt[k]==0 );
                        k = pt[k];
                        map[k] = mark;
                    } while( k != j );
                    mark++;
                }
            }

            map[0] = mark;
            return map;
        }
        int* scan_junctions( short* pt, int* lt )
        {
            int len = pt[0];
            int* map = (int*)calloc( 1+len, sizeof(int) );
            memmove( map, lt, (1+len)*sizeof(int) );

            int mark = lt[0];

            while( --mark ) {
                int j;
                int count = 0;
                for( j = 1; j <= len; j++ ) if( map[j] == mark ) count++;
                if( count <= 2 ) {
                    for( j = 1; j <= len; j++ ) if( map[j] == mark ) map[j] = 0;
                }
            }

            return map;
        }



// The "strength map" encodes a rough measure of the "fragility" of the
// structure at the considered index. Long helices are "strong", isolated
// base pairs, junctions, etc, are "weak".
//
// The map comes into play when deciding which side in a misfolded pair
// should be mutated. The heuristic is that a "strong" domain should be
// able to take a mutation more easily than a weaker one.
//
        int* make_strength_map( short* pt, short* mt )
        {
            int len = pt[0];
            int* map = (int*)calloc( 1+len, sizeof(int) );
            int j;

            for( j = 1; j <= len; j++ ) {
                int i, k;
                if( pt[j] ) {
                    for( i = 1; i < 3; i++ ) {
                        k = j - i;
                        if( k >= 1 ) map[k] += 3 - i;
                        k = j + i;
                        if( k <= len ) map[k] += 3 - i;
                    }
                } else if( mt[j] ) {
                    int up = j;
                    while( up <= len && pt[up]==0 ) up++;
                    int down = j;
                    while( down >= 1 && pt[down]==0 ) down--;

                    int minus = 2;
                    if( up > len || down < 1 ) { // dangling end
                        minus = 4;
                    } else if ( pt[down] > down && pt[up] < up ) { // hairpin
                        minus = 3;
                    } else if ( pt[down] < down && pt[up] > up ) { // multiloop
                        minus = 4;
                    } else { // internal loop I guess...
                        if( pt[down]-pt[up] == 1 ) { // bulge
                            minus = up-down==2 ? 5 : 4;
                        } else {
                            minus = up-down==2 ? 3 : 2;
                        }
                    }
                    for( i = 1; i < minus; i++ ) {
                        k = j - i;
                        if( k >= 1 ) map[k] -= minus - i;
                        k = j + i;
                        if( k <= len ) map[k] -= minus - i;
                    }
                }
            }

            // normalize the array to make the minimum exactly 1
            int m = 999;
            for( j = 1; j <= len; j++ ) if( map[j] < m ) m = map[j];
            for( j = 1; j <= len; j++ ) map[j] += 1 - m;

            return map;
        }
    };
}


#endif // __RNAMAKE_NEMO_SAMPLER_H__
