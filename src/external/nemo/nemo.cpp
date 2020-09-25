// Copyright (C) 2018 Fernando Portela <nando8888@gmail.com>


#include <cstdlib>
#include <ctime>
#include <cstring>
#include <ctype.h>
#include <cmath>
#include <cstdio>
#include <limits.h>
#include <unistd.h>

// ViennaRNA stuff
extern "C" {
#include "RNAstruct.h"
#include "fold_vars.h"
#include "fold.h"
#include "params.h"
#include "part_func.h"
#include "utils.h"
#include "convert_epars.h"
#include "read_epars.h"
#include "MEA.h"
}

// set to 0 to prevent use of base pair distances in the scoring function
#define USE_BPD 1

// set to 0 to prevent use of free energy differences in the scoring function
#define USE_DDG 1

// set to 0 to disable heuristics in the sampling phase
#define USE_DOMAIN_KNOWLEDGE 1

// A few globals
// FIXME: globals are ugly, refactor whenever possible
//
int    verbosity = 0;
bool   standard_nmcs = false;
int    iter = 2500;
char*  target = NULL;
char*  start = NULL;
char*  seed = NULL;
short* pt = NULL;       // Pair Table
short* mt = NULL;       // Mismatch Table
int*   lt = NULL;       // Loop Table
int*   jct = NULL;      // Junction Table
int*   smap = NULL;     // Strength Map
int    npairs;
int*   shuffle = NULL;
bool   found = false;

bool   force_boost = false;

// --------------------------------------------------------------------------
// Utilities
//
int my_urn( int from, int to )
{
    return ( ( (int) (drand48()*(to-from+1)) ) + from );
}

int pair_map( char a, char b )
{
    if( (a ^ b ^ 'G' ^ 'C') == 0 ) return 3;
    if( (a ^ b ^ 'A' ^ 'U') == 0 ) return 2;
    if( (a ^ b ^ 'G' ^ 'U') == 0 ) return 1;
    return 0;
}

// FIXME: A mismatch table should be 2D, not just 1D
//        In the grand scheme of things, it probably doesn't affect results
//        that much, but there's a clear flaw about using a simple 1D array:
//        all 1-N internal loops (except 1-1 ones) are missing the fact that
//        the lone unpaired base has 2 mismatch partners, not just one.
//
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

// Initialize helper arrays marking the closing pairs in various loops
//
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

// When resetting a position to N, take the pair partner as well, if any
//
void set_N( char* seq, short* pt, int pos )
{
    seq[pos] = 'N';
    if( pt[1+pos] ) seq[pt[1+pos]-1] = 'N';
}


// --------------------------------------------------------------------------
// Playout scoring
//
// The first experimentations with Nested Monte-Carlo searches were done
// with normal_score() and were painfully slow.
// quick_score() was an attempt to use a (much) faster approximation.
// But solving rate dropped by about a 1/3
//
// Dead code, currently
//
double quick_score( char* position )
{
    int length = strlen( target );
    int* map = (int*) calloc( length*length, sizeof(int) );
    
    int i,j,k;
    
    // fill
    for( i = 0; i < length-1; i++ ) {
        for( j = i+1; j < length; j++ ) {
            map[i*length+j] = pair_map( position[i], position[j] );
        }
    }

    double sum = 1;
    int o;
    for( o = 0; o < 2; o++ ) {
        for( i = 2; i < length - 5; i++ ) {
            int prev_len, len = 0;
            int prev_weight, weight = 0;
            int expected = 0;
            int prev_j = -length;
            int last = 0;
            for( j=i, k=i+o+4; j >= 0 && k < length; j--, k++ ) {
                int v = map[j*length+k];
                if( v==0 && last==0 ) continue;
                if( v==0 ) {
                    if( len > 2 && len != expected ) {
                        sum += log(len) * weight;
                    }
                    prev_len = len;
                    len = 0;
                    prev_weight = weight;
                    weight = 0;
                    prev_j = j;
                    expected = 0;
                } else {
                    if( last == 0 ) {
                        if( j - prev_j <= 6 ) {
                            sum -= log(prev_len) * prev_weight;
                            len = prev_len;
                            weight = prev_weight - (j - prev_j - 2);
                        }
                    }
                    if( pt[j+1] == k+1 ) expected++;
                    weight += v;
                    len++;
                    
                    if( j >= 3 && k < length - 3 ) {
                        int x = j*length+k;
                        if( map[x-(length-1)]==1 && map[x-2*(length-1)]==1 && map[x-3*(length-1)]>0
                            && position[j-1] != position[j-2] ) {
                            // crossed UG
                            map[x-(length-1)] = 0;
                            map[x-2*(length-1)] = 0;
                        }
                    }
                }
                last = v;
            }
            if( len > 2 && len != expected ) {
                sum += log(len) * weight;
            }
        }
    }
    
    free( map );
    double score = 1. / (1. + log( sum ) );
    if( 0 ) printf("QS: %s %f\n", position, score );
    return score;
}

    
double normal_score( char* position )
{
    char* secstr = strdup( position );
    secstr[0] = '\0';
    double e = fold( position, secstr );
    int bpd = bp_distance( target, secstr );
#if USE_DDG
    double es = energy_of_structure( position, target, 0 );
#endif
    free( secstr );
    if( bpd == 0 ) {
        found = true;
        if( verbosity > 1 ) printf( "Found match:\n%s\n", position );
        return 1.0;
    }
    double score;
#if USE_BPD
    score = (npairs==0 ? 1.0 / (1.0 + bpd) : 1.0 - (0.5 * bpd) / npairs);
#endif
#if USE_DDG
    double e_factor = 1.01 + es - e;
#if USE_BPD
    score *= (score < 0 ? e_factor : 1.0/e_factor);
#else
    score = 1.0/e_factor;
#endif
#endif
    return score;
}


// --------------------------------------------------------------------------
// Playouts related functions
//

// we may need these some day, or something alike
//
const char* maxCC = "CCC";
const char* maxGG = "GGG";

bool is_legal( char* position )
{
#if 0
    if( maxCC && strstr( position, maxCC ) != NULL ) return false;
    if( maxGG && strstr( position, maxGG ) != NULL ) return false;
#endif
    return true;
}


bool undecided( char x )
{
    return (strchr( "AUGC", x ) == NULL);
}


bool locked( int pos )
{
    return !undecided( start[pos-1] );
}


bool test_move( char* position, int o, char base )
{
    bool ok = true;
    char save = position[o];
    position[o] = base;
    ok = is_legal( position );
    position[o] = save;
    return ok;
}


bool play_move( char* position, const char* bases, int z = -1 )
{
    int l = strspn( position, bases );
    int len = pt[0];
    
    if( pt[l+1] ) { // it's a pair
        if( z < 0 ) { // and we're instructed to choose a random one
            int w[] = { 16, 16, 9, 9, 2, 2 };

#if USE_DOMAIN_KNOWLEDGE
            // Tuning of weights for adjacent stacks in junctions.
            // When looking at adjcent stacks in junctions from the point of view of an
            // observer at the center of the loop, it is usually better to make sure that
            // the left one has a high probability of being GC/CG, while the rightmost one
            // can usually afford to be demoted to AU/UA.
            //
            if( l < len-1 && mt[1+l+1] && pt[1+l+1] && pt[1+l+1] < len && mt[pt[1+l+1]+1] == 1+l) {
                // ok, looks like adjacent stacks, but we need to eliminate the possibility
                // that it could be a 0-N bulge
                if( jct[1+l] || jct[pt[1+l]] ) {
                    w[0] -= 6; w[1] -= 6; w[2] += 6; w[3] += 6;
                }
            }
            if( l > 0 && mt[1+l-1] && pt[1+l-1] && pt[1+l-1] > 1 && mt[pt[1+l-1]-1] == 1+l) {
                if( jct[1+l] || jct[pt[1+l]] ) {
                    w[0] += 6; w[1] += 6; w[2] -= 6; w[3] -= 6;
                }
            }

            // Tuning for triloops
            // Only GC/CG closing pairs work if the supporting pair also is GC/CG
            int d, j, k;
            if( pt[1+l] > (1+l) ) {
                d = pt[1+l] - (1+l);
                j = 1+l + 1;
                k = 1+l - 1;
            } else {
                d = (1+l) - pt[1+l];
                j = pt[1+l] + 1;
                k = pt[1+l] - 1;
            }
            if( d == 4 && pair_map( position[k-1], position[pt[1+k-1]-1] ) == 3 ) {
                w[2] = 0; w[3] = 0; w[4] = 0; w[5] = 0;
            }
            if( d == 6 && pt[j] && pt[j]-j == 4 && !undecided( position[j-1] )
                && pair_map( position[j-1], position[pt[j]-1] ) != 3 ) {
                w[0] = 0; w[1] = 0;
            }
#endif

            if( !test_move( position, l, 'C' ) ) w[0] = 0;
            if( !test_move( position, l, 'G' ) ) { w[1] = 0; w[5] = 0; }
            if( !test_move( position, l, 'U' ) ) { w[2] = 0; w[4] = 0; }
            if( !test_move( position, l, 'A' ) ) w[3] = 0;
            if( !test_move( position, pt[1+l]-1, 'C' ) ) w[1] = 0;
            if( !test_move( position, pt[1+l]-1, 'G' ) ) { w[0] = 0; w[4] = 0; }
            if( !test_move( position, pt[1+l]-1, 'U' ) ) { w[3] = 0; w[5] = 0; }
            if( !test_move( position, pt[1+l]-1, 'A' ) ) w[2] = 0;

            // FIXME: what if all w[]==0 ?            
            int dice = int_urn( 0, w[0]+w[1]+w[2]+w[3]+w[4]+w[5]-1 );
            for( z = 0; dice >= w[z]; z++ ) dice -= w[z];
            // should leave us with 0 <= z <= 5 randomly selected according to weights
            // FIXME: isn't this a perfect spot for asserting?
        }
        switch( z ) {
        case 0:
            position[l] = 'C'; position[pt[l+1]-1] = 'G';
            break;
        case 1:
            position[l] = 'G'; position[pt[l+1]-1] = 'C';
            break;
        case 2:
            position[l] = 'U'; position[pt[l+1]-1] = 'A';
            break;
        case 3:
            position[l] = 'A'; position[pt[l+1]-1] = 'U';
            break;
        case 4:
            position[l] = 'U'; position[pt[l+1]-1] = 'G';
            break;
        case 5:
            position[l] = 'G'; position[pt[l+1]-1] = 'U';
            break;
        default:
            return false;
        }
    } else { // unpaired, so choose a random base
        int w[] = { 93, 1, 5, 1 };
        if( z < 0 && mt[l+1] ) { // is mismatched
            if( pt[mt[l+1]] ) { // mismatch is paired
                // FIXME: this would probably be a good place to use some (ML-based?) optimization
                //        for weights fine-tuning
                switch( position[mt[l+1]-1] ) {
                case 'A': w[0]=5; w[1]=0; w[2]=2; w[3]=1; break;
                case 'U': w[0]=0; w[1]=6; w[2]=1; w[3]=4; break;
                case 'G': w[0]=2; w[1]=1; w[2]=5; w[3]=0; break;
                case 'C': w[0]=6; w[1]=4; w[2]=0; w[3]=1; break;
                }
            } else { // mismatch is unpaired
                int dice;
                int up, down, close;
                if( l > 0 && pt[1+l-1] == mt[1+l]+1 ) {
                    up = l+1;
                    down = mt[1+l];
                } else {
                    down = l+1;
                    up = mt[1+l];
                }
                close = up - 1;
                int u, d;
                for( u = 0; up+u <= len && pt[up+u]==0; u++ );
                for( d = 0; down > d && pt[down-d]==0; d++ );
                bool internal_loop = ( up+u <= len && down > d && pt[up+u] == down-d );
                //
                // once computations to identify the motif are done,
                // we try to apply some domain-related knowledge
                //
                // FIXME: the above structural analysis should probably be done once
                //        at the start since it only depends on the pt and mt arrays
                //        which stay static throughout the execution
                //
                // FIXME: these "tricks" below should be learned by some ML algorithm
                //
                if( internal_loop && down-d == close ) { // unbranched & single segment = hairpin
#if USE_DOMAIN_KNOWLEDGE
                    // A test for a potential 'slide' from a triloop to a GAAA tetraloop
                    // If matched, have a good chance to apply the "anti-boost" (U/C in the middle
                    // of the triloop)
                    if( l > 1 && d == 3 && undecided(position[l+1]) && position[l-1] == 'G' 
                        && pair_map(position[l-2], position[l+3]) > 0
                        && (force_boost || int_urn(0,1)==0) ) {
                        dice = int_urn(0, 9);
                        if( dice < 5 ) {
                            position[l+1] = 'U';
                        } else if( dice < 9 ) {
                            position[l+1] = 'C';
                        }
                    }
                    // Simple apical loop G/A boosting
                    if( undecided(position[mt[1+l]-1]) && d > 3
                        && (force_boost || int_urn(0,5)!=0) ) {
                        if( test_move( position, l, 'G' )
                            && test_move( position, mt[1+l]-1, 'A' ) ) {
                            z = 2;
                            position[mt[1+l]-1] = 'A';
                        }
                    }
#endif
                } else if( internal_loop ) {
                    if( mt[l+1] > l+1 ) {
                        w[0]=5; w[1]=1; w[2]=20; w[3]=1;
                    } else {
                        w[0]=21; w[1]=1; w[2]=4; w[3]=1;
                    }
                    if( undecided(position[mt[1+l]-1]) ) {
#if USE_DOMAIN_KNOWLEDGE
                        // 1-1 internal loop
                        if( u == 1 && d == 1 ) {
                            dice = int_urn(0, 9);
                            if( force_boost || dice < 8 ) {
                                if( test_move( position, l, 'G' )
                                    && test_move( position, mt[1+l]-1, 'G' ) ) {
                                    z = 2;
                                    position[mt[1+l]-1] = 'G';
                                }
                            }
                        }
                        // 2-2 internal loop
                        if( u == 2 && d == 2 ) {
                            dice = int_urn(0, 1);
                            if( (force_boost || dice == 0) && undecided(position[l+1])
                                && undecided(position[mt[1+l]-1-1]) ) {
                                if( test_move( position, l, 'U' )
                                    && test_move( position, mt[1+l]-1, 'G' ) 
                                    && test_move( position, l+1, 'G' ) 
                                    && test_move( position, mt[1+l]-1-1, 'U' ) ) {
                                    z = 1;
                                    position[mt[1+l]-1] = 'G';
                                    position[l+1] = 'G';
                                    position[mt[1+l]-1-1] = 'U';
                                }
                            }
                        }
                        // if a selection hasn't been made yet, try typical boosts for
                        // internal loops (G/A, A/G, U/U)
                        if( z < 0 ) {
                            dice = int_urn(0, 9);
                            if( dice < 3 ) {
                                z = 2;
                                position[mt[1+l]-1] = 'A';
                            } else if( dice < 6 ) {
                                z = 0;
                                position[mt[1+l]-1] = 'G';
                            } else if( dice < 7 ) {
                                z = 1;
                                position[mt[1+l]-1] = 'U';
                            }
                        }
#endif
                    } else {
                        switch( position[mt[1+l]-1] ) {
                        case 'A': w[0]=4; w[1]=0; w[2]=4; w[3]=1; break;
                        case 'U': w[0]=0; w[1]=6; w[2]=1; w[3]=2; break;
                        case 'G': w[0]=6; w[1]=1; w[2]=2; w[3]=0; break;
                        case 'C': w[0]=4; w[1]=1; w[2]=0; w[3]=1; break;
                        }
                    }
                } else { // a mismatch in a junction or external loop
                    w[0]=97; w[1]=1; w[2]=1; w[3]=1;
#if USE_DOMAIN_KNOWLEDGE
                    if( down == 1+l ) {
                        if( (l < 1 || pt[1+l-1]==0) && position[l+1] == 'G' && position[pt[1+l+1]-1] == 'C' ) {
                            w[3]=48; // increase chance of C boost
                        }
                    } else {
                        if( (l >= len-1 || pt[1+l+1]==0) && position[l-1] == 'G' && position[pt[1+l-1]-1] == 'C' ) {
                            w[2]=48; // increase chance of G boost
                        }
                    }
#endif
                }
            }
        }

        if( z < 0 ) {
            if( !test_move( position, l, 'A' ) ) w[0] = 0;
            if( !test_move( position, l, 'U' ) ) w[1] = 0;
            if( !test_move( position, l, 'G' ) ) w[2] = 0;
            if( !test_move( position, l, 'C' ) ) w[3] = 0;
            int dice = int_urn( 0, w[0]+w[1]+w[2]+w[3]-1 );
            for( z = 0; dice >= w[z]; z++ ) dice -= w[z];
        }
        switch( z ) {
        case 0: position[l] = 'A'; break;
        case 1: position[l] = 'U'; break;
        case 2: position[l] = 'G'; break;
        case 3: position[l] = 'C'; break;
        }
        
    }
    return is_legal( position );
}


// --------------------------------------------------------------------------
// The "core" of the Nested Monte-Carlo algorithm
//

double sample( char* position )
{
    int k;
    for( k = 0; k < pt[0]; k++ ) {
        if( position[k]=='N' && pt[1+k] ) position[k] = 'P';
    }
    // first fill the target pairs
    char bases[] = "AUGCN";
    while( strspn( position, bases ) != strlen( position ) ) {
        play_move( position, bases, -1 );
    }
    // now fill the rest, i.e. the unpaired bases
    bases[4] = 0;
    while( strspn( position, bases ) != strlen( position ) ) {
        play_move( position, bases, -1 );
    }
    return normal_score( position );
}


#define WORST_SCORE -100000.0

double nested( char* position, int level )
{
    if( strspn( position, "AUGC" ) == strlen( position ) ) {
        return normal_score( position );
    }

    double best = WORST_SCORE;
    char* best_playout = strdup( position );
    char* best_local = strdup( position );

    while( strspn( position, "AUGC" ) != strlen( position ) ) {
        int l = strspn( position, "AUGC" );
        double max = WORST_SCORE;
        int z;
        int zm = pt[l+1]==0? 4 : 6;
        if( level == 1 ) {
            if (verbosity > 3) printf("==== %s\n", position);
            #pragma omp parallel for schedule(dynamic)
            for( z = 0; z < zm; z++ ) {
                char* playout = strdup( position );
                if( !play_move( playout, "AUGC", z ) ) continue;
                double v = sample( playout );
                if (verbosity > 3) printf("---- %s %f\n", playout, v);
                #pragma omp critical(max_update)
                {
                    if( v > max ) {
                        max = v;
                        strcpy( best_local, playout );
                    }
                }
                free( playout );
            }
        } else {
            for( z = 0; z < zm; z++ ) {
                char* playout = strdup( position );
                if( !play_move( playout, "AUGC", z ) ) continue;
                double v = nested( playout, level - 1 );
                if( v > max ) {
                    max = v;
                    strcpy( best_local, playout );
                }
                free( playout );
            }
        }

        if( standard_nmcs || (max > best) ) {
            best = max;
            strcpy( best_playout, best_local );
        }
        
        if( found ) {
            strcpy( position, best_playout );
            break;
        }
        
        position[l] = best_playout[l];
        if( pt[l+1] ) {
            position[pt[l+1]-1] = best_playout[pt[l+1]-1];
        } else if( mt[l+1] && position[mt[l+1]-1]=='N' && pt[mt[l+1]-1]==0 ) {
            position[mt[l+1]-1] = best_playout[mt[l+1]-1];
        }
    }

    free( best_playout );
    free( best_local );
    return best;
}

// --------------------------------------------------------------------------
// Main body of the program
//

void config_eterna( void )
{
    // convert_parameter_file("vrna185.par", "vrna185x2.par", VRNA_CONVERT_OUTPUT_ALL);
    
    read_parameter_file("vrna185x.par");
    dangles = 1;
    update_fold_params(); // useful ?
}


void usage( const char* cmd )
{
    printf( "Usage: %s [-E] [-v[v[v]]] [-i <n_iter>] <target_struct> [<start_sequence>]\n", cmd );
}


// FIXME: use something proper like getopt (?)
//
bool parse_arguments( int argc, char** argv )
{
    char* command = argv[0];

    while( argc > 1 && argv[1][0] == '-' ) {
        if( strcmp( argv[1], "-E" )==0 ) {
            config_eterna();
            argc--;
            argv++;
        } else if( strcmp( argv[1], "-S" )==0 ) {
            standard_nmcs = true;
            argc--;
            argv++;
        } else if( strcmp( argv[1], "-i" )==0 ) {
            iter = atoi( argv[2] );
            argc -= 2;
            argv += 2;
        } else if( strncmp( argv[1], "-v", 2 )==0 ) {
            char* p = &(argv[1][1]);
            for( /* */; (*p)=='v'; p++ ) verbosity++;
            argc--;
            argv++;
        } else {
            printf( "Unknow argument '%s'.\n", argv[1] );
            usage( command );
            return false;
        }
    }

    if( argc > 1 ) {
        target = strdup( argv[1] );
        if( strspn( target, "(.)" ) != strlen(target) ) {
            printf( "Invalid character '%c' in structure.\n",  target[strspn( target, "(.)" )] );
            return false;
        }
        if( argc > 2 ) {
            start = strdup( argv[2] );
            if( strlen(start) != strlen(target) ) {
                printf( "Sequence length doesn't match target structure.\n" );
                return false;
            }
            if( strspn( start, "AUGCN" ) != strlen(start) ) {
                printf( "Invalid character '%c' in sequence.\n",  start[strspn( start, "AUGCN" )] );
                return false;
            }
            if( argc > 3 ) {
                seed = strdup( argv[3] );
                if( strlen(seed) != strlen(target) ) {
                    printf( "Seed sequence length doesn't match target structure.\n" );
                    return false;
                }
                if( strspn( seed, "AUGC" ) != strlen(start) ) {
                    printf( "Invalid character '%c' in seed sequence.\n",  seed[strspn( seed, "AUGC" )] );
                    return false;
                }
            }
        } else {
            start = strdup( argv[1] );
            memset( start, 'N', strlen( target ) );
        }
    } else {
        usage( command );
        return false;
    }

    return true;
}


void init_globals( void )
{
    init_rand();
    time_t now = time( NULL );
    srand48( now ^ 0x5A5A5A5A );

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

    if( verbosity > 2 ) {
        printf( "   " );
        for( k = 1; k <= len; k++ ) {
            printf( "%c", pt[k]? ( pt[k] > k? (mt[k]? '[' : '(') : (mt[k]? ']' : ')' ) ) : mt[k]? '*' : '.' );
        }
        printf("\n   ");
        for( k = 1; k <= len; k++ ) {
            printf( "%c", lt[k]? (jct[k]? 'Y' : '|') : '.' );
        }
        printf("\n");
        fflush( stdout );
    }
}

/*
int main( int argc, char** argv )
{
    if( !parse_arguments( argc, argv ) ) return 1;

    init_globals();

    int len = strlen( start );
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

    for( ; iter; iter-- ) {
        strcpy( position, copy );
        // half the time, use the best boosts we know of during playouts
        force_boost = int_urn(0,1)==0;

        // (try to) solve
        final = nested( position, 1 );

        // score the search result
        secstr[0] = 0;
        e = fold( position, secstr );
        bpd = bp_distance( target, secstr );
        if( verbosity > 0 ) {
            printf( " %f %.2f %d (%d)\n", final, e, bpd, iter ); fflush( stdout );
        }
        if( bpd == 0 ) break;
        
        if( bpd < closest_bpd ) {
            closest_bpd = bpd;
            strcpy( closest_seq, position );
            strcpy( closest_struct, secstr );
            closest_fe = e;
        }

        // Our previous attempt failed. Rather than simply increase depth, or just
        // try again with the original starting point, we build a new starting point based
        // on what seems to have worked and what didn't (misfolds)

        int i, j, k, c;

        // let's get a pair map of the misfolded structure
        short* pa = make_pair_table( secstr );

        bool* retry = (bool*)calloc( 1+len, sizeof(bool) );
        for( j = 1; j <= len; j++ ) retry[j] = (pt[j] != pa[j]);

        for( j = 1; j <= len; j++ ) {
            // add mismatches of the misfolded bases
            if( !retry[j] && mt[j] && retry[mt[j]] ) retry[j] = true;
        }

        short* ma = make_mismatch_table( pa );
        int*   la = scan_loops( pa );

        for( j = 1; j < la[0]; j++ ) {
            for( i = 1; i <= len && la[i] != j; i++ ) ;
            bool has_opened_pair = false;
            bool closing_good = true;
            k = i;
            do {
                do {
                    if( ++k > len ) k = 1;
                    if( pa[k]==0 && pt[k] ) has_opened_pair = true;
                } while( pa[k]==0 );
                if( pa[k] != pt[k] ) closing_good = false;
                k = pa[k];
            } while( k != i && closing_good );

            if( closing_good && has_opened_pair ) {
                do {
                    do {
                        if( ++k > len ) k = 1;
                        if( pa[k] || ma[k] ) retry[k] = true;
                    } while( pa[k]==0 );
                    k = pa[k];
                    retry[k] = true;
                } while( k != i );
            }
        }

        // we're ready

        do {

            strcpy( copy, position );

            for( k = 0; k < len; k++ ) {
                // choose random r, with r != k
                int r = int_urn(0, len-2);
                if( r >= k ) r++;
                // cool exercise for job interviews: ask the candidate what these 3 lines do.
                shuffle[k] ^= shuffle[r];
                shuffle[r] ^= shuffle[k];
                shuffle[k] ^= shuffle[r];
                // full points for the correct answer
                // bonus point for finding the answer without writing anything down
                // bonus point if the candidate declares "this is evil / bad style"
                // bonus point if the candidate can cite a context where this trick
                //             might be justified or useful
            }
                    
            for( k = 0, c = 0; k < len; k++ ) {
                i = shuffle[k];

                // not a misfolded (or otherwise interesting) spot? let's keep it
                if( !retry[1+i] ) continue;

                // we don't want to reset and retry all misfolded bases and pairs
                // so we use probabilities, first 1/1, then 1/2, 1/3 and so on
                // this guarantess that we will reset at least one base or pair,
                // and maybe a few more.
                if( int_urn(0,c++)!=0 ) continue;

                // FIXME: Ad hoc rules, based on personal experience.
                // Could probably use some ML improvements as well...

                if( copy[i] != 'N' && pa[1+i] && copy[pa[1+i]-1] != 'N' ) {
                    if( pt[1+i]==0 && pt[pa[1+i]]==0 ) {
                        set_N( copy, pt, i );
                        set_N( copy, pt, pa[1+i]-1 );
                    } else {
                        if( locked( pa[1+i] ) || int_urn(0, smap[1+i]+smap[pa[1+i]]-1) < smap[1+i] ) {
                            set_N( copy, pt, i );
                        } else {
                            set_N( copy, pt, pa[1+i]-1 );
                        }
                    }
                } else if( copy[i] != 'N' ) {
                    set_N( copy, pt, i );
                }
                
                if( mt[1+i] ) {
                    set_N( copy, pt, mt[1+i]-1 );
                }

                if( i > 0 && pt[1+i] && pt[1+i-1]==0 && mt[1+i-1] ) {
                    j = i - 1;
                    copy[j] = 'N';
                    set_N( copy, pt, mt[1+j]-1 );
                    if( pa[1+i]==0 ) {
                        do {
                          if( --j < 0 ) break;
                          set_N( copy, pt, j );
                        } while( pt[1+j]==0 );
                    }
                }
                if( i < len-1 && pt[1+i] && pt[1+i+1]==0 && mt[1+i+1] ) {
                    j = i + 1;
                    copy[j] = 'N';
                    set_N( copy, pt, mt[1+j]-1 );
                    if( pa[1+i]==0 ) {
                        do {
                          if( ++j >= len ) break;
                          set_N( copy, pt, j );
                        } while( pt[1+j]==0 );
                    }
                }
                if( pt[1+i]>1 && pt[pt[1+i]-1]==0 && mt[pt[1+i]-1] ) {
                    j = pt[1+i]-1 - 1;
                    copy[j] = 'N';
                    set_N( copy, pt, mt[1+j]-1 );
                }
                if( pt[1+i] && pt[1+i]<len && pt[pt[1+i]+1]==0 && mt[pt[1+i]+1] ) {
                    j = pt[1+i]-1 + 1;
                    copy[j] = 'N';
                    set_N( copy, pt, mt[1+j]-1 );
                }
            }

            // make sure we haven't messed with the locked bases from the program input        
            for( i = 0; i < len; i++ ) {
                if( start[i] != 'N' ) copy[i] = start[i];
            }

        } while( strspn( copy, "AUGC" ) == strlen( copy ) );

        if( verbosity > 2 ) {
            printf( "C: %s\n", position );
            printf( "S: " );
            for( k = 1; k <= len; k++ ) printf( "%c", retry[k] ? 'X' : '.' );
            printf( "\n" );
            printf( "N: %s\n", copy );
            fflush( stdout );
        }

        free( retry );
        free( la );
        free( ma );
        free( pa );

        // if we seem to be chasing our own tail, reset completely
        if( strcmp( last_copy, copy ) == 0 ) {
            if( ++stuck > 10 ) {
                strcpy( copy, start );
                stuck = 0;
            }
        } else {
            strcpy( last_copy, copy );
            stuck = 0;
        }
    }

    free( shuffle );
    free( smap );
    free( jct );
    free( lt );
    free( mt );
    free( pt );

    if( bpd == 0 ) {
        printf( "NMC: %s %f\n", position, final );
        printf( "STR: %s %.2f %d %d\n", secstr, e, bpd, n_iter - iter );
    } else {
        printf( "NMC: %s\n", closest_seq );
        printf( "STR: %s %.2f %d %d\n", closest_struct, closest_fe, closest_bpd, n_iter - iter );
        printf( "CPY: %s\n", copy );
    }

    return 0;
}

*/