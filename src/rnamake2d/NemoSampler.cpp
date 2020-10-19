// Copyright (C) 2018 Fernando Portela <nando8888@gmail.com>
#include <rnamake2d/NemoSampler.h>

// set to 0 to prevent use of base pair distances in the scoring function
#define USE_BPD 1

// set to 0 to prevent use of free energy differences in the scoring function
#define USE_DDG 1

// set to 0 to disable heuristics in the sampling phase
#define USE_DOMAIN_KNOWLEDGE 1

namespace rnamake2d {

    int
    my_urn(int from, int to) {
        return (((int) (drand48() * (to - from + 1))) + from);
    }

    int
    pair_map(char a, char b) {
        if ((a ^ b ^ 'G' ^ 'C') == 0) return 3;
        if ((a ^ b ^ 'A' ^ 'U') == 0) return 2;
        if ((a ^ b ^ 'G' ^ 'U') == 0) return 1;
        return 0;
    }

    std::vector<short>
    wrapped_make_pair_table(  char * target )  {
       auto table = std::vector<short>(strlen(target));
       short* data = make_pair_table(target);

       for(int ii = 0; ii < strlen(target) ; ++ii) {
           table[ii] = data[ii];
       }
       free( data );
        return table;
    }

    std::vector<short>
    make_mismatch_table(std::vector<short> const & pt) {
        //short *mt = (short *) calloc(1 + pt[0], sizeof(short));
        auto mt = std::vector<short>( 1 + pt[0],0 );
        for (int j = 1; j <= pt[0]; j++) {
            if (pt[j]) {
                if (j < pt[0] && pt[j] > 1 && pt[j + 1] && pt[pt[j] - 1]
                    && pt[j + 1] != pt[j] - 1
                    ) {
                    mt[j + 1] = pt[j] - 1;
                    mt[pt[j] - 1] = j + 1;
                }
                continue;
            }
            if (j < pt[0] && pt[j + 1] && pt[j + 1] + 1 <= pt[0]) {
                mt[j] = pt[j + 1] + 1;
                mt[pt[j + 1] + 1] = j;
            } else if (j > 1 && pt[j - 1] && pt[j - 1] - 1 >= 1) {
                mt[j] = pt[j - 1] - 1;
                mt[pt[j - 1] - 1] = j;
            }
        }

        return mt;
    }

// Initialize helper arrays marking the closing pairs in various loops
//
    std::vector<int>
    scan_loops(std::vector<short> pt) {
        int len = pt[0];
        auto map = std::vector<int>(len+1) ;
        //int *map = (int *) calloc(1 + len, sizeof(int));
        int mark = 1;

        for (int j = 1; j < len; j++) {
            if (pt[j] <= 1) continue;
            if (map[j]) continue;
            if (pt[pt[j] - 1] != j + 1) {
                map[j] = mark;
                int k = j;
                do {
                    do {
                        if (++k > len) k = 1;
                    } while (pt[k] == 0);
                    k = pt[k];
                    map[k] = mark;
                } while (k != j);
                mark++;
            }
        }

        map[0] = mark;
        return map;
    }

    std::vector<int>
    scan_junctions(std::vector<short> const& pt, std::vector<int> const & lt) {
        int len = pt[0];
        //int *map = (int *) calloc(1 + len, sizeof(int));
        auto map = std::vector<int>(len + 1);
        //memmove(map, lt, (1 + len) * sizeof(int));

        for(int ii = 0; ii < lt.size(); ++ii) {
           map[ii]  = lt[ii];
        }

        int mark = lt[0];

        while (--mark) {
            int j;
            int count = 0;
            for (j = 1; j <= len; j++) if (map[j] == mark) count++;
            if (count <= 2) {
                for (j = 1; j <= len; j++) if (map[j] == mark) map[j] = 0;
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
    std::vector<int>
    make_strength_map(std::vector<short> const & pt  , std::vector<short> const & mt) {
        int len = pt[0];
        //int *map = (int *) calloc(1 + len, sizeof(int));
        auto map = std::vector<int>(1+len) ;
        int j;

        for (j = 1; j <= len; j++) {
            int i, k;
            if (pt[j]) {
                for (i = 1; i < 3; i++) {
                    k = j - i;
                    if (k >= 1) map[k] += 3 - i;
                    k = j + i;
                    if (k <= len) map[k] += 3 - i;
                }
            } else if (mt[j]) {
                int up = j;
                while (up <= len && pt[up] == 0) up++;
                int down = j;
                while (down >= 1 && pt[down] == 0) down--;

                int minus = 2;
                if (up > len || down < 1) { // dangling end
                    minus = 4;
                } else if (pt[down] > down && pt[up] < up) { // hairpin
                    minus = 3;
                } else if (pt[down] < down && pt[up] > up) { // multiloop
                    minus = 4;
                } else { // internal loop I guess...
                    if (pt[down] - pt[up] == 1) { // bulge
                        minus = up - down == 2 ? 5 : 4;
                    } else {
                        minus = up - down == 2 ? 3 : 2;
                    }
                }
                for (i = 1; i < minus; i++) {
                    k = j - i;
                    if (k >= 1) map[k] -= minus - i;
                    k = j + i;
                    if (k <= len) map[k] -= minus - i;
                }
            }
        }

        int m = 999;
        for (j = 1; j <= len; j++) if (map[j] < m) m = map[j];
        for (j = 1; j <= len; j++) map[j] += 1 - m;

        return map;
    }

    void set_N(char *seq, std::vector<short> const & pt, int pos) {
        seq[pos] = 'N';
        if (pt[1 + pos]) seq[pt[1 + pos] - 1] = 'N';
    }

    double normal_score(const char *position, const char *target, bool found, int npairs) {
        char *secstr = strdup(position);
        secstr[0] = '\0';
        double e = fold(position, secstr);
        int bpd = bp_distance(target, secstr);
#if USE_DDG
        double es = energy_of_structure(position, target, 0);
#endif
        if( secstr != nullptr) free(secstr);
        if (bpd == 0) {
            found = true;
            return 1.0;
        }
        double score;
#if USE_BPD
        score = (npairs == 0 ? 1.0 / (1.0 + bpd) : 1.0 - (0.5 * bpd) / npairs);
#endif
#if USE_DDG
        double e_factor = 1.01 + es - e;
#if USE_BPD
        score *= (score < 0 ? e_factor : 1.0 / e_factor);
#else
        score = 1.0/e_factor;
#endif
#endif
        return score;
    }

    const char *maxCC = "CCC";
    const char *maxGG = "GGG";

    bool is_legal(char *position) {
#if 0
        if( maxCC && strstr( position, maxCC ) != nullptr ) return false;
        if( maxGG && strstr( position, maxGG ) != nullptr ) return false;
#endif
        return true;
    }

    bool undecided(char x) {
        return (strchr("AUGC", x) == nullptr);
    }

    bool locked(char *start, int pos) {
        return !undecided(start[pos - 1]);
    }

    bool test_move(char *position, int o, char base) {
        bool ok = true;
        char save = position[o];
        position[o] = base;
        ok = is_legal(position);
        position[o] = save;
        return ok;
    }

    bool
    NemoSampler::play_move(char *position, const char *bases, int z = -1) {
        int l = strspn(position, bases);
        int len = pt[0];

        if (pt[l + 1]) { // it's a pair
            if (z < 0) { // and we're instructed to choose a random one
                int w[] = {16, 16, 9, 9, 2, 2};

#if USE_DOMAIN_KNOWLEDGE
                // Tuning of weights for adjacent stacks in junctions.
                // When looking at adjcent stacks in junctions from the point of view of an
                // observer at the center of the loop, it is usually better to make sure that
                // the left one has a high probability of being GC/CG, while the rightmost one
                // can usually afford to be demoted to AU/UA.
                //
                if (l < len - 1 && mt[1 + l + 1] && pt[1 + l + 1] && pt[1 + l + 1] < len &&
                    mt[pt[1 + l + 1] + 1] == 1 + l) {
                    // ok, looks like adjacent stacks, but we need to eliminate the possibility
                    // that it could be a 0-N bulge
                    if (jct[1 + l] || jct[pt[1 + l]]) {
                        w[0] -= 6;
                        w[1] -= 6;
                        w[2] += 6;
                        w[3] += 6;
                    }
                }
                if (l > 0 && mt[1 + l - 1] && pt[1 + l - 1] && pt[1 + l - 1] > 1 && mt[pt[1 + l - 1] - 1] == 1 + l) {
                    if (jct[1 + l] || jct[pt[1 + l]]) {
                        w[0] += 6;
                        w[1] += 6;
                        w[2] -= 6;
                        w[3] -= 6;
                    }
                }

                // Tuning for triloops
                // Only GC/CG closing pairs work if the supporting pair also is GC/CG
                int d, j, k;
                if (pt[1 + l] > (1 + l)) {
                    d = pt[1 + l] - (1 + l);
                    j = 1 + l + 1;
                    k = 1 + l - 1;
                } else {
                    d = (1 + l) - pt[1 + l];
                    j = pt[1 + l] + 1;
                    k = pt[1 + l] - 1;
                }
                if (d == 4 && pair_map(position[k - 1], position[pt[1 + k - 1] - 1]) == 3) {
                    w[2] = 0;
                    w[3] = 0;
                    w[4] = 0;
                    w[5] = 0;
                }
                if (d == 6 && pt[j] && pt[j] - j == 4 && !undecided(position[j - 1])
                    && pair_map(position[j - 1], position[pt[j] - 1]) != 3) {
                    w[0] = 0;
                    w[1] = 0;
                }
#endif

                if (!test_move(position, l, 'C')) w[0] = 0;
                if (!test_move(position, l, 'G')) {
                    w[1] = 0;
                    w[5] = 0;
                }
                if (!test_move(position, l, 'U')) {
                    w[2] = 0;
                    w[4] = 0;
                }
                if (!test_move(position, l, 'A')) w[3] = 0;
                if (!test_move(position, pt[1 + l] - 1, 'C')) w[1] = 0;
                if (!test_move(position, pt[1 + l] - 1, 'G')) {
                    w[0] = 0;
                    w[4] = 0;
                }
                if (!test_move(position, pt[1 + l] - 1, 'U')) {
                    w[3] = 0;
                    w[5] = 0;
                }
                if (!test_move(position, pt[1 + l] - 1, 'A')) w[2] = 0;

                // FIXME: what if all w[]==0 ?
                int dice = int_urn(0, w[0] + w[1] + w[2] + w[3] + w[4] + w[5] - 1);
                for (z = 0; dice >= w[z]; z++) dice -= w[z];
                // should leave us with 0 <= z <= 5 randomly selected according to weights
                // FIXME: isn't this a perfect spot for asserting?
            }
            switch (z) {
                case 0:
                    position[l] = 'C';
                    position[pt[l + 1] - 1] = 'G';
                    break;
                case 1:
                    position[l] = 'G';
                    position[pt[l + 1] - 1] = 'C';
                    break;
                case 2:
                    position[l] = 'U';
                    position[pt[l + 1] - 1] = 'A';
                    break;
                case 3:
                    position[l] = 'A';
                    position[pt[l + 1] - 1] = 'U';
                    break;
                case 4:
                    position[l] = 'U';
                    position[pt[l + 1] - 1] = 'G';
                    break;
                case 5:
                    position[l] = 'G';
                    position[pt[l + 1] - 1] = 'U';
                    break;
                default:
                    return false;
            }
        } else { // unpaired, so choose a random base
            int w[] = {93, 1, 5, 1};
            if (z < 0 && mt[l + 1]) { // is mismatched
                if (pt[mt[l + 1]]) { // mismatch is paired
                    // FIXME: this would probably be a good place to use some (ML-based?) optimization
                    //        for weights fine-tuning
                    switch (position[mt[l + 1] - 1]) {
                        case 'A':
                            w[0] = 5;
                            w[1] = 0;
                            w[2] = 2;
                            w[3] = 1;
                            break;
                        case 'U':
                            w[0] = 0;
                            w[1] = 6;
                            w[2] = 1;
                            w[3] = 4;
                            break;
                        case 'G':
                            w[0] = 2;
                            w[1] = 1;
                            w[2] = 5;
                            w[3] = 0;
                            break;
                        case 'C':
                            w[0] = 6;
                            w[1] = 4;
                            w[2] = 0;
                            w[3] = 1;
                            break;
                    }
                } else { // mismatch is unpaired
                    int dice;
                    int up, down, close;
                    if (l > 0 && pt[1 + l - 1] == mt[1 + l] + 1) {
                        up = l + 1;
                        down = mt[1 + l];
                    } else {
                        down = l + 1;
                        up = mt[1 + l];
                    }
                    close = up - 1;
                    int u, d;
                    for (u = 0; up + u <= len && pt[up + u] == 0; u++);
                    for (d = 0; down > d && pt[down - d] == 0; d++);
                    bool internal_loop = (up + u <= len && down > d && pt[up + u] == down - d);
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
                    if (internal_loop && down - d == close) { // unbranched & single segment = hairpin
#if USE_DOMAIN_KNOWLEDGE
                        // A test for a potential 'slide' from a triloop to a GAAA tetraloop
                        // If matched, have a good chance to apply the "anti-boost" (U/C in the middle
                        // of the triloop)
                        if (l > 1 && d == 3 && undecided(position[l + 1]) && position[l - 1] == 'G'
                            && pair_map(position[l - 2], position[l + 3]) > 0
                            && (force_boost || int_urn(0, 1) == 0)) {
                            dice = int_urn(0, 9);
                            if (dice < 5) {
                                position[l + 1] = 'U';
                            } else if (dice < 9) {
                                position[l + 1] = 'C';
                            }
                        }
                        // Simple apical loop G/A boosting
                        if (undecided(position[mt[1 + l] - 1]) && d > 3
                            && (force_boost || int_urn(0, 5) != 0)) {
                            if (test_move(position, l, 'G')
                                && test_move(position, mt[1 + l] - 1, 'A')) {
                                z = 2;
                                position[mt[1 + l] - 1] = 'A';
                            }
                        }
#endif
                    } else if (internal_loop) {
                        if (mt[l + 1] > l + 1) {
                            w[0] = 5;
                            w[1] = 1;
                            w[2] = 20;
                            w[3] = 1;
                        } else {
                            w[0] = 21;
                            w[1] = 1;
                            w[2] = 4;
                            w[3] = 1;
                        }
                        if (undecided(position[mt[1 + l] - 1])) {
#if USE_DOMAIN_KNOWLEDGE
                            // 1-1 internal loop
                            if (u == 1 && d == 1) {
                                dice = int_urn(0, 9);
                                if (force_boost || dice < 8) {
                                    if (test_move(position, l, 'G')
                                        && test_move(position, mt[1 + l] - 1, 'G')) {
                                        z = 2;
                                        position[mt[1 + l] - 1] = 'G';
                                    }
                                }
                            }
                            // 2-2 internal loop
                            if (u == 2 && d == 2) {
                                dice = int_urn(0, 1);
                                if ((force_boost || dice == 0) && undecided(position[l + 1])
                                    && undecided(position[mt[1 + l] - 1 - 1])) {
                                    if (test_move(position, l, 'U')
                                        && test_move(position, mt[1 + l] - 1, 'G')
                                        && test_move(position, l + 1, 'G')
                                        && test_move(position, mt[1 + l] - 1 - 1, 'U')) {
                                        z = 1;
                                        position[mt[1 + l] - 1] = 'G';
                                        position[l + 1] = 'G';
                                        position[mt[1 + l] - 1 - 1] = 'U';
                                    }
                                }
                            }
                            // if a selection hasn't been made yet, try typical boosts for
                            // internal loops (G/A, A/G, U/U)
                            if (z < 0) {
                                dice = int_urn(0, 9);
                                if (dice < 3) {
                                    z = 2;
                                    position[mt[1 + l] - 1] = 'A';
                                } else if (dice < 6) {
                                    z = 0;
                                    position[mt[1 + l] - 1] = 'G';
                                } else if (dice < 7) {
                                    z = 1;
                                    position[mt[1 + l] - 1] = 'U';
                                }
                            }
#endif
                        } else {
                            switch (position[mt[1 + l] - 1]) {
                                case 'A':
                                    w[0] = 4;
                                    w[1] = 0;
                                    w[2] = 4;
                                    w[3] = 1;
                                    break;
                                case 'U':
                                    w[0] = 0;
                                    w[1] = 6;
                                    w[2] = 1;
                                    w[3] = 2;
                                    break;
                                case 'G':
                                    w[0] = 6;
                                    w[1] = 1;
                                    w[2] = 2;
                                    w[3] = 0;
                                    break;
                                case 'C':
                                    w[0] = 4;
                                    w[1] = 1;
                                    w[2] = 0;
                                    w[3] = 1;
                                    break;
                            }
                        }
                    } else { // a mismatch in a junction or external loop
                        w[0] = 97;
                        w[1] = 1;
                        w[2] = 1;
                        w[3] = 1;
#if USE_DOMAIN_KNOWLEDGE
                        if (down == 1 + l) {
                            if ((l < 1 || pt[1 + l - 1] == 0) && position[l + 1] == 'G' &&
                                position[pt[1 + l + 1] - 1] == 'C') {
                                w[3] = 48; // increase chance of C boost
                            }
                        } else {
                            if ((l >= len - 1 || pt[1 + l + 1] == 0) && position[l - 1] == 'G' &&
                                position[pt[1 + l - 1] - 1] == 'C') {
                                w[2] = 48; // increase chance of G boost
                            }
                        }
#endif
                    }
                }
            }

            if (z < 0) {
                if (!test_move(position, l, 'A')) w[0] = 0;
                if (!test_move(position, l, 'U')) w[1] = 0;
                if (!test_move(position, l, 'G')) w[2] = 0;
                if (!test_move(position, l, 'C')) w[3] = 0;
                int dice = int_urn(0, w[0] + w[1] + w[2] + w[3] - 1);
                for (z = 0; dice >= w[z]; z++) dice -= w[z];
            }
            switch (z) {
                case 0:
                    position[l] = 'A';
                    break;
                case 1:
                    position[l] = 'U';
                    break;
                case 2:
                    position[l] = 'G';
                    break;
                case 3:
                    position[l] = 'C';
                    break;
                default: break;
            }

        }
        return is_legal(position);
    }

    double
    NemoSampler::sample( char *position )  {
        int k;
        for (k = 0; k < pt[0]; k++) {
            if (position[k] == 'N' && pt[1 + k]) position[k] = 'P';
        }
        // first fill the target pairs
        char bases[] = "AUGCN";
        while (strspn(position, bases) != strlen(position)) {
            play_move(position, bases, -1);
        }
        // now fill the rest, i.e. the unpaired bases
        bases[4] = 0;
        while (strspn(position, bases) != strlen(position)) {
            play_move(position, bases,  -1);
        }
        return normal_score(position, target, found, npairs);
    }


    double
    NemoSampler::nested( char *position , int level ) {
        if (strspn(position, "AUGC") == strlen(position)) {
            return normal_score(position, target, found, npairs);
        }

        double best = worst_score;
        char *best_playout = strdup(position);
        char *best_local = strdup(position);

        while (strspn(position, "AUGC") != strlen(position)) {
            int l = strspn(position, "AUGC");
            double max = worst_score;
            int z;
            int zm = pt[l + 1] == 0 ? 4 : 6;
            if (level == 1) {
#pragma omp parallel for schedule(dynamic)
                for (z = 0; z < zm; z++) {
                    char *playout = strdup(position);
                    if (!play_move(playout, "AUGC",  z)) continue;

                    double v = sample( playout );
#pragma omp critical(max_update)
                    {
                        if (v > max) {
                            max = v;
                            strcpy(best_local, playout);
                        }
                    }
                    if( playout != nullptr) free(playout);
                }
            } else {
                for (z = 0; z < zm; z++) {
                    char *playout = strdup(position);
                    if (!play_move(playout, "AUGC",  z)) continue;
                    double v = nested(playout, level - 1 );

                    if (v > max) {
                        max = v;
                        strcpy(best_local, playout);
                    }
                    if( playout != nullptr) free(playout);
                }
            }

            if (standard_nmcs || (max > best)) {
                best = max;
                strcpy(best_playout, best_local);
            }

            if (found) {
                strcpy(position, best_playout);
                break;
            }

            position[l] = best_playout[l];
            if (pt[l + 1]) {
                position[pt[l + 1] - 1] = best_playout[pt[l + 1] - 1];
            } else if (mt[l + 1] && position[mt[l + 1] - 1] == 'N' && pt[mt[l + 1] - 1] == 0) {
                position[mt[l + 1] - 1] = best_playout[mt[l + 1] - 1];
            }
        }

        if( best_playout != nullptr ) free(best_playout);
        if( best_local != nullptr) free(best_local);
        return best;
    }

    int
    NemoSampler::nemo_main(Design& design) {

        char *copy = strdup(start);
        char *secstr = strdup(target);
        secstr[0] = '\0';

        int closest_bpd = 999999;
        char *closest_struct = strdup(target);
        double closest_fe;

        int stuck = 0;
        char *last_copy = strdup(start);

        double final, e;
        int n_iter = iter;

        if (seed) strcpy(copy, seed);

        for (; iter; iter--) {
            strcpy(position, copy);
            // half the time, use the best boosts we know of during playouts
            force_boost = int_urn(0, 1) == 0;
            // (try to) solve
            final = nested( position,  1 );
            // score the search result
            secstr[0] = 0;
            e = fold(position, secstr);
            bpd = bp_distance(target, secstr);
            if (bpd == 0) break;

            if (bpd < closest_bpd) {
                closest_bpd = bpd;
                strcpy(closest_seq, position);
                strcpy(closest_struct, secstr);
                closest_fe = e;
            }

            // Our previous attempt failed. Rather than simply increase depth, or just
            // try again with the original starting point, we build a new starting point based
            // on what seems to have worked and what didn't (misfolds)

            int i, j, k, c;

            // let's get a pair map of the misfolded structure
            auto pa = wrapped_make_pair_table(secstr);

            bool *retry = (bool *) calloc(1 + len, sizeof(bool));
            for (j = 1; j <= len; j++) retry[j] = (pt[j] != pa[j]);

            for (j = 1; j <= len; j++) {
                // add mismatches of the misfolded bases
                if (!retry[j] && mt[j] && retry[mt[j]]) retry[j] = true;
            }

            auto ma = make_mismatch_table(pa);
            auto la = scan_loops(pa);

            for (j = 1; j < la[0]; j++) {
                for (i = 1; i <= len && la[i] != j; i++);
                bool has_opened_pair = false;
                bool closing_good = true;
                k = i;
                do {
                    do {
                        if (++k > len) k = 1;
                        if (pa[k] == 0 && pt[k]) has_opened_pair = true;
                    } while (pa[k] == 0);
                    if (pa[k] != pt[k]) closing_good = false;
                    k = pa[k];
                } while (k != i && closing_good);

                if (closing_good && has_opened_pair) {
                    do {
                        do {
                            if (++k > len) k = 1;
                            if (pa[k] || ma[k]) retry[k] = true;
                        } while (pa[k] == 0);
                        k = pa[k];
                        retry[k] = true;
                    } while (k != i);
                }
            }

            // we're ready

            do {

                strcpy(copy, position);

                for (k = 0; k < len; k++) {
                    // choose random r, with r != k
                    int r = int_urn(0, len - 2);
                    if (r >= k) r++;
                    shuffle[k] ^= shuffle[r];
                    shuffle[r] ^= shuffle[k];
                    shuffle[k] ^= shuffle[r];
                }

                for (k = 0, c = 0; k < len; k++) {
                    i = shuffle[k];

                    // not a misfolded (or otherwise interesting) spot? let's keep it
                    if (!retry[1 + i]) continue;

                    // we don't want to reset and retry all misfolded bases and pairs
                    // so we use probabilities, first 1/1, then 1/2, 1/3 and so on
                    // this guarantess that we will reset at least one base or pair,
                    // and maybe a few more.
                    if (int_urn(0, c++) != 0) continue;

                    // FIXME: Ad hoc rules, based on personal experience.
                    // Could probably use some ML improvements as well...

                    if (copy[i] != 'N' && pa[1 + i] && copy[pa[1 + i] - 1] != 'N') {
                        if (pt[1 + i] == 0 && pt[pa[1 + i]] == 0) {
                            set_N(copy, pt, i);
                            set_N(copy, pt, pa[1 + i] - 1);
                        } else {
                            if (locked(start, pa[1 + i]) ||
                                int_urn(0, smap[1 + i] + smap[pa[1 + i]] - 1) < smap[1 + i]) {
                                set_N(copy, pt, i);
                            } else {
                                set_N(copy, pt, pa[1 + i] - 1);
                            }
                        }
                    } else if (copy[i] != 'N') {
                        set_N(copy, pt, i);
                    }

                    if (mt[1 + i]) {
                        set_N(copy, pt, mt[1 + i] - 1);
                    }

                    if (i > 0 && pt[1 + i] && pt[1 + i - 1] == 0 && mt[1 + i - 1]) {
                        j = i - 1;
                        copy[j] = 'N';
                        set_N(copy, pt, mt[1 + j] - 1);
                        if (pa[1 + i] == 0) {
                            do {
                                if (--j < 0) break;
                                set_N(copy, pt, j);
                            } while (pt[1 + j] == 0);
                        }
                    }
                    if (i < len - 1 && pt[1 + i] && pt[1 + i + 1] == 0 && mt[1 + i + 1]) {
                        j = i + 1;
                        copy[j] = 'N';
                        set_N(copy, pt, mt[1 + j] - 1);
                        if (pa[1 + i] == 0) {
                            do {
                                if (++j >= len) break;
                                set_N(copy, pt, j);
                            } while (pt[1 + j] == 0);
                        }
                    }
                    if (pt[1 + i] > 1 && pt[pt[1 + i] - 1] == 0 && mt[pt[1 + i] - 1]) {
                        j = pt[1 + i] - 1 - 1;
                        copy[j] = 'N';
                        set_N(copy, pt, mt[1 + j] - 1);
                    }
                    if (pt[1 + i] && pt[1 + i] < len && pt[pt[1 + i] + 1] == 0 && mt[pt[1 + i] + 1]) {
                        j = pt[1 + i] - 1 + 1;
                        copy[j] = 'N';
                        set_N(copy, pt, mt[1 + j] - 1);
                    }
                }

                // make sure we haven't messed with the locked bases from the program input
                for (i = 0; i < len; i++) {
                    if (start[i] != 'N') copy[i] = start[i];
                }

            } while (strspn(copy, "AUGC") == strlen(copy));

            if( retry != nullptr ) free(retry);

            // if we seem to be chasing our own tail, reset completely
            if (strcmp(last_copy, copy) == 0) {
                if (++stuck > 10) {
                    strcpy(copy, start);
                    stuck = 0;
                }
            } else {
                strcpy(last_copy, copy);
                stuck = 0;
            }
        }

        if (bpd == 0) {
            design.sequence.reserve(strlen(position));
            strncpy(&(design.sequence[0]),position,strlen(position));

            //printf("NMC: %s %f\n", position, final);
            //printf("STR: %s %.2f %d %d\n", secstr, e, bpd, n_iter - iter);
        } else {
            design.sequence.reserve(strlen(closest_seq));
            strncpy(&(design.sequence[0]),closest_seq,strlen(closest_seq));

            //printf("NMC: %s\n", closest_seq);
            //printf("STR: %s %.2f %d %d\n", closest_struct, closest_fe, closest_bpd, n_iter - iter);
            //printf("CPY: %s\n", copy);
        }
        return 0;
    }

    void
    NemoSampler::initialize_design(Design& design) {

        if(current_) {
            reset_();
        }

        current_ = design.id;
        const int size = design.target.size();
        target = (char*)malloc(sizeof(char)*(size+1)); target[size] = '\0';
        start = (char*)malloc(sizeof(char)*(size+1)); start[size] = '\0';
        seed = (char*)malloc(sizeof(char)*(size+1)); seed[size] = '\0';

        position = (char*)malloc(sizeof(char)*(size+1)); position[size] = '\0';
        closest_seq = (char*)malloc(sizeof(char)*(size+1)); closest_seq[size] = '\0';

        strncpy( closest_seq, start, strlen( start ) );
        strncpy( position, start,  strlen( start ) );

        strncpy( target, design.target.c_str(), design.target.size() );
        strncpy( start, design.sequence.c_str(), design.sequence.size() );
        strncpy( seed, design.sequence.c_str(), design.sequence.size() );

        len = strlen( seed );
        pt = wrapped_make_pair_table( target );
        mt = make_mismatch_table( pt );
        lt = scan_loops( pt );
        jct = scan_junctions( pt , lt );
        smap = make_strength_map( pt , mt );

        char* p;
        for( npairs = 0, p = target; p[npairs]; p[npairs]=='('? (void)npairs++ : (void)p++ );

        shuffle.reserve(len);

        for(int  k = 0; k < len; k++ ) shuffle[k] = k;
        nemo_main(design);
        design.initialize_features();
    }

    void
    NemoSampler::mutate(Design & design)  {
        if( current_ != design.id ) {
            std::cerr<<"Design id is NOT equal to the current Sampler id Number. Exiting...\n";
            exit(0);
        }
        char *copy = strdup(start);
        char *secstr = strdup(target);
        secstr[0] = '\0';

        int closest_bpd = 999999;
        char *closest_struct = strdup(target);
        double closest_fe;

        int stuck = 0;
        char *last_copy = strdup(start);

        double final, e;
        int n_iter = iter;

        if (seed) strcpy(copy, seed);

       strcpy(position, copy);
       // half the time, use the best boosts we know of during playouts
       force_boost = int_urn(0, 1) == 0;
       // (try to) solve
       final = nested( position,  1 );
       // score the search result
       secstr[0] = 0;
       e = fold(position, secstr);
       bpd = bp_distance(target, secstr);

       design.candiate = position;
       design.initialize_mutant();
       //design.update(true);
       return ;
       if (bpd < closest_bpd) {
           closest_bpd = bpd;
           strcpy(closest_seq, position);
           strcpy(closest_struct, secstr);
           closest_fe = e;
       }

       int i, j, k, c;

       // let's get a pair map of the misfolded structure
       auto pa = wrapped_make_pair_table(secstr);

       bool *retry = (bool *) calloc(1 + len, sizeof(bool));
       for (j = 1; j <= len; j++) retry[j] = (pt[j] != pa[j]);

       for (j = 1; j <= len; j++) {
           // add mismatches of the misfolded bases
           if (!retry[j] && mt[j] && retry[mt[j]]) retry[j] = true;
       }

       auto ma = make_mismatch_table(pa);
       auto la = scan_loops(pa);

       for (j = 1; j < la[0]; j++) {
           for (i = 1; i <= len && la[i] != j; i++);
           bool has_opened_pair = false;
           bool closing_good = true;
           k = i;
           do {
               do {
                   if (++k > len) k = 1;
                   if (pa[k] == 0 && pt[k]) has_opened_pair = true;
               } while (pa[k] == 0);
               if (pa[k] != pt[k]) closing_good = false;
               k = pa[k];
           } while (k != i && closing_good);

           if (closing_good && has_opened_pair) {
               do {
                   do {
                       if (++k > len) k = 1;
                       if (pa[k] || ma[k]) retry[k] = true;
                   } while (pa[k] == 0);
                   k = pa[k];
                   retry[k] = true;
               } while (k != i);
           }
       }

       do {

           strcpy(copy, position);

           for (k = 0; k < len; k++) {
               // choose random r, with r != k
               int r = int_urn(0, len - 2);
               if (r >= k) r++;
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

           for (k = 0, c = 0; k < len; k++) {
               i = shuffle[k];

               // not a misfolded (or otherwise interesting) spot? let's keep it
               if (!retry[1 + i]) continue;

               // we don't want to reset and retry all misfolded bases and pairs
               // so we use probabilities, first 1/1, then 1/2, 1/3 and so on
               // this guarantess that we will reset at least one base or pair,
               // and maybe a few more.
               if (int_urn(0, c++) != 0) continue;

               // FIXME: Ad hoc rules, based on personal experience.
               // Could probably use some ML improvements as well...

               if (copy[i] != 'N' && pa[1 + i] && copy[pa[1 + i] - 1] != 'N') {
                   if (pt[1 + i] == 0 && pt[pa[1 + i]] == 0) {
                       set_N(copy, pt, i);
                       set_N(copy, pt, pa[1 + i] - 1);
                   } else {
                       if (locked(start, pa[1 + i]) ||
                           int_urn(0, smap[1 + i] + smap[pa[1 + i]] - 1) < smap[1 + i]) {
                           set_N(copy, pt, i);
                       } else {
                           set_N(copy, pt, pa[1 + i] - 1);
                       }
                   }
               } else if (copy[i] != 'N') {
                   set_N(copy, pt, i);
               }

               if (mt[1 + i]) {
                   set_N(copy, pt, mt[1 + i] - 1);
               }

               if (i > 0 && pt[1 + i] && pt[1 + i - 1] == 0 && mt[1 + i - 1]) {
                   j = i - 1;
                   copy[j] = 'N';
                   set_N(copy, pt, mt[1 + j] - 1);
                   if (pa[1 + i] == 0) {
                       do {
                           if (--j < 0) break;
                           set_N(copy, pt, j);
                       } while (pt[1 + j] == 0);
                   }
               }
               if (i < len - 1 && pt[1 + i] && pt[1 + i + 1] == 0 && mt[1 + i + 1]) {
                   j = i + 1;
                   copy[j] = 'N';
                   set_N(copy, pt, mt[1 + j] - 1);
                   if (pa[1 + i] == 0) {
                       do {
                           if (++j >= len) break;
                           set_N(copy, pt, j);
                       } while (pt[1 + j] == 0);
                   }
               }
               if (pt[1 + i] > 1 && pt[pt[1 + i] - 1] == 0 && mt[pt[1 + i] - 1]) {
                   j = pt[1 + i] - 1 - 1;
                   copy[j] = 'N';
                   set_N(copy, pt, mt[1 + j] - 1);
               }
               if (pt[1 + i] && pt[1 + i] < len && pt[pt[1 + i] + 1] == 0 && mt[pt[1 + i] + 1]) {
                   j = pt[1 + i] - 1 + 1;
                   copy[j] = 'N';
                   set_N(copy, pt, mt[1 + j] - 1);
               }
           }

           // make sure we haven't messed with the locked bases from the program input
           for (i = 0; i < len; i++) {
               if (start[i] != 'N') copy[i] = start[i];
           }

       } while (strspn(copy, "AUGC") == strlen(copy));


       if( retry != nullptr ) free(retry);

       // if we seem to be chasing our own tail, reset completely
       if (strcmp(last_copy, copy) == 0) {
           if (++stuck > 10) {
               strcpy(copy, start);
               stuck = 0;
           }
       } else {
           strcpy(last_copy, copy);
           stuck = 0;
       }

        if (bpd == 0 && !strcmp(design.sequence.c_str() ,position )) {
            //design.candiate = position;
            design.candiate.reserve( strlen(position) ) ;
            strncpy( &(design.candiate[0]), position, strlen(position));
            //printf("NMC: %s %f\n", position, final);
            //printf("STR: %s %.2f %d %d\n", secstr, e, bpd, n_iter - iter);
        } else if (!strcmp(design.sequence.c_str(),closest_seq)) {
            //design.candiate = closest_seq;
            design.candiate.reserve( strlen(closest_seq) ) ;
            strncpy( &(design.candiate[0]), closest_seq, strlen(closest_seq));            //printf("NMC: %s\n", closest_seq);
            //printf("STR: %s %.2f %d %d\n", closest_struct, closest_fe, closest_bpd, n_iter - iter);
            //printf("CPY: %s\n", copy);
        } else {
            std::cout<<"Warning: no mutation occurred for desin: "<<current_<<"\n";
        }
    }
}