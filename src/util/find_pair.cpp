#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wwritable-strings"
#pragma GCC diagnostic ignored "-Wformat-security"


#include <util/x3dna_src.h>
#include <util/find_pair.h>
#include <util/x3dna.h>
#include <math/numerical.h>
#include <memory>
#include <string>
#include <stdexcept>

using namespace std;

namespace util {
    PairFinder::PairFinder(string pdb) {
        char *c = &pdb[0];

        sprintf(args.pdbfile, c);
        sprintf(args.outfile, "stdout");
        
        X3dna::X3BPInfo *bp_info;
    }

    /* clean up files with fixed name */
    static void clean_files(void) {
        remove_file(MUL_FILE);
        remove_file(ALLP_FILE);
        remove_file(BPORDER_FILE);
        remove_file(BESTP_FILE);
        remove_file(REF_FILE);
        remove_file(MREF_FILE);
        remove_file(HLXREG_FILE);
        remove_file(COLCHN_FILE);
        remove_file(COLHLX_FILE);
        remove_file(MULBP_FILE);
        remove_file(TMP_FILE);
    }

    template <typename... Args>
    std::string string_format(const std::string &format, Args... args) {
        // Extra space for '\0'
        size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1;
        if (size <= 0) {
            throw std::runtime_error("Error during formatting.");
        }
        std::unique_ptr<char[]> buf(new char[size]);
        snprintf(buf.get(), size, format.c_str(), args...);
        return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
    }

    /* a user reported confusion regarding the directionality of the reference
            * frame x-, y- and z-axis: row-wise or column-wise? Here adding more info to
            * make this point clear: it is ROW-wise. Keven (U. Penn); Oct 9, 2007 */
    void 
    PairFinder::_write_fpmst(double *morg, double *morien, FILE *rframe, X3dna::X3BPInfo *bp_info) {

        int i, j;
        // frprintf(rframe, "%10.4f %10.4f %10.4f  # origin\n",
        //         morg[1], morg[2], morg[3]);

        auto d = math::Point(morg[1], morg[2], morg[3]);
        auto rs = math::Points();
        for (i = 1; i <= 3; i++) {
            j = (i - 1) * 3;
            // frprintf(rframe, "%10.4f %10.4f %10.4f  # %c-axis\n", morien[j + 1],
            //         morien[j + 2], morien[j + 3],
            //         (i == 1) ? 'x' : (i == 2) ? 'y': 'z');
            string num_1 = string_format("%10.4f", morien[j + 1]);
            string num_2 = string_format("%10.4f", morien[j + 2]);
            string num_3 = string_format("%10.4f", morien[j + 3]);
            auto p = math::Point(stod(num_1), stod(num_2), stod(num_3));
            rs.push_back(p);
        }
        auto reg = std::regex(
            "#\\s+(?:\\.+\\d+\\>)*(\\w+):\\.*(-*\\d+)\\S:\\[\\.*(\\S+)\\](\\w+)\\s+\\-\\s+(?:\\.+\\d+\\>)*(\\w+):\\.*(-*\\d+)\\S:\\[\\.*(\\S+)\\](\\w+)");
        auto r = math::Matrix(rs[0].x(), rs[0].y(), rs[0].z(),
                                  rs[1].x(), rs[1].y(), rs[1].z(),
                                  rs[2].x(), rs[2].y(), rs[2].z());
        auto res1 = X3dna::X3Residue{bp_info->res1_num, bp_info->res1_chain_id, ' '};
        auto res2 = X3dna::X3Residue{bp_info->res2_num, bp_info->res2_chain_id, ' '};
        auto bp = X3dna::X3Basepair{res1, res2, d, r, X3dnaBPType::cWUW};
        bps.push_back(bp);
    }

    /* print out selected list for checking: temporary */
    static void print_list(long num_residue, long **mylist, char *my_str, FILE *fp)
    {
        long i, j, numb;

        // frprintf(fp, "====================== %s ======================\n", my_str);
        for (i = 1; i <= num_residue; i++) {
            numb = mylist[i][0];
            if (!numb)
                continue;
            // frprintf(fp, "%4ld[%+2ld]:", i, numb);
            // for (j = 1; j <= labs(numb); j++)
            //     // frprintf(fp, " %5ld", mylist[i][j]);
            // frprintf(fp, "\n");
        }
        // frprintf(fp, "\n");
    }

    /* check if two lists are identical: same # & residue order */
    static int isequal_list(int *listA, int *listB) {
        int i, is_equal = 0;

        if (listA[0] == listB[0]) {
            for (i = 1; i <= listA[0]; i++)
                if (listA[i] != listB[i])
                    break;
            if (i > listA[0])
                is_equal = 1;
        }
        return is_equal;
    }

    /* find the best-paired residue id#
 * pair_stat[PSTNUM]:
 *     j, bpid, d, dv, angle, dNN, dsum, bp-org,  x1,   y1,   z1,   x2,   y2,   z2
 *col# 1   2    3   4    5     6     7     8-10  11-13 14-16 17-19 20-22 23-25 26-28
 */
    static void best_pair(long i, long num_residue, long *RY, long **seidx, double **xyz,
                          long *idx, double **NC1xyz, long *matched_idx, double **orien,
                          double **org, long **ring_atom, char **AtomName, char *bseq,
                          miscPars *misc_pars, long *pair_stat) {
        double ddmin = XBIG, rtn_val[RTNNUM];
        long bpid, j, k, nout = PSTNUM - 1;

        init_lvector(pair_stat, 1, nout, 0);

        for (j = 1; j <= num_residue; j++) {
            if (j == i || RY[j] < 0 || matched_idx[j])
                continue;
            check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx, AtomName, misc_pars,
                       rtn_val, &bpid, ring_atom, 0);
            if (bpid && rtn_val[5] < ddmin) {
                ddmin = rtn_val[5];
                pair_stat[1] = j;
                pair_stat[2] = bpid;
                for (k = 1; k <= nout - 2; k++)
                    pair_stat[2 + k] = lround(MFACTOR * rtn_val[k]);
            }
        }
    }

    /* check if 2 bps co-planar based on simple geometry criterion */
    static long bp_coplanar(long i, double d, double d2, double *txyz, double *txyz2, long n,
                            long *ddidx, double **bp_xyz) {
        double *dp = NULL; /* use pointer for efficiency */
        long j = 0, co_planar = 0;

        if (fabs(d) < OLCRT) {
            j = ddidx[1];
            dp = txyz;
        }
        else if (fabs(d2) < OLCRT) {
            j = ddidx[n];
            dp = txyz2;
        }
        if (j && fabs(dot(&bp_xyz[i][9], dp)) < OLCRT &&
            fabs(dot(&bp_xyz[i][18], dp)) < OLCRT &&
            fabs(dot(&bp_xyz[j][9], dp)) < OLCRT && fabs(dot(&bp_xyz[j][18], dp)) < OLCRT)
            co_planar = 1;
        return co_planar;
    }

    static long is_circular(long num_bp, long **bp_order) {
        long i;

        if (num_bp <= 2)
            return false;

        for (i = 1; i <= num_bp; i++)
            if (bp_order[i][1] != -1)
                return false;

        return true;
    }

    /* find base-pair neighbors using simple geometric criterion for re_ordering */
    static void bp_context(long num_bp, miscPars *misc_pars, double **bp_xyz,
                           long **bp_order, long **end_list, long *num_ends, FILE *tfp) {
        double helix_break = misc_pars->helix_break;
        double d = EMPTY_NUMBER, d2, d3, ddmin[9], txyz[4], txyz2[4], txyz3[4], zave[4];
        long i, j, k, m, n, overlap = 0, quadruple = 0, cnum = 8, ddidx[9];

        // frprintf(tfp, "\nBase-pair context information\n");

        for (i = 1; i <= num_bp; i++) {
            init_dvector(ddmin, 1, cnum, XBIG);
            init_lvector(ddidx, 1, cnum, 0);
            d = dot(&bp_xyz[i][9], &bp_xyz[i][18]); /* between base normals */
            (d <= 0.0) ? ddxyz(bp_xyz[i] + 18, bp_xyz[i] + 9, zave) : sumxyz(bp_xyz[i] + 18, bp_xyz[i] + 9, zave);
            vec_norm(zave);
            for (j = 1; j <= num_bp; j++) {
                if (j == i)
                    continue;
                ddxyz(bp_xyz[j], bp_xyz[i], txyz);
                d = veclen(txyz);
                for (k = 1; k <= cnum; k++)
                    if (d < ddmin[k]) {
                        for (m = cnum; m > k; m--)
                            if (ddidx[n = m - 1]) {
                                ddmin[m] = ddmin[n];
                                ddidx[m] = ddidx[n];
                            }
                        ddmin[k] = d;
                        ddidx[k] = j;
                        break;
                    }
            }
            if (ddidx[1] && ddidx[2]) {
            /* at least 2 nearest neighbors */
                if (ddmin[1] > helix_break) {
                    /* isolated bp */
                    end_list[++*num_ends][1] = i; /* [i 0 0] */
                    n = 2;
                } else {
                    if (!overlap && ddmin[1] < OLCRT)
                        overlap = 1;
                    ddxyz(bp_xyz[ddidx[1]], bp_xyz[i], txyz); /* i's nearest neighbor */
                    d = dot(zave, txyz);

                    /* check for swapped bps in deformed pdr026/ur0009 */
                    if (ddidx[3] && ddmin[2] <= helix_break && ddmin[3] <= helix_break) {
                        ddxyz(bp_xyz[ddidx[2]], bp_xyz[i], txyz2);
                        ddxyz(bp_xyz[ddidx[3]], bp_xyz[i], txyz3);
                        d2 = dot(zave, txyz2);
                        d3 = dot(zave, txyz3);
                        if (d * d2 < 0.0 && d * d3 < 0.0 && fabs(d2) > fabs(d3)) {
                            lval_swap(&ddidx[2], &ddidx[3]);
                            dval_swap(&ddmin[2], &ddmin[3]);
                            // // frprintf(stderr, "[swapping 2nd & 3rd] %4ld %8.2f %8.2f "
                            //                 "%8.2f %8.2f\n",
                            //         i, ddmin[2], ddmin[3], d2, d3);
                        }
                    }
                    n = 0;
                    for (j = 2; j <= cnum && ddidx[j]; j++) {
                        if (ddmin[j] > helix_break)
                            break;
                        ddxyz(bp_xyz[ddidx[j]], bp_xyz[i], txyz2);
                        d2 = dot(zave, txyz2);
                        if (d * d2 < 0.0) {
                            /* vertical direction */
                            n = j;
                            bp_order[i][1] = -1; /* middle base-pair */
                            bp_order[i][2] = ddidx[1];
                            bp_order[i][3] = ddidx[j];
                            break;
                        }
                    }
                    if (!n) { 
                        /* terminal bp */
                        n = 2;
                        end_list[++*num_ends][1] = i;
                        end_list[*num_ends][2] = ddidx[1];
                        bp_order[i][2] = ddidx[1];
                        ddxyz(bp_xyz[ddidx[1]], bp_xyz[ddidx[n]], txyz2);
                        d2 = dot(zave, txyz2); /* check for normal terminal */
                        if (d * d2 < 0.0 && veclen(txyz2) <= helix_break) {
                            end_list[*num_ends][3] = ddidx[n];
                            bp_order[i][3] = ddidx[n];
                        }
                    }
                }
                // frprintf(tfp, "%4ld: %4ld %4ld %4ld %8.2f %8.2f%c", i, bp_order[i][1],
                //         bp_order[i][2], bp_order[i][3], ddmin[1], ddmin[n],
                //         (ddmin[n] > helix_break) ? '*' : ' ');
                if (!bp_order[i][2]) {
                    // frprintf(tfp, "  isolated base-pairs\n");
                } else {
                    // (!bp_order[i][3]) ? fprintf(tfp, " (%4ld)", ddidx[n]) : fprintf(tfp, "       ");
                    ddxyz(bp_xyz[ddidx[n]], bp_xyz[i], txyz2);
                    d2 = dot(zave, txyz2);
                    if (!overlap && !quadruple)
                        quadruple = bp_coplanar(i, d, d2, txyz, txyz2, n, ddidx, bp_xyz);
                    // fprintf(tfp, " ==> %8.2f %8.2f %c", d, d2, (d * d2 > 0.0) ? '*' : ' ');
                    d = magang(txyz, txyz2); /* magang normalize txyz & txyz2 */
                    // fprintf(tfp, " (%8.2f%c)\n", d, (d <= 90.0) ? '*' : ' ');
                }
            }
        }

        if (!*num_ends) {
            /* circular (as in 4nnu), or num_bp == 1 || 2 */
            end_list[++*num_ends][1] = 1;
            if (is_circular(num_bp, bp_order)) {
                /* normally m=2, n=3 */
                m = lval_min(bp_order[1][2], bp_order[1][3]); /* [1]: -1 2 44 */
                end_list[*num_ends][2] = m;                   /* [2]: -1 3 1 */
                n = (bp_order[m][2] == 1) ? bp_order[m][3] : bp_order[m][2];
                end_list[*num_ends][3] = n;
            } else if (num_bp == 2) {
                if (d <= helix_break) {
                    end_list[*num_ends][2] = 2; /* 1 2 0 && 2 1 0 */
                    end_list[++*num_ends][1] = 2;
                    end_list[*num_ends][2] = 1;
                } else {
                    end_list[++*num_ends][1] = 2; /* 1 0 0 && 2 0 0 */
                }
            }
        }
        // frprintf(tfp, "\nEnd base-pair list\n");
        // for (i = 1; i <= *num_ends; i++)
        //     frprintf(tfp, "%4ld: %4ld %4ld %4ld\n", i, end_list[i][1], end_list[i][2], end_list[i][3]);
    }

    /* locate all possible helical regions, including isolated base-pairs */
    static void locate_helix(long num_bp, long **helix_idx, long num_ends, long *num_helix,
                             long **end_list, long **bp_order, long *bp_idx, long *helix_marker)
    {
        long i, ip = 0, j, k, k0, k2, k3, m;
        long *matched_idx;

        helix_idx[*num_helix][1] = 1;

        matched_idx = lvector(1, num_bp); /* indicator for used bps */

        for (i = 1; i <= num_ends && ip < num_bp; i++) {
            k = 0;
            k0 = 0;
            for (j = 1; j <= 3; j++)
                if (end_list[i][j]) {
                    k += matched_idx[end_list[i][j]];
                    k0++;
                }
            if (k == k0)
                continue; /* end point of a processed helix */
            for (j = 1; j <= 3 && ip < num_bp; j++) {
                k = end_list[i][j];
                if (k && !matched_idx[k]) {
                    bp_idx[++ip] = k;
                    matched_idx[k] = 1;
                }
            }
            for (j = 1; j <= num_bp; j++) {
                k = bp_idx[ip];
                k2 = bp_order[k][2];
                k3 = bp_order[k][3];
                if (!bp_order[k][1]) {
                    /* end-point */
                    if (k2 && !matched_idx[k2] && !k3) { 
                        /* e.g. 0 10 0 */
                        bp_idx[++ip] = k2;
                        matched_idx[k2] = 1;
                    }
                    break; /* normal case */
                }
                m = matched_idx[k2] + matched_idx[k3];
                if (m == 2 || m == 0)
                    break; /* chain terminates */
                if (k2 == bp_idx[ip - 1]) {
                    bp_idx[++ip] = k3;
                    matched_idx[k3] = 1;
                } else if (k3 == bp_idx[ip - 1]) {
                    bp_idx[++ip] = k2;
                    matched_idx[k2] = 1;
                } else {
                    break; /* no direct connection */
                }
            }
            helix_idx[*num_helix][2] = ip;
            helix_marker[ip] = 1; /* helix_marker & helix_idx are parallel */
            if (ip < num_bp)
                helix_idx[++*num_helix][1] = ip + 1;
        }

        if (ip < num_bp) { 
            /* all un-classified bps */
            // // frprintf(stderr, "[%ld %ld]: complicated structure, left over"
            //                 " base-pairs put into the last region [%ld]\n",
            //         ip, num_bp, *num_helix);
            helix_idx[*num_helix][2] = num_bp;
            helix_marker[num_bp] = 1;
            for (j = 1; j <= num_bp; j++)
                if (!matched_idx[j])
                    bp_idx[++ip] = j;
        }
        free_lvector(matched_idx, 1, num_bp);
    }

    /* get the two base indices [*n1 and *n2] of base-pair m */
    static void get_ij(long m, long *swapped, long **base_pairs, long *n1, long *n2) {
        if (swapped[m]) {
            *n1 = base_pairs[m][2];
            *n2 = base_pairs[m][1];
        } else {
            *n1 = base_pairs[m][1];
            *n2 = base_pairs[m][2];
        }
    }

    /* make strand I of the first step in 5'--->3' direction */
    static void first_step(long i, long **helix_idx, long *bp_idx, long *swapped,
                           long **base_pairs, double **o3_p) {
        long i1, i2, j, j1_osx, j2, k, m, n;

        if (helix_idx[i][3] == 1)
            return;
        j = helix_idx[i][1];
        m = bp_idx[j];
        n = bp_idx[j + 1];
        get_ij(m, swapped, base_pairs, &i1, &j1_osx);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        k = is_linked(i1, i2, o3_p);
        if (k == -1)
            swapped[m] = !swapped[m];
        else if (!k) { 
            /* no linkage: reverse the strand direction */
            lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx);
            j = helix_idx[i][1];
            m = bp_idx[j];
            n = bp_idx[j + 1];
            get_ij(m, swapped, base_pairs, &i1, &j1_osx);
            get_ij(n, swapped, base_pairs, &i2, &j2);
            k = is_linked(i1, i2, o3_p);
            if (k == -1)
                swapped[m] = !swapped[m];
            else if (!k)                                            /* still no connection */
                lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx); /* reverse back */
        }
    }

    /* make sure strand I is in the 5'--->3' direction */
    static long chain1dir(long m, long n, long *swapped, long **base_pairs, double **o3_p) {
        long i1, i2, j1_osx, j2, k, irev0 = 0, irev1 = 1;

        get_ij(m, swapped, base_pairs, &i1, &j1_osx);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        k = is_linked(i1, i2, o3_p);
        return (k == -1) ? irev1 : irev0;
    }

    /* get base normal position index in bp_xyz */
    static void get_bidx(long m, long *swapped, long *idx1, long *idx2) {
        if (swapped[m]) {
            *idx1 = 18;
            *idx2 = 9;
        } else {
            *idx1 = 9;
            *idx2 = 18;
        }
    }

    /* get the dot product between base-pair m & n normals in WC geometry */
    static double wcbp_zdir(long m, long n, long idxm1, long idxm2, long idxn1, long idxn2,
                            double **bp_xyz) {
        double dm[4], dn[4];

        ddxyz(bp_xyz[m] + idxm2, bp_xyz[m] + idxm1, dm); /* opposite direction */
        ddxyz(bp_xyz[n] + idxn2, bp_xyz[n] + idxn1, dn);
        vec_norm(dm);
        vec_norm(dn);

        return dot(dm, dn);
    }

    /* get the dot product between base-pair m & n x-axes in WC geometry */
    static double wcbp_xang(long m, long n, double **bp_xyz) {
        double dm[4], dn[4];

        sumxyz(bp_xyz[m] + 3, bp_xyz[m] + 12, dm); /* same direction */
        sumxyz(bp_xyz[n] + 3, bp_xyz[n] + 12, dn);

        return magang(dm, dn);
    }

    /* get the relative base-pair orientations of WC-geometry step */
    static long wc_bporien(long m, long n, long *swapped, long **base_pairs, double **bp_xyz,
                           double **o3_p) {
        long i1, i2, j1_osx, j2, idxm1, idxm2, idxn1, idxn2, irev = 0;

        if (base_pairs[m][3] > 0 && base_pairs[n][3] > 0) { 
            /* WC geometry */
            get_ij(m, swapped, base_pairs, &i1, &j1_osx);
            get_ij(n, swapped, base_pairs, &i2, &j2);
            if (wcbp_xang(m, n, bp_xyz) > END_STACK_XANG || /* bdl070: 105 degrees */
                is_linked(i1, i2, o3_p) || is_linked(j1_osx, j2, o3_p))
                return irev; /* apposite bp orientation */
            get_bidx(m, swapped, &idxm1, &idxm2);
            get_bidx(n, swapped, &idxn1, &idxn2);
            if (wcbp_zdir(m, n, idxm1, idxm2, idxn1, idxn2, bp_xyz) < 0.0 &&
                wcbp_zdir(m, n, idxm1, idxm2, idxn2, idxn1, bp_xyz) > 0.0)
                irev = 1;
        }
        return irev;
    }

    /* use O3* distance criterion to decide if the top bp should be flipped:
 *    I       II
 *    i2 ---- j2
 *    |       |
 *    i1 ---- j1 */
    static long check_o3dist(long m, long n, long *swapped, long **base_pairs, double **o3_p) {
        double di1_i2, di1_j2, dj1_i2, dj1_j2;
        long i1, i2, j1_osx, j2, irev = 0;

        get_ij(m, swapped, base_pairs, &i1, &j1_osx);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        di1_i2 = distance_ab(o3_p, i1, i2, 4, 4);
        di1_j2 = distance_ab(o3_p, i1, j2, 4, 4);
        dj1_i2 = distance_ab(o3_p, j1_osx, i2, 4, 4);
        dj1_j2 = distance_ab(o3_p, j1_osx, j2, 4, 4);
        if ((di1_i2 > 0.0 && di1_j2 > 0.0 && di1_i2 > di1_j2) &&
            (dj1_i2 > 0.0 && dj1_j2 > 0.0 && dj1_j2 > dj1_i2))
            irev = 1;
        return irev;
    }

    /* use continuous single chain (if any) to control if to flip the top bp */
    static long check_schain(long m, long n, long *swapped, long **base_pairs, double **o3_p) {
        long i1, i2, j1_osx, j2, irev = 0;

        get_ij(m, swapped, base_pairs, &i1, &j1_osx);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        if (!is_linked(i1, i2, o3_p) && !is_linked(j1_osx, j2, o3_p) &&
            (is_linked(i1, j2, o3_p) || is_linked(j1_osx, i2, o3_p)))
            irev = 1;
        return irev;
    }

    /* check for situations where wc_bporien, check_o3dist & check_schain do not apply */
    static long check_others(long m, long n, long *swapped, long **base_pairs, double **o3_p,
                             double **bp_xyz) {
        double d[5], a1[4], a2[4], r1[4], r2[4];
        long i, j, i1, i2, j1_osx, j2, idxm1, idxm2, idxn1, idxn2;
        long irev0 = 0, irev1 = 1;

        /* check if any two residues are directly connected */
        get_ij(m, swapped, base_pairs, &i1, &j1_osx);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        if (is_linked(i1, i2, o3_p) || is_linked(j1_osx, j2, o3_p) ||
            is_linked(i1, j2, o3_p) || is_linked(j1_osx, i2, o3_p))
            return irev0;

        get_bidx(m, swapped, &idxm1, &idxm2);
        get_bidx(n, swapped, &idxn1, &idxn2);

        /* original orientation */
        for (i = 1; i <= 3; i++) {
            j = (3 - i) * 3;
            a1[i] = dot(&bp_xyz[m][idxm1 - j], &bp_xyz[n][idxn1 - j]);
            a2[i] = dot(&bp_xyz[m][idxm2 - j], &bp_xyz[n][idxn2 - j]);
        }
        i1 = (a1[1] > 0.0 && a1[2] > 0.0 && a1[3] > 0.0);
        i2 = (a2[1] > 0.0 && a2[2] > 0.0 && a2[3] > 0.0);
        if (i1 && i2) /* as is */
            return irev0;

        /* swap the two bases in the top base-pair */
        for (i = 1; i <= 3; i++) {
            j = (3 - i) * 3;
            r1[i] = dot(&bp_xyz[m][idxm1 - j], &bp_xyz[n][idxn2 - j]);
            r2[i] = dot(&bp_xyz[m][idxm2 - j], &bp_xyz[n][idxn1 - j]);
        }
        j1_osx = (r1[1] > 0.0 && r1[2] > 0.0 && r1[3] > 0.0);
        j2 = (r2[1] > 0.0 && r2[2] > 0.0 && r2[3] > 0.0);
        if (!i1 && !i2) {
            if (j1_osx || j2) /* reverse */
                return irev1;
            else if (!j1_osx && !j2) /* keep as is */
                return irev0;
        }
        d[1] = dot2ang(a1[1]) + dot2ang(a1[2]) + dot2ang(a1[3]); /* i1 */
        d[2] = dot2ang(a2[1]) + dot2ang(a2[2]) + dot2ang(a2[3]); /* i2 */
        d[3] = dot2ang(r1[1]) + dot2ang(r1[2]) + dot2ang(r1[3]); /* j1 */
        d[4] = dot2ang(r2[1]) + dot2ang(r2[2]) + dot2ang(r2[3]); /* j2 */

        if (i1 && j1_osx)
            return (d[1] > d[3]) ? irev1 : irev0;
        if (i1 && j2)
            return (d[1] > d[4]) ? irev1 : irev0;
        if (i2 && j1_osx)
            return (d[2] > d[3]) ? irev1 : irev0;
        if (i2 && j2)
            return (d[2] > d[4]) ? irev1 : irev0;

        return irev0;
    }

    /* check directions of two chains, mark and reverse as necessary a */
    static void check_direction(long i, long **helix_idx, long *bp_idx, long *swapped,
                                long **base_pairs, double **o3_p, long *direction) {
        long i1, i2, j, j1_osx, j2, k, m, n;

        /* check if strand I is in 5'-->3' direction */
        init_lvector(direction, 1, 6, 0);
        for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
            m = bp_idx[j];
            n = bp_idx[j + 1];
            get_ij(m, swapped, base_pairs, &i1, &j1_osx);
            get_ij(n, swapped, base_pairs, &i2, &j2);
            k = is_linked(i1, i2, o3_p);
            (k == 1) ? ++direction[1] : (k == -1) ? ++direction[2]
                                                  : ++direction[3];
            k = is_linked(j1_osx, j2, o3_p);
            (k == 1) ? ++direction[4] : (k == -1) ? ++direction[5]
                                                  : ++direction[6];
        }
        if ((direction[1] && direction[2]) || (direction[4] && direction[5])) {
            helix_idx[i][7] = 1;
            return;
        }
        if (direction[1] + direction[2] + direction[4] + direction[5] == 0)
            return; /* all empty: no direct connection */

        m = bp_idx[helix_idx[i][1]]; /* start bp */
        n = bp_idx[helix_idx[i][2]]; /* end bp */
        get_ij(m, swapped, base_pairs, &i1, &j1_osx);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        if (direction[3] || direction[6])
            helix_idx[i][5] = 1; /* broken O3'[i] to P[i+1] linkage */
        if (direction[1] && !direction[2]) {
            /* 5'--->3' direction */
            if (!direction[4] && direction[5]) {
                /* anti-parallel */
                if (i1 > j2) {
                    /* swap and reverse the two strands */
                    for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++)
                        swapped[bp_idx[j]] = !swapped[bp_idx[j]];
                    lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx);
                }
            }
            else if (direction[4] && !direction[5]) { 
                /* parallel */
                helix_idx[i][6] = 1;
                if (i1 > j1_osx) { /* swapped the two strands */
                    for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++)
                        swapped[bp_idx[j]] = !swapped[bp_idx[j]];
                }
            }
        }
    }

    /* direction of strand II: either anti-parallel or parallel to strand I */
    static void check_strand2(long i, long **helix_idx, long *bp_idx, double **bp_xyz,
                              long *swapped, long **base_pairs, double **o3_p,
                              long *direction, FILE *tfp) {
        long i1, i2, j, j1_osx, j2, k, m, n, anti_p, parallel;

        if (!helix_idx[i][7]) {
            if (direction[1] + direction[2] + direction[4] + direction[5] == 0)
                return;
            init_lvector(helix_idx[i], 5, 7, 0);
            for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
                m = bp_idx[j];
                n = bp_idx[j + 1];
                if (wc_bporien(m, n, swapped, base_pairs, bp_xyz, o3_p))
                    continue; /* in good WC-geometry orientation */
                get_ij(m, swapped, base_pairs, &i1, &j1_osx);
                get_ij(n, swapped, base_pairs, &i2, &j2);
                if (!is_linked(i1, i2, o3_p) && !is_linked(j1_osx, j2, o3_p) &&
                    ((is_linked(i1, j2, o3_p) == 1) ||
                     (is_linked(i1, j2, o3_p) && is_linked(j1_osx, i2, o3_p)))) {
                    swapped[n] = !swapped[n]; /* pr0004, ptr012, trna12 */
                    // frprintf(tfp, "                  000    [%ld-%ld]\n", m, n);
                    // fprintf(stderr, "000:    [%ld-%ld]\n", m, n);
                }
            }
        } else {
            helix_idx[i][7] = 0;
            anti_p = (direction[1] > direction[2]) && (direction[4] < direction[5]);
            parallel = (direction[1] > direction[2]) && (direction[4] > direction[5]);
            for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
                m = bp_idx[j];
                n = bp_idx[j + 1];
                get_ij(m, swapped, base_pairs, &i1, &j1_osx);
                get_ij(n, swapped, base_pairs, &i2, &j2);
                k = is_linked(j1_osx, j2, o3_p);
                if (!is_linked(i1, i2, o3_p) && /* strand I not linked */
                    ((anti_p && k == 1) || (parallel && k == -1))) {
                    swapped[n] = !swapped[n];
                    // frprintf(tfp, "                  2nd %2ld [%ld-%ld]\n", k, m, n);
                }
                /* another round of post-processing */
                get_ij(n, swapped, base_pairs, &i2, &j2);
                if (!is_linked(i1, i2, o3_p) && !is_linked(j1_osx, j2, o3_p)) {
                    if ((anti_p && is_linked(j1_osx, i2, o3_p) == 1) ||
                        (parallel && is_linked(i1, j2, o3_p) == -1)) {
                        // frprintf(tfp, "                  3rdL   [%ld-%ld]\n", m, n);
                        swapped[m] = !swapped[m]; /* lower bp */
                    }
                    else if ((anti_p && is_linked(i1, j2, o3_p) == 1) ||
                             (parallel && is_linked(j1_osx, i2, o3_p) == -1)) {
                        swapped[n] = !swapped[n]; /* upper bp */
                        // frprintf(tfp, "                  3rdU    [%ld-%ld]\n", m, n);
                    }
                }
            }
        }
        check_direction(i, helix_idx, bp_idx, swapped, base_pairs, o3_p, direction);
    }

    /* check WC-geometry step with negative rise */
    static void check_rise(long i, long **helix_idx, long *bp_idx, long *swapped,
                           long **base_pairs, double **bp_xyz, double **o3_p) {
        double d1, d2, rise, dorg[4], mn[4];
        long i1, j1_osx, i2, j2, idxm1, idxm2, idxn1, idxn2, j, k, m, n, num = 0;

        for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
            m = bp_idx[j];
            n = bp_idx[j + 1];
            if (base_pairs[m][3] > 0 && base_pairs[n][3] > 0) { 
                /* WC geometry */
                get_bidx(m, swapped, &idxm1, &idxm2);
                get_bidx(n, swapped, &idxn1, &idxn2);
                if (wcbp_zdir(m, n, idxm1, idxm2, idxn1, idxn2, bp_xyz) < 0.0 &&
                    wcbp_xang(m, n, bp_xyz) > END_STACK_XANG)
                    // fprintf(stderr, "^^vv opposite bp direction: %ld(%ld)"
                    //                 " %ld(%ld)-%ld(%ld)\n",
                    //         i, helix_idx[i][3], m, j, n, j + 1);
                ddxyz(bp_xyz[m], bp_xyz[n], dorg);
                for (k = 1; k <= 3; k++){
                    mn[k] = bp_xyz[m][k + idxm1] - bp_xyz[m][k + idxm2] + /* bp m */
                            bp_xyz[n][k + idxn1] - bp_xyz[n][k + idxn2];  /* bp n */
                }
                vec_norm(mn);
                rise = dot(dorg, mn);
                if (rise < 0.0) {
                    get_ij(m, swapped, base_pairs, &i1, &j1_osx);
                    get_ij(n, swapped, base_pairs, &i2, &j2);
                    d1 = distance_ab(o3_p, i2, i1, 4, 8);
                    d2 = distance_ab(o3_p, j1_osx, j2, 4, 8);
                    if (dval_in_range(d1, 0.0, O3P_UPPER) && dval_in_range(d2, 0.0, O3P_UPPER)) {
                        num++;
                        // fprintf(stderr, "===> %ld(%ld) %ld-%ld [%ld-%ld (%ld)]:"
                        //                 " %8.2f%8.2f%8.2f\n",
                        //         i, helix_idx[i][3], m, n,
                        //         j, j + 1, helix_idx[i][2], rise, d1, d2);
                    }
                }
            }
        }

    }

    /* make sure the leading strand is in the 5' to 3' direction
 * helix_idx[][7]: start#, ending#, total#, Z-DNA, break, parallel, weird
 *           col#     1       2       3       4      5       6        7
 */
    static void five2three(long num_bp, long *num_helix, long **helix_idx, long *bp_idx,
                           double **bp_xyz, long **base_pairs, double **o3_p, FILE *tfp) {
        double do3_p;
        long i, j, k, m, n, rev_wc, rev_o3d, rev_csc, rev_oth, rev_s1;
        long direction[7], *swapped;

        /* check if O3'[i] is wrongly connected to P[i] */
        for (i = 1; i <= *num_helix; i++)
            for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++)
            {
                k = base_pairs[bp_idx[j]][1];
                do3_p = distance_ab(o3_p, k, k, 4, 8);
                if (dval_in_range(do3_p, 0.0, O3P_UPPER))
                {
                    // fprintf(stderr, "Wrong: O3'[i] connected to P[i] (%ld/%ld/%g)\n", j, k, do3_p);
                    return; /* ignore 5'-->3' re-arrangement */
                }
            }

        // frprintf(tfp, "\nHelix region information\n");
        swapped = lvector(1, num_bp);
        for (i = 1; i <= *num_helix; i++) {
            helix_idx[i][3] = helix_idx[i][2] - helix_idx[i][1] + 1;
            print_sep(tfp, '-', 84);
            // frprintf(tfp, "Helix #%4.4ld\n", i);

            first_step(i, helix_idx, bp_idx, swapped, base_pairs, o3_p);

            /* The first base-pair is used as the global reference */
            for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
                m = bp_idx[j];
                n = bp_idx[j + 1];
                rev_wc = wc_bporien(m, n, swapped, base_pairs, bp_xyz, o3_p);
                rev_o3d = check_o3dist(m, n, swapped, base_pairs, o3_p);
                rev_csc = check_schain(m, n, swapped, base_pairs, o3_p);
                rev_oth = check_others(m, n, swapped, base_pairs, o3_p, bp_xyz);
                // frprintf(tfp, "          %4ld: %2ld %2ld %2ld %2ld", j, rev_wc, rev_o3d,
                //         rev_oth, rev_csc);
                if (rev_wc)
                    swapped[n] = !swapped[n];
                else if (rev_o3d || rev_csc || rev_oth) /* rev_o3d: udf027 */
                    swapped[n] = !swapped[n];
                /* 5'--->3' direction of strand I as final proof */
                rev_s1 = chain1dir(m, n, swapped, base_pairs, o3_p);
                if (rev_s1)
                    swapped[n] = !swapped[n];
                // frprintf(tfp, " %2ld [%ld-%ld]\n", rev_s1, m, n);
            }

            /* one more round checking for WC geometry steps: 09-23-2002
         * works for Luger pol_end1.pdb & pol_end2.pdb, and bdl070.pdb */
            // frprintf(tfp, "\n              ===> 2nd around checking or WC geometry steps\n");
            for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
                m = bp_idx[j];
                n = bp_idx[j + 1];
                rev_wc = wc_bporien(m, n, swapped, base_pairs, bp_xyz, o3_p);
                if (rev_wc) {
                    swapped[m] = !swapped[m]; /* m here */
                    // frprintf(tfp, "          %4ld: [%ld-%ld]\n", j, m, n);
                }
            }
            // frprintf(tfp, "\n");

            check_direction(i, helix_idx, bp_idx, swapped, base_pairs, o3_p, direction);
            check_strand2(i, helix_idx, bp_idx, bp_xyz, swapped, base_pairs, o3_p, direction, tfp);
            check_rise(i, helix_idx, bp_idx, swapped, base_pairs, bp_xyz, o3_p);

            for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
                m = bp_idx[j];
                if (swapped[m]) {
                    lval_swap(&base_pairs[m][1], &base_pairs[m][2]);
                    for (k = 1; k <= 9; k++) { 
                        /* swap base I and II xyz-axes */
                        dval_swap(&bp_xyz[m][k + 3], &bp_xyz[m][k + 12]);
                        lval_swap(&base_pairs[m][k + 11], &base_pairs[m][k + 20]);
                    }
                }
            }

            // frprintf(tfp, "\n%4ld [%c%c%c]:", helix_idx[i][3], helix_idx[i][5] ? 'b' : '-',
                    // helix_idx[i][6] ? 'p' : '-', helix_idx[i][7] ? '?' : '-');
            // for (j = 1; j <= 6; j++)
            //     frprintf(tfp, "%6ld", direction[j]);
            // frprintf(tfp, "\n           ");
            k = 0;
            for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
                m = bp_idx[j];
                // frprintf(tfp, "%6ld", swapped[m] ? -m : m);
                // if (++k % 10 == 0 && k != helix_idx[i][3])
                //     frprintf(tfp, "\n           ");
            }
            // frprintf(tfp, "\n");
        }

        free_lvector(swapped, 1, num_bp);
    }

    /* check if a helical region is in left-handed Z-form */
    static void check_zdna(long *num_helix, long **helix_idx, long *bp_idx, double **bp_xyz,
                           long **base_pairs, FILE *tfp) {
        double txyz[4];
        long i, j, m, n, nweird = 0, nrev, mixed_rl = 0;

        // frprintf(tfp, "\nZ-DNA helical region if any\n");
        for (i = 1; i <= *num_helix; i++) {
            if (helix_idx[i][5] || helix_idx[i][6] || helix_idx[i][7] || helix_idx[i][3] <= 1) {
                nweird++;
                continue; /* break/parallel/weird/only one pair */
            }
            nrev = 0;
            for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
                m = bp_idx[j];
                if (helix_idx[i][3] == 1)
                    continue;
                if (j < helix_idx[i][2]) {
                    n = bp_idx[j + 1];
                    ddxyz(bp_xyz[m], bp_xyz[n], txyz);
                }
                /* with base 1 normal & WC-geometry */
                if (dot(txyz, &bp_xyz[m][9]) < 0.0 && base_pairs[m][3] > 0)
                    nrev++;
                else
                    break;
            }
            if (nrev == helix_idx[i][3]) {
                helix_idx[i][4] = 1;
                mixed_rl++;
                // frprintf(tfp, "Helix #%4.4ld (%4ld) is a Z-DNA\n", i, helix_idx[i][3]);
            }
        }

        // if (!nweird && mixed_rl && mixed_rl != *num_helix)
        //     frprintf(stderr, "This structure has right-/left-handed helical regions\n");
    }

    /* Kth base-pair classification: WC, WC-geometry, Anti-parallel */
    static void set_wc3(long *pair_k, char *wc) {
        char zdir;
        double z1[4], z2[4];
        long i, j;

        for (i = 18; i <= 20; i++) {
            j = i - 17;
            z1[j] = pair_k[i] / MFACTOR;
            z2[j] = pair_k[i + 9] / MFACTOR;
        }
        zdir = (dot(z1, z2) < 0.0) ? '-' : '+';
        get_bp_3char_symbols(pair_k[3], zdir, wc);
    }

    /* re-order base-pairs into separate helical regions */
    static void re_ordering(long num_bp, long **base_pairs, long *bp_idx, long *helix_marker,
                            long **helix_idx, miscPars *misc_pars, long *num_helix,
                            double **o3_p, char *bseq, long **seidx, char **ResName,
                            char *ChainID, long *ResSeq, char **Miscs) {
        char b1[BUF512], b2[BUF512], wc[4];
        double **bp_xyz;
        long i, i_order, j, j_order, num_ends = 0;
        long **bp_order, **end_list;
        FILE *tfp;

        tfp = open_file(BPORDER_FILE, "w");
        print_bp_crit(misc_pars, tfp);
        // frprintf(tfp, "Base-pair information BEFORE re-ordering\n");
        for (i = 1; i <= num_bp; i++) {
            set_wc3(base_pairs[i], wc);
            i_order = base_pairs[i][1];
            j_order = base_pairs[i][2];
            j = seidx[i_order][1];
            base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[i_order], 1, b1);
            j = seidx[j_order][1];
            base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[j_order], 2, b2);
            // frprintf(tfp, "%5ld: %5ld %5ld %s-%s-%s", i, i_order, j_order, b1, wc, b2);
            // for (j = 4; j <= 8; j++) /* d, dv, angle, dNN */
            //     // frprintf(tfp, " %6.2f", base_pairs[i][j] / MFACTOR);
            // for (j = 9; j <= 11; j++) /* bp origin */
            //     frprintf(tfp, " %8.2f", base_pairs[i][j] / MFACTOR);
            // frprintf(tfp, "\n");
        }

        /* bp origin, base I/II xyz-axes copied from columns 9-29 of base_pairs
     *       bp_orgin,   x1,    y1,   z1,    x2,    y2,    z2
     * col#     1-3     4-6    7-9   10-12  13-15  16-18  19-21 */
        bp_xyz = dmatrix(1, num_bp, 1, 21);
        for (i = 1; i <= num_bp; i++)
            for (j = 1; j <= 21; j++)
                bp_xyz[i][j] = base_pairs[i][j + 8] / MFACTOR;

        bp_order = lmatrix(1, num_bp, 1, 3);
        end_list = lmatrix(1, num_bp, 1, 3);

        bp_context(num_bp, misc_pars, bp_xyz, bp_order, end_list, &num_ends, tfp);

        locate_helix(num_bp, helix_idx, num_ends, num_helix, end_list, bp_order, bp_idx, helix_marker);

        five2three(num_bp, num_helix, helix_idx, bp_idx, bp_xyz, base_pairs, o3_p, tfp);

        check_zdna(num_helix, helix_idx, bp_idx, bp_xyz, base_pairs, tfp);

        close_file(tfp);
        free_dmatrix(bp_xyz, 1, num_bp, 1, 21);
        free_lmatrix(bp_order, 1, num_bp, 1, 3);
        free_lmatrix(end_list, 1, num_bp, 1, 3);
    }

    /* print out a summary of helix information: Z-DNA, break, parallel, weird? */
    static void helix_info(long **helix_idx, long idx, FILE *fp) {
        // frprintf(fp, "%s%s%s%s\n", (helix_idx[idx][4]) ? "  ***Z-DNA***" : "",
        //         (helix_idx[idx][5]) ? "  ***broken O3'[i] to P[i+1] linkage***"
        //                             : "",
        //         (helix_idx[idx][6]) ? "  ***parallel***" : "",
        //         (helix_idx[idx][7]) ? "  ***intra-chain direction reverse***" : "");
    }

    /* write out a PDB file containing all best pairs oriented with mean base-pair normal */
    void
    PairFinder::_write_bestpairs(long num_bp, long **base_pairs, long *bp_idx, char *bseq,
                                long **seidx, char **AtomName, char **ResName, char *ChainID,
                                long *ResSeq, char **Miscs, double **xyz, double **orien,
                                double **org, long **htm_water, miscPars *misc_pars) {

        char b1[BUF512], b2[BUF512], idmsg[BUF512], wc[4];
        double morg[4], morien[10];
        long i, ia, ib, j, k;
        long inum_base = 2, ivec[3]; /* for middle frame orientation */
        FILE *mfp, *rframe;

        mfp = open_file(BESTP_FILE, "w");
        rframe = open_file(REF_FILE, "w");

        // frprintf(rframe, "%5ld base-pairs\n", num_bp);

        for (i = 1; i <= num_bp; i++) {
            k = bp_idx[i];
            ia = base_pairs[k][1];
            ib = base_pairs[k][2];
            j = seidx[ia][1];
            base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[ia], 1, b1);
            j = seidx[ib][1];
            base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[ib], 2, b2);
            set_wc3(base_pairs[k], wc);
            sprintf(idmsg, "%s-%s-%s", b1, wc, b2);
            // frprintf(mfp, "%6s    %4ld\n", "MODEL ", i);
            // frprintf(mfp, "REMARK    Section #%4.4ld %s\n", i, idmsg);
            // frprintf(mfp, "REMARK    %s\n", Gvars.X3DNA_VER);
            ivec[1] = ia;
            ivec[2] = ib;
            pair2mst(inum_base, ivec, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                     orien, org, seidx, morien, morg, htm_water, misc_pars, mfp);
            // frprintf(mfp, "ENDMDL\n");

            // frprintf(rframe, "... %5ld %c%c%c   # %s - %s\n", i, bseq[ia], wc[2], bseq[ib], nt_info[ia], nt_info[ib]);
            auto r = std::regex(
            "#\\s+(?:\\.+\\d+\\>)*(\\w+):\\.*(-*\\d+)\\S:\\[\\.*(\\S+)\\](\\w+)\\s+\\-\\s+(?:\\.+\\d+\\>)*(\\w+):\\.*(-*\\d+)\\S:\\[\\.*(\\S+)\\](\\w+)");
            auto s = string_format("... %5ld %c%c%c   # %s - %s\n", i, bseq[ia], wc[2], bseq[ib], nt_info[ia], nt_info[ib]);
            auto m = std::smatch();
            std::regex_search(s, m, r);
            auto bp_info = new X3dna::X3BPInfo(m);
            PairFinder::_write_fpmst(morg, morien, rframe, bp_info);
        }
    }

    /* write out helical regions in a multiple structure PDB file */
    static void write_helix(long num_helix, long **helix_idx, long *bp_idx, long **seidx,
                            char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                            char **Miscs, double **xyz, long **base_pairs, long **htm_water,
                            miscPars *misc_pars) {
        long i, inum, j, k, m, n, tnum_res, num_residue;
        long *ivec, *ivect;
        FILE *mfp;

        num_residue = htm_water[2][0];
        ivec = lvector(1, num_residue);
        ivect = lvector(1, num_residue);
        mfp = open_file(HLXREG_FILE, "w");
        for (i = 1; i <= num_helix; i++) {
            inum = 0;
            // frprintf(mfp, "%6s    %4ld\n", "MODEL ", i);
            // frprintf(mfp, "REMARK    Section #%4.4ld %ld base-pairs", i, helix_idx[i][3]);
            helix_info(helix_idx, i, mfp);
            // frprintf(mfp, "REMARK    %s\n", Gvars.X3DNA_VER);
            k = 0; /* residue index in the helical region */
            for (n = 1; n <= 2; n++) { 
                /* two strands */
                for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
                    if (n == 2 && !helix_idx[i][6]) /* anti-parallel */
                        m = base_pairs[bp_idx[helix_idx[i][2] - j + helix_idx[i][1]]][n];
                    else
                        m = base_pairs[bp_idx[j]][n];
                    k++;
                    ivec[k] = m;
                }
            }
            tnum_res = attached_residues(k, ivec, ivect, seidx, xyz, htm_water, misc_pars);
            for (j = 1; j <= tnum_res; j++) {
                m = ivect[j];
                pdb_record(seidx[m][1], seidx[m][2], &inum, 0, AtomName, ResName,
                           ChainID, ResSeq, xyz, Miscs, mfp);
            }
            // frprintf(mfp, "ENDMDL\n");
        }

        close_file(mfp);
        free_lvector(ivec, 1, num_residue);
        free_lvector(ivect, 1, num_residue);
    }

    /* for non-base-pairs found in structure "pdbfile" */
    static void no_basepairs(char *pdbfile, char *outfile, char *parfile) {
        FILE *fp;
        // // frprintf(stderr, "no base-pairs found for this structure\n");
        fp = open_file(outfile, "w");
        // frprintf(fp, "%s\n", pdbfile);
        // frprintf(fp, "%s.out\n", parfile);
        // frprintf(fp, "    2         # duplex\n");
        // frprintf(fp, "    0         # number of base-pairs\n");
        close_file(fp);
    }

    /* get the best pair listing */
    static long find_bestpair(long nout, long **base_pairs, long num_residue, char *bseq,
                              long **seidx, long *RY, char **AtomName, double **xyz,
                              long *idx, double **orien, double **org, double **NC1xyz,
                              long **ring_atom, miscPars *misc_pars) {
        long i, j, num1 = 0, num2 = 1, num_bp = 0;
        long pair_istat[PSTNUM], pair_jstat[PSTNUM];
        long *matched_idx;

        matched_idx = lvector(1, num_residue);
        while (num1 < num2) {
            num1 = num2;
            for (i = 1; i <= num_residue; i++) {
                if (RY[i] < 0 || matched_idx[i])
                    continue; /* non-base or paired base */
                best_pair(i, num_residue, RY, seidx, xyz, idx, NC1xyz, matched_idx,
                          orien, org, ring_atom, AtomName, bseq, misc_pars, pair_istat);
                if (pair_istat[1]) { 
                    /* with paired base */
                    best_pair(pair_istat[1], num_residue, RY, seidx, xyz, idx, NC1xyz,
                              matched_idx, orien, org, ring_atom, AtomName, bseq,
                              misc_pars, pair_jstat);
                    if (i == pair_jstat[1]) { 
                        /* best match between i && pair_istat[1] */
                        matched_idx[i] = 1;
                        matched_idx[pair_istat[1]] = 1;
                        base_pairs[++num_bp][1] = i;
                        for (j = 1; j <= nout; j++)
                            base_pairs[num_bp][j + 1] = pair_istat[j];
                    }
                }
            }
            num2 = 0;
            for (i = 1; i <= num_residue; i++)
                if (matched_idx[i])
                    num2++;
        }
        free_lvector(matched_idx, 1, num_residue);

        return num_bp;
    }

    /* RasMol script for coloring the structure by chains or helical regions */
    static void col_helices(long num_helix, long **helix_idx, long *bp_idx, long **base_pairs,
                            long **seidx, char *pdbfile, char *ChainID, long *ResSeq) {
        static char *col_code[9] = {"violet", "red", "green", "blue", "yellow",
                                    "cyan", "magenta", "orange", "purple"};
        long i, ia, ib, ic, j, k;
        FILE *fpc, *fph;


        fpc = open_file(COLCHN_FILE, "w");

        // frprintf(fpc, "zap\nload nmrpdb hel_regions.pdb\n");
        // frprintf(fpc, "# load %s\n", pdbfile);
        // frprintf(fpc, "# restrict not (protein or water)\n");
        // frprintf(fpc, "\n");


        fph = open_file(COLHLX_FILE, "w");
        // frprintf(fph, "zap\nload nmrpdb hel_regions.pdb\n");
        // frprintf(fph, "# load %s\n", pdbfile);
        // frprintf(fph, "# restrict not (protein or water)\n");
        // frprintf(fph, "\n");

        for (i = 1; i <= num_helix; i++) {
            ic = i % 9; /* color code */
            // frprintf(fph, "\n#------Helix #%ld, color: %s------\n", i, col_code[ic]);
            for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
                k = bp_idx[j];
                ia = seidx[base_pairs[k][1]][1];
                ib = seidx[base_pairs[k][2]][1];
                // frprintf(fpc, "select %ld:%c\n", ResSeq[ia], ChainID[ia]);
                // frprintf(fpc, "color %s\n", col_code[1]); /* strand I */
                // frprintf(fpc, "select %ld:%c\n", ResSeq[ib], ChainID[ib]);
                // frprintf(fpc, "color %s\n", col_code[2]); /* strand II */
                // frprintf(fph, "select %ld:%c, %ld:%c\n", ResSeq[ia], ChainID[ia],
                //         ResSeq[ib], ChainID[ib]);
                // frprintf(fph, "color %s\n", col_code[ic]);
            }
        }

        // frprintf(fpc, "\nselect all\n");
        close_file(fpc);
        // frprintf(fph, "\nselect all\n");
        close_file(fph);
    }

    static void set_nmarkers(long idx, long ib, long ie, long *helix_marker, long **helix_idx,
                             long *num_helix, long *num_1bp, long *nmarkers) {
        long i, k;

        *num_helix = 0;
        *num_1bp = 0;

        for (i = ib; i <= ie; i++) {
            k = i - ib + 1; /* 1-index */

            if (helix_marker[i]) { 
                /* change of helix regions */
                (*num_helix)++;
                if ((!idx && helix_idx[*num_helix][3] == 1) || (idx && helix_idx[idx][3] == 1)) {
                    (*num_1bp)++;
                    nmarkers[k] = 1;
                }
                else if (i != ie)
                    nmarkers[k] = 9;
            }
        }
    }

    /* print out base-pairing information in a format suitable for analyze/cehs */
    static void x3dna_input(long idx, long start_num, long end_num, long nbp, char *pdbfile,
                            char *outfile, char *parfile, long hetatm, long *bp_idx,
                            long *helix_marker, long **helix_idx, long **base_pairs,
                            long **seidx, char **ResName, char *ChainID, long *ResSeq,
                            char **Miscs, char *bseq, miscPars *misc_pars, long detailed) {
        char b1[BUF512], b2[BUF512], wc[4], *cmarkers;
        char outfile_new[BUF512], parfile_new[BUF512];
        double x[4];
        long i, i_order, j, j_order, k, m, *nmarkers;
        long num_bp, num_helix, num_1bp, num_nwc = 0;
        FILE *fp;

        if (idx) { 
            /* structure taken as helical regions */
            sprintf(outfile_new, "%s_%4.4ld", outfile, idx);
            sprintf(parfile_new, "%s_%4.4ld", parfile, idx);
        }
        else { 
            /* structure taken as a whole */
            sprintf(outfile_new, outfile);
            sprintf(parfile_new, parfile);
        }

        fp = open_file(outfile_new, "w");

        // frprintf(fp, "%s\n", pdbfile);
        // frprintf(fp, "%s.out\n", parfile_new);
        // frprintf(fp, "    2         # duplex\n");
        // frprintf(fp, "%5ld         # number of base-pairs\n", nbp);
        // frprintf(fp, "    1 %5ld    # explicit bp numbering/hetero atoms\n", hetatm);

        num_bp = end_num - start_num + 1;
        cmarkers = cvector(1, num_bp);
        nmarkers = lvector(1, num_bp);
        set_nmarkers(idx, start_num, end_num, helix_marker, helix_idx, &num_helix, &num_1bp, nmarkers);
        set_chain_nmarkers019_to_symbols(num_bp, nmarkers, cmarkers);

        for (i = start_num; i <= end_num; i++) {
            k = bp_idx[i];
            i_order = base_pairs[k][1];
            j_order = base_pairs[k][2];
            if (base_pairs[k][3] != 2)
                num_nwc++;
            set_wc3(base_pairs[k], wc);
            j = seidx[i_order][1];
            base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[i_order], 1, b1);
            j = seidx[j_order][1];
            base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[j_order], 2, b2);

            m = i - start_num + 1;
            // frprintf(fp, "%5ld %5ld %3ld #%5ld %c %s-%s-%s", i_order, j_order, nmarkers[m], m,
            //         cmarkers[m], b1, wc, b2);

            for (j = 4; j <= 8; j++)
                // frprintf(fp, " %6.2f", base_pairs[k][j] / MFACTOR);
            if (detailed) { 
                /* print out origin + x-axis */
                for (j = 9; j <= 11; j++)
                    // frprintf(fp, " %8.2f", base_pairs[k][j] / MFACTOR);
                for (j = 1; j <= 3; j++)
                    x[j] = (base_pairs[k][j + 11] + (base_pairs[k][j + 20])) / MFACTOR;
                vec_norm(x);
                // for (j = 1; j <= 3; j++)
                //     frprintf(fp, " %8.2f", x[j]);
            }
            // frprintf(fp, "\n");
        }

        // frprintf(fp, "##### ");
        print_bp_crit(misc_pars, fp);
        // frprintf(fp, "##### %ld non-Watson-Crick base-pair%s", num_nwc, (num_nwc == 1) ? "" : "s");
        (num_helix == 1) ? strcpy(b1, "x") : strcpy(b1, "ces");
        // frprintf(fp, ", and %ld heli%s", num_helix, b1);
        // frprintf(fp, " (%ld isolated bp%s)\n", num_1bp, (num_1bp == 1) ? "" : "s");

        if (!idx) { 
            /* structure taken as a whole */
            for (i = 1; i <= num_helix; i++) {
                // if (helix_idx[i][3] == 1)
                //     // frprintf(fp, "##### Helix #%ld (%ld): %ld", i, helix_idx[i][3], helix_idx[i][1]);
                // else
                //     // frprintf(fp, "##### Helix #%ld (%ld): %ld - %ld", i,
                //             helix_idx[i][3], helix_idx[i][1], helix_idx[i][2]);
                helix_info(helix_idx, i, fp);
            }
        }
        else { 
            /* for each helical region */
            // if (nbp == 1)
            //     // frprintf(fp, "##### Helix #1 (%ld): %ld", nbp, nbp);
            // else
            //     // frprintf(fp, "##### Helix #1 (%ld): 1 - %ld", nbp, nbp);
            helix_info(helix_idx, idx, fp);
        }

        close_file(fp);
        free_cvector(cmarkers, 1, DUMMY);
        free_lvector(nmarkers, 1, DUMMY);
    }

    /* print out base-pairing information in a format suitable for Curves */
    static void curves_input(long idx, long start_num, long end_num, long nbp, char *pdbfile,
                             char *outfile, char *parfile, long *bp_idx, long **base_pairs,
                             long zdna, long parallel) {
        char outfile_new[BUF512], parfile_new[BUF512];
        long i, j;
        FILE *fp;

        if (idx) { 
            /* structure taken as helical regions */
            sprintf(outfile_new, "%s_%4.4ld", outfile, idx);
            sprintf(parfile_new, "%s_%4.4ld", parfile, idx);
        }
        else { 
            /* structure taken as a whole */
            strcpy(outfile_new, outfile);
            strcpy(parfile_new, parfile);
        }
        fp = open_file(outfile_new, "w");
        // frprintf(fp, "&inp file=%s, comb=.t., fit=.t., grv=.t., %s\n"
        //             "     lis=%s, pdb=%s_grp, &end\n",
        //         pdbfile,
        //         zdna ? "dinu=.t.," : "", parfile_new, parfile_new);
        // frprintf(fp, "2 %ld %ld 0 0\n", nbp, (parallel) ? nbp : -nbp);
        // for (i = 1; i <= 2; i++) {
        //     for (j = start_num; j <= end_num; j++)
        //         frprintf(fp, " %ld", base_pairs[bp_idx[j]][i]);
        //     frprintf(fp, "\n");
        // }
        // frprintf(fp, "0.0 0.0 0.0 %.1f\n", zdna ? 180.0 : 0.0);
        close_file(fp);
    }

    /* print out base-pairing information in a format suitable for Curves+ */
    static void curves_plus_input(long idx, long start_num, long end_num, long nbp,
                                  char *pdbfile, char *outfile, char *parfile, long *bp_idx,
                                  long **base_pairs, long zdna, long parallel) {
        char outfile_new[BUF512], parfile_new[BUF512];
        char *p, str[BUF512];
        char *cmnfiles[] = {".cda", ".lis", "_X.pdb", "_b.pdb"};
        long num = sizeof cmnfiles / sizeof cmnfiles[0];
        long i, j;
        FILE *fp;

        UNUSED_PARAMETER(nbp); /* not used, but kept here for consistency */

        if (idx) { 
            /* structure taken as helical regions */
            sprintf(outfile_new, "%s_%4.4ld", outfile, idx);
            sprintf(parfile_new, "%s_%4.4ld", parfile, idx);
        }
        else { 
            /* structure taken as a whole */
            strcpy(outfile_new, outfile);
            strcpy(parfile_new, parfile);
        }

        for (i = 0; i < num; i++) {
            sprintf(str, "%s%s", parfile_new, cmnfiles[i]);
            remove_file(str);
        }

        p = getenv("CURVES_PLUS_STDLIB");
        if (p) {
            strcpy(str, p);
            check_slash(str); /* add slash at the end */
        }
        else
            strcpy(str, "./");

        fp = open_file(outfile_new, "w");
        // frprintf(fp, "&inp file=%s,\n", pdbfile);
        // frprintf(fp, "     lis=%s,\n", parfile_new);
        // frprintf(fp, "     fit=.t.,\n");
        // frprintf(fp, "     lib=%sstandard,\n", str);
        // frprintf(fp, "     isym=%d,\n", zdna ? 2 : 1);
        // frprintf(fp, "&end\n");

        // frprintf(fp, "2 %d %d 0 0\n", 1, (parallel) ? 1 : -1);
        // for (i = 1; i <= 2; i++) {
        //     for (j = start_num; j <= end_num; j++)
        //         frprintf(fp, " %ld", base_pairs[bp_idx[j]][i]);
        //     frprintf(fp, "\n");
        // }
        close_file(fp);
    }

    void
    PairFinder::_duplex(long num, long num_residue, char *bseq, long **seidx, long *RY,
                       char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                       char **Miscs, double **xyz, char *parfile,
                       miscPars *misc_pars) {
        double **orien, **org, **NC1xyz, **o3_p;
        long i, num_bp, num_helix = 1;
        long nout = PSTNUM - 1, nout_p1 = PSTNUM;
        long *idx, *bp_idx, *helix_marker, **htm_water;
        long **base_pairs, **helix_idx, **ring_atom;

        orien = dmatrix(1, num_residue, 1, 9);
        org = dmatrix(1, num_residue, 1, 3);
        NC1xyz = dmatrix(1, num_residue, 1, 7); /* RN9/YN1 & C1' atomic coordinates */
        o3_p = dmatrix(1, num_residue, 1, 8);   /* O3'/P atomic coordinates */

        idx = lvector(1, num);
        atom_idx(num, AtomName, NULL, idx);

        htm_water = lmatrix(1, 4, 0, num); /* HETATM and water index */
        init_htm_water(args.waters, num, num_residue, idx, htm_water);
        identify_htw(num_residue, seidx, RY, AtomName, ResName, ChainID, ResSeq,
                     Miscs, xyz, htm_water);

        base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
                  ResSeq, Miscs, xyz, orien, org, NC1xyz, o3_p);

        /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
        ring_atom = lmatrix(1, num_residue, 1, 19);
        ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);

        /* find best base-pairs */
        base_pairs = lmatrix(1, num_residue, 1, nout_p1);
        /* base_pairs[][PSTNUM]:
     *      i, j, bpid, d, dv, angle, dNN, dsum, bp-org,
     * col# 1  2    3   4  5     6     7     8    9-11
     *       x1,   y1,   z1,   x2,   y2,   z2
     *      12-14 15-17 18-20 21-23 24-26 27-29
     * i.e., add one more column for i from pair_stat [best_pair]
     */
        num_bp = find_bestpair(nout, base_pairs, num_residue, bseq, seidx, RY, AtomName,
                               xyz, idx, orien, org, NC1xyz, ring_atom, misc_pars);
        if (!num_bp) {
            no_basepairs(args.pdbfile, args.outfile, parfile);
            goto NO_BASE_PAIR; /* to clean up */
        }

        bp_idx = lvector(1, num_bp);
        helix_marker = lvector(1, num_bp);
        helix_idx = lmatrix(1, num_bp, 1, 7);
        re_ordering(num_bp, base_pairs, bp_idx, helix_marker, helix_idx, misc_pars,
                    &num_helix, o3_p, bseq, seidx, ResName, ChainID, ResSeq, Miscs);
        PairFinder::_write_bestpairs(num_bp, base_pairs, bp_idx, bseq, seidx, AtomName, ResName, ChainID,
                        ResSeq, Miscs, xyz, orien, org, htm_water, misc_pars);
        write_helix(num_helix, helix_idx, bp_idx, seidx, AtomName, ResName, ChainID, ResSeq,
                    Miscs, xyz, base_pairs, htm_water, misc_pars);

        if (args.curves)
        { /* generate Curves input file */
            if (args.divide && num_helix > 1)
                for (i = 1; i <= num_helix; i++)
                    curves_input(i, helix_idx[i][1], helix_idx[i][2], helix_idx[i][3],
                                 args.pdbfile, args.outfile, parfile, bp_idx, base_pairs,
                                 helix_idx[i][4], helix_idx[i][6]);
            else
                curves_input(0, 1, num_bp, num_bp, args.pdbfile, args.outfile, parfile,
                             bp_idx, base_pairs, helix_idx[1][4], helix_idx[1][6]);
        }
        else if (args.curves_plus)
        {
            if (args.divide && num_helix > 1)
                for (i = 1; i <= num_helix; i++)
                    curves_plus_input(i, helix_idx[i][1], helix_idx[i][2], helix_idx[i][3],
                                      args.pdbfile, args.outfile, parfile, bp_idx,
                                      base_pairs, helix_idx[i][4], helix_idx[i][6]);
            else
                curves_plus_input(0, 1, num_bp, num_bp, args.pdbfile, args.outfile, parfile,
                                  bp_idx, base_pairs, helix_idx[1][4], helix_idx[1][6]);
        }
        else { 
            /* 3DNA/CEHS input file */

            col_helices(num_helix, helix_idx, bp_idx, base_pairs, seidx, args.pdbfile,
                        ChainID, ResSeq);

            if (args.divide && num_helix > 1)
                for (i = 1; i <= num_helix; i++) {

                    x3dna_input(i, helix_idx[i][1], helix_idx[i][2], helix_idx[i][3],
                                args.pdbfile, args.outfile, parfile, args.hetatm, bp_idx,
                                helix_marker, helix_idx, base_pairs, seidx, ResName, ChainID,
                                ResSeq, Miscs, bseq, misc_pars, args.detailed);
                }

            else /* overall */

                x3dna_input(0, 1, num_bp, num_bp, args.pdbfile, args.outfile, parfile,
                            args.hetatm, bp_idx, helix_marker, helix_idx, base_pairs, seidx,
                            ResName, ChainID, ResSeq, Miscs, bseq, misc_pars, args.detailed);

        }
        free_lvector(bp_idx, 1, num_bp);
        free_lvector(helix_marker, 1, num_bp);
        free_lmatrix(helix_idx, 1, num_bp, 1, 7);

    NO_BASE_PAIR:
        free_lmatrix(base_pairs, 1, num_residue, 1, nout_p1);
    }

    static long read_mapping_table(char *cvt_table, char n1[][5], char n2[][5])
    {
        char t1[BUF512], t2[BUF512], *p0, *line;
        long num = 0;
        FILE *fp;

        fp = open_file(cvt_table, "r");

        while ((p0 = my_getline(fp)) != NULL)
        {
            line = trim(p0); /* keep the original value of p0 */
            if (line[0] == '#' || line[0] == '\0')
            {
                free(p0);
                continue; /* skip empty and commented lines */
            }

            upperstr(line);
            if (sscanf(line, "%s %s", t1, t2) != 2 || strlen(t1) != 4 || strlen(t2) != 4)
            {
                // frprintf(stderr, "invalid line: <%s>\n", p0);
                free(p0);
                continue;
            }

            cvtstr_c1toc2(t1, '_', ' ');
            cvtstr_c1toc2(t2, '_', ' ');
            num++;
            strcpy(n1[num], t1);
            strcpy(n2[num], t2);

            free(p0);
        }

        close_file(fp);

        return num; /* 1-indexed, actual number */
    }

    static void write_atom_coordinates(FILE *fp, long *serial, long idx, char **AtomName,
                                       char **ResName, char *ChainID, long *ResSeq,
                                       char **Miscs, double **xyz)
    {
        char str[BUF512];

        (*serial)++;
        deduce_misc(Miscs, AtomName, idx, str);
        // frprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
        //         (str[0] == 'A') ? "ATOM  " : "HETATM", *serial, AtomName[idx],
        //         str[1], ResName[idx], ChainID[idx], ResSeq[idx], str[2], xyz[idx][1],
        //         xyz[idx][2], xyz[idx][3], str + 3);
    }

    static char get_standard_one_letter_nt(char *rname)
    {
        if (is_equal_string(rname, "  A") || is_equal_string(rname, " DA") ||
            is_equal_string(rname, "ADE"))
            return 'A';
        if (is_equal_string(rname, "  C") || is_equal_string(rname, " DC") ||
            is_equal_string(rname, "CYT"))
            return 'C';
        if (is_equal_string(rname, "  G") || is_equal_string(rname, " DG") ||
            is_equal_string(rname, "GUA"))
            return 'G';
        if (is_equal_string(rname, "  T") || is_equal_string(rname, " DT") ||
            is_equal_string(rname, "THY"))
            return 'T';
        if (is_equal_string(rname, "  U") || is_equal_string(rname, " DU") ||
            is_equal_string(rname, "URA"))
            return 'U';
        fatal("unrecognized residue name: '%s'\n", rname);
        return 'X';
    }

    void
    PairFinder::_handle_str()
    {
        char parfile[BUF512], *ChainID, *bseq, **AtomName, **ResName, **Miscs;
        double **xyz;
        long num, num_residue, *ResSeq, *RY, **seidx;

        del_extension(args.pdbfile, parfile);
        /* read in the PDB file */
        num = number_of_atoms(args.pdbfile, args.hetatm, Gvars.misc_pars.alt_list);
        AtomName = cmatrix(1, num, 0, 4);
        ResName = cmatrix(1, num, 0, 3);
        ChainID = cvector(1, num);
        ResSeq = lvector(1, num);
        xyz = dmatrix(1, num, 1, 3);
        Miscs = cmatrix(1, num, 0, NMISC);
        read_pdb(args.pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
                 args.hetatm, Gvars.misc_pars.alt_list);

        /* get the numbering information of each residue */
        seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

        /* get base sequence, RY identification */
        bseq = cvector(1, num_residue);
        RY = lvector(1, num_residue);
        get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

        nt_info = cmatrix(1, num_residue, 0, BUF32);
        populate_nt_info(num_residue, seidx, ResName, ChainID, ResSeq, Miscs, bseq, nt_info);

        PairFinder::_duplex(num, num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
               ResSeq, Miscs, xyz, parfile, &Gvars.misc_pars);
        free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
        free_lmatrix(seidx, 1, num_residue, 1, 2);
        free_cvector(bseq, 1, num_residue);
        free_lvector(RY, 1, num_residue);
        free_cmatrix(nt_info, 1, num_residue, 0, BUF32);
    }

    void
    PairFinder::find_pair(X3dna::X3Basepairs &basepairs)
    {
        time_t time0;

        time(&time0);

        set_my_globals("find_pair");

        // frprintf(stderr, "\nhandling file <%s> \n \n", args.pdbfile);

        PairFinder::_handle_str();

        clear_my_globals();

        basepairs = bps;

    }
}