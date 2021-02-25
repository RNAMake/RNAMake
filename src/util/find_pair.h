#ifndef __FIND_PAIR_H__
#define __FIND_PAIR_H__

#include <bits/stdc++.h>
#include <math/numerical.h>
#include <util/x3dna_src.h>
#include <util/x3dna.h>

namespace util
{

    class PairFinder
    {
    public:
        PairFinder(std::string pdb);

        std::vector<X3dna::X3Basepair>
        find_pair();
        
    private:
        struct Args
        {
            char pdbfile[BUF512];
            char outfile[BUF512];
            long ds = 2;
            long curves = FALSE;
            long curves_plus = FALSE;
            long divide = FALSE;
            long hetatm = TRUE;
            long pairs = FALSE;
            long detailed = FALSE;
            long waters = FALSE;
            long hjb = FALSE;
        };

        Args args = Args();

        // X3dna::X3BPInfo *bp_info;

        void
        handle_str();

        void
        write_fpmst(double *morg, double *morien, FILE *rframe, X3dna::X3BPInfo *bp_info);

        void
        write_bestpairs(long num_bp, long **base_pairs, long *bp_idx, char *bseq,
                                long **seidx, char **AtomName, char **ResName, char *ChainID,
                                long *ResSeq, char **Miscs, double **xyz, double **orien,
                                double **org, long **htm_water, miscPars *misc_par);

        void
        duplex(long num, long num_residue, char *bseq, long **seidx, long *RY,
                       char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                       char **Miscs, double **xyz, char *parfile,
                       miscPars *misc_pars);

        char **nt_info;

        std::vector<X3dna::X3Basepair> bps;
    };
}

#endif