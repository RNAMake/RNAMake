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
        PairFinder();

        std::vector<X3dna::X3Basepair>
        find_pair(String pdb);
        
    private:
        typedef struct
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
        } args;

        String bp_info;

        std::vector<X3dna::X3Basepair> bps;
    };
}

#endif