#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wwritable-strings"

#ifndef __FIND_PAIR_H__
#define __FIND_PAIR_H__

#include <math/numerical.h>
#include <util/x3dna.h>
#include <util/x3dna_src.h>

namespace util {

class PairFinder {
 public:
  PairFinder(std::string pdb);

  void find_pair(X3dna::X3Basepairs &basepairs);

  // private functions
 private:
  void _handle_str();

  void _write_fpmst(double *morg, double *morien, FILE *rframe,
                    X3dna::X3BPInfo *bp_info);

  void _write_bestpairs(long num_bp, long **base_pairs, long *bp_idx,
                        char *bseq, long **seidx, char **AtomName,
                        char **ResName, char *ChainID, long *ResSeq,
                        char **Miscs, double **xyz, double **orien,
                        double **org, long **htm_water, miscPars *misc_par);

  void _duplex(long num, long num_residue, char *bseq, long **seidx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, char *parfile, miscPars *misc_pars);

  // Private variables
 private:
  char **nt_info;

  std::vector<X3dna::X3Basepair> bps;

  struct Args {
    char pdbfile[BUF512];
    char outfile[BUF512];
    long ds = 2;
    long curves = false;
    long curves_plus = false;
    long divide = false;
    long hetatm = true;
    long pairs = false;
    long detailed = false;
    long waters = false;
    long hjb = false;
  };

  Args args = Args();
};

}  // namespace util

#endif
