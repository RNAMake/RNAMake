#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wwritable-strings"

#ifndef __FIND_PAIR_H__
#define __FIND_PAIR_H__

#include <math/numerical.hpp>
#include <util/x3dna/x3dna.h>

namespace util {

class PairFinder {

public:
  explicit PairFinder(std::string pdb);

  void find_pair(util::x3dna::X3dna::X3Basepairs &basepairs);

  // private functions
private:
  void _handle_str();

  void _write_fpmst(double *morg, double *morien, FILE *rframe,
                    util::x3dna::X3dna::X3BPInfo *bp_info);

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

  std::vector<x3dna::X3dna::X3Basepair> _bps;

  std::map<std::pair<int, std::string>, math::Vector3> _atoms;

  struct Args {
    String pdbfile;
    char outfile[BUF512];
    long ds = 2;
    long divide = false;
    long hetatm = true;
    long pairs = false;
    long detailed = false;
    long waters = false;
  };

  math::Vector3s vectors;

  Args _args = Args();
};

} // namespace util

#endif
