//
//  x3dna.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__x3dna__
#define __RNAMake__x3dna__

// standard headers
#include <cstring>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// RNAMake Headers
#include <base/paths.hpp>
#include <base/string.hpp>
#include <math/matrix_3x3.hpp>
#include <math/numerical.hpp>
#include <math/vector_3.hpp>

namespace util {
namespace x3dna {

// lowercase x3dna is the namespace
// uppercase X3dna is the object class

class X3dnaException : public std::runtime_error {
public:
  X3dnaException(String const &message) : std::runtime_error(message) {}
};

//- = U
//+ = P
//. = D
enum class X3dnaBPType {
  cmU = 0,   // cm-
  cMUM = 1,  // cM-M
  tWPW = 2,  // tW+W
  cDPM = 3,  // c.+M
  DWPW = 4,  //.W+W
  tWUM = 5,  // tW-M
  tmUM = 6,  // tm-M
  cWPM = 7,  // cW+M
  DWUW = 8,  //.W-W
  cMPD = 9,  // cM+.
  cDUm = 10, // c.-m
  cMPW = 11, // cM+W
  tMPm = 12, // tM+m
  tMUW = 13, // tM-W
  cmUm = 14, // cm-m
  cMUW = 15, // cM-W
  cWUW = 16, // cW-W
  cDUM = 17, // c.-M
  cmPM = 18, // cm+M
  cmUM = 19, // cm-M
  DDDD = 20, //....
  cmUW = 21, // cm-W
  tMUm = 22, // tM-m
  cDUW = 23, // c.-W
  cMPm = 24, // cM+m
  cMUm = 25, // cM-m
  cDDD = 26, // c...
  tWPm = 27, // tW+m
  cDPm = 28, // c.+m
  tmPm = 29, // tm+m
  tWPD = 30, // tW+.
  tmPW = 31, // tm+W
  tDDD = 32, // t...
  cWUD = 33, // cW-.
  cWUM = 34, // cW-M
  tDUW = 35, // t.-W
  tMPM = 36, // tM+M
  tDUM = 37, // t.-M
  cMUD = 38, // cM-.
  cWUm = 39, // cW-m
  tDPm = 40, // t.+m
  tMUD = 41, // tM-.
  cmPW = 42, // cm+W
  cMPM = 43, // cM+M
  cmPD = 44, // cm+.
  cmUD = 45, // cm-.
  cDUD = 46, // c.-.
  cWPW = 47, // cW+W
  tDUD = 48, // t.-.
  tDPW = 49, // t.+W
  tmUm = 50, // tm-m
  cWPD = 51, // cW+.
  tmPD = 52, // tm+.
  tDPD = 53, // t.+.
  cDPD = 54, // c.+.
  tDUm = 55, // t.-m
  tDPM = 56, // t.+M
  // Added by CJ
  tWUD = 57, // tW-.
  tmUW = 58, // tm-W
  tMUM = 59, // tM-M
  tMPD = 60, // tM+.
  cDPW = 61, // c.+W
  tmPM = 62, // tm+M
  tWUm = 63, // tW-m
  cWPm = 64, // cW+m
  tmUD = 65, // tm-.
  tWPM = 66, // tW+m
  DWPm = 67, //.W+m
  tMPW = 68, // tM+W
  DDPm = 69, //..+m
  tWUW = 70, // tW-W
  cmPm = 71, // cm+m
  DWUm = 72, //.W-m
  DMPm = 73, //.M+m
  DWPM = 74, //.W+M
  DMPM = 75, //.M+M
  DmPW = 76, //.m+W
  DWUM = 77, //.W+M
  DmPm = 78, //.m+m
  DDUM = 79, //..-M
  DMUm = 80, //.M-m
  DDUm = 81, //..-m
  DMPW = 82, //.M+W
  DMPD = 83, //.M+.
  DMUM = 84, //.M-M
  DmUm = 85, //.m-m
  DMUW = 86, //.M-W
  DWUD = 87, //.W-.
};

class X3dna {
public:
  X3dna();

  ~X3dna() {
    // delete s_;
  }

public:
  static void set_envs() {
    // check x3dna path
    auto os_name = getprogname();
    auto x3dna_path = base::path::resources_path() + "/x3dna/" + os_name + "/";
    String env = "X3DNA=" + x3dna_path;
    auto s = strdup(env.c_str());
    putenv(s);
  }

public:
  struct X3Residue {
    inline X3Residue(int nnum = -1, char nchain_id = ' ', char ni_code = ' ')
        : num(nnum), chain_id(nchain_id), i_code(ni_code) {}

    inline bool operator==(X3Residue const &r) const {
      if (num == r.num && chain_id == r.chain_id && i_code == r.i_code) {
        return 1;
      } else {
        return 0;
      }
    }
    // added by CJ ... designed to chaeck that all data fields were filled out
    // appropriately
    bool valid() const {
      return num != std::numeric_limits<int>::min() &&
             chain_id != ' '; // don't include ni_code bc it's usually none...
                              // change this?
    }
    int num;
    char chain_id, i_code;
  };

  typedef std::vector<X3Residue> X3Residues;

  struct X3Basepair {
    X3Residue res1, res2;
    math::Vector3 d;
    math::Matrix3x3 r;
    X3dnaBPType bp_type;
    // added by CJ... designed to check that all data fields were filled out
    // appropriately
    bool valid() const {
      return res1.valid() && res1.valid() &&
             !roughly_equal(r, math::Matrix3x3{-1., -1., -1., -1., -1., -1.,
                                               -1., -1., -1.}) &&
             !roughly_equal(d, math::Vector3{-1., -1., -1.});
    }
    String to_string() const;

    unsigned int key() const {
      const auto start = std::min(res1.num, res2.num);
      const auto end = std::max(res1.num, res2.num);
      return start * 256 * 256 + end * 256 + static_cast<int>(bp_type);
    }
  };

  typedef std::vector<X3Basepair> X3Basepairs;

  struct X3Motif {
    X3Residues residues;
    String mtype;
  };

  typedef std::vector<X3Motif> X3Motifs;

  struct X3BPInfo {
    inline X3BPInfo(std::smatch const &match) {
      auto i = -1;
      for (auto const &v : match) {
        i++;
        if (i == 0) {
          continue;
        } else if (i == 1) {
          res1_chain_id = String(v)[0];
        } else if (i == 2) {
          res1_num = std::stoi(String(v));
        } else if (i == 3) {
          res1_type_name = String(v);
        } else if (i == 4) {
          res1_name = String(v)[0];
        } else if (i == 5) {
          res2_chain_id = String(v)[0];
        } else if (i == 6) {
          res2_num = std::stoi(String(v));
        } else if (i == 7) {
          res2_type_name = String(v);
        } else if (i == 8) {
          res2_name = String(v)[0];
        }
      }
      if (i != 8) {
        throw X3dnaException("could not properly parse X3dna ref frame line:");
      }
    }
    String res1_type_name, res2_type_name;
    char res1_name, res2_name;
    char res1_chain_id, res2_chain_id;
    int res1_num, res2_num;
  };

public:
  X3Basepairs get_basepairs(String const &) const;

  X3Basepairs get_basepairs_json(String const &) const;

  X3Motifs get_motifs(String const &) const;

public:
  void set_rebuild_files(bool rebuild_files) const {
    rebuild_files_ = rebuild_files;
  }

private:
  void _generate_ref_frame(String const &) const;

  void generate_dssr_file(String const &) const;

private:
  void _delete_files(Strings const &) const;

  void _delete_file(String const &) const;

  math::Vector3 _convert_string_to_point(String const &) const;

  void _parse_ref_frame_file(String const &, X3Basepairs &) const;

  std::map<String, Strings>
  _parse_dssr_file_into_sections(String const &) const;

  Strings _split_over_white_space(String const &) const;

  X3Residue *_parse_dssr_res_str(String const &) const;

  X3Motifs _parse_dssr_section(Strings const &, String const &) const;

  X3Motifs _parse_dssr_helix_section(Strings const &) const;

private:
  String bin_path_;
  char *s_;
  Strings ref_frame_files_to_delete_, dssr_files_to_delete_;
  // flags to decide whether to build files or not
  mutable bool rebuild_files_;
  mutable bool generated_ref_frames_, generated_dssr_;
  mutable bool no_ref_frames_;
};

X3dnaBPType get_x3dna_by_type(String const &);

String get_str_from_x3dna_type(X3dnaBPType);

String compare_bps(X3dna::X3Basepairs &, X3dna::X3Basepairs &);

void json_cleanup();
}
} // namespace util

template <typename T1, typename T2, typename T3>

struct triplet {
  T1 first;
  T2 second;
  T3 third;

  bool operator<(const triplet &other) const {
    if (first == other.first) {
      if (second == other.second) {
        return third < other.third;
      } else {
        return second < other.second;
      }
    } else {
      return first < other.first;
    }
  }
};

#endif /* defined(__RNAMake__x3dna__) */

/// @brief - x3dna_src.h

#ifndef _X3DNA_H
#define _X3DNA_H

#define NR_END 1 /* for NRC related functions */
#define FREE_ARG char *

#define NO_MATCH -1L
#define DUMMY -1L

#define BUF32 32
#define BUF512 512
#define BUF1K 1024
#define BUF2K 2048
#define BUFBIG 8192

#define UNUSED_PARAMETER(x) (void)(x)
#define NELEMS(x) ((sizeof(x)) / (sizeof((x)[0])))

/* ********* SLRE: http://code.google.com/p/slre/ */
enum slre_option { SLRE_CASE_SENSITIVE = 0, SLRE_CASE_INSENSITIVE = 1 };
enum slre_capture { SLRE_STRING, SLRE_INT, SLRE_FLOAT };


typedef struct {
  long min_base_hb;
  double hb_lower;
  double hb_dist1;
  double hb_dist2;
  char hb_atoms[BUF512];
  long hb_idx[BUF512];
  char alt_list[BUF512];

  double max_dorg;
  double min_dorg;
  double max_dv;
  double min_dv;
  double max_plane_angle;
  double min_plane_angle;
  double max_dNN;
  double min_dNN;

  double helix_break;
  double std_curved;

  double water_dist;
  double water_dlow;
  char water_atoms[BUF512];
  long water_idx[BUF512];

  double o3p_dist;
} miscPars;

typedef struct {
  std::vector<math::Vector3> vect;
  std::vector<String> info_vect;
} bp_vectors;

typedef struct {
  double origin[3];
  double x_axis[3];
  double y_axis[3];
  double z_axis[3];
} base_pair;

typedef struct {
  long DEBUG;
  long VERBOSE;
  long NUM_ELE;
  long CHAIN_CASE;
  long ALL_MODEL;
  long ATTACH_RESIDUE;   /* following Pascal's request */
  long THREE_LETTER_NTS; /* output 3-letter nucleotides: ADE/CYT/GUA/THY/URA */
  long PDBV3;            /* PDB v3.x; OP1/OP2/C7; DA/DC/DG/DT etc */
  long ORIGINAL_COORDINATE; /* use original PDB coordinates in bestpairs.pdb etc
                             */
  long OCCUPANCY;           /* if to check for atom occupancy property */
  long HEADER;              /* how much header info to output */

  long mmcif;       /* output structure in PDBx/mmCIF format */
  double NT_CUTOFF; /* RMSD cutoff for identify a nucleotide */

  char X3DNA_VER[BUF512];
  char X3DNA_HOMEDIR[BUF512];
  char CHAIN_MARKERS[BUF512];
  char REBUILD_CHAIN_IDS[BUF512];

  char *PROGNAME;
  char **ATOM_NAMES;

  long NUM_SATOM;
  char **ATOMLIST;

  long NUM_SBASE;
  char **BASELIST;

  char **AtomName0;
  char **ResName0;
  long Name0;

  long label_RC8_YC6;

  miscPars misc_pars;
} struct_Gvars;

/* global variables declaration */
extern struct_Gvars Gvars;

#define DEBUG_LEVEL 6

#define SEPC '\t'        /* tab-delimited fields */
#define WSPACES " ,\t\n" /* space, comma, tab & newline */
#define SKIPS "#\0"      /* lines to be skipped */

#define UNKATM "XX"
#define DUMSTR "XXXXXX"
#define SFACTOR 2L /* scale factor for increasing memory */
#define PI 3.141592653589793
#define XEPS 1.0e-7
#define XBIG 1.0e+18
#define XBIG_CUTOFF 1.0e+16
#define MFACTOR 10000.0
#define NMISC 34    /* number of characters of miscellaneous items */
#define NATOMCOL 11 /* number of atom with defined colors */
#define NBASECOL 7  /* number of base with defined colors */

#define PS_DFTSIZE 500   /* 500 points PS default size */
#define PS_BOUND 10      /* boundary offset */
#define FIG_DFTSIZE 8333 /* 8333 units XFIG default size */
#define FIG_BOUND 166    /* boundary offset: 10/72*1200 */

#define PAR_FILE "misc_3dna.par" /* miscellaneous parameters */
#define BASE_FILE "baselist.dat" /* 3-letter to 1-letter base residue */
#define ATOM_FILE "atomlist.dat" /* 2-letter to atomic symbol */
#define HELP3DNA "help3dna.dat"  /* help file for 3DNA */

#define REF_FILE "ref_frames.dat" /* reference frames */
#define MREF_FILE                                                              \
  "mref_frames.dat" /* multiplet reference frames in <find_pair> */
#define POC_FILE                                                               \
  "poc_haxis.r3d"                  /* P, O4' and C1' atom radius in r3d format \
                                    */
#define MUL_FILE "multiplets.pdb"  /* multiplets */
#define ALLP_FILE "allpairs.pdb"   /* all base-pairs */
#define BESTP_FILE "bestpairs.pdb" /* best base-pairs */
#define STACK_FILE "stacking.pdb"  /* for stacking diagram */
#define HSTACK_FILE                                                            \
  "hstacking.pdb" /* for stacking w.r.t. to middle helical frame */
#define HLXREG_FILE "hel_regions.pdb" /* helical regions in <find_pair> */
#define BPORDER_FILE "bp_order.dat"   /* base-pair ordering in <find_pair> */
#define COLCHN_FILE "col_chains.scr"  /* color chains in <find_pair> */
#define COLHLX_FILE "col_helices.scr" /* color helices in <find_pair> */
#define AUX_FILE "auxiliary.par"      /* auxiliary parameters in <analyze> */
#define BPSTEP_FILE "bp_step.par" /* base-pair step parameters in <analyze> */
#define HLXSTEP_FILE                                                           \
  "bp_helical.par"                   /* helical step parameters in <analyze>   \
                                      */
#define SEVEN_FILE "cf_7methods.par" /* compare seven methods in <analyze> */
#define HB_FILE                                                                \
  "hbonds_info.dat" /* H-bonding information <r3d_atom/stack2img> */
#define LKG_FILE "bonds_lkg.dat"    /* bond linkage, as in HB_FILE */
#define SNUM_FILE "serial_num.pdb"  /* PDB file with serial atom numbers */
#define ATOMALC_FILE "atom_lkg.alc" /* atomic linkage <r3d_atom/stack2img> */
#define BBLKALC_FILE "bblk_lkg.alc" /* base block linkage <pdb2img> */
#define TMP_FILE "tmp_file"    /* temporary multiplets input file <find_pair> */
#define MULBP_FILE "mulbp.inp" /* multiplets input file <find_pair> */
#define ROTMAT_FILE "rotmat.dat" /* rotation matrix <rotate_mol> */
#define VIEW1_FILE "pmiview1"    /* PMI view #1 file <rotate_mol> */
#define VIEW2_FILE "pmiview2"    /* PMI view #2 file <rotate_mol> */
#define VIEW3_FILE "pmiview3"    /* PMI view #3 file <rotate_mol> */

#define NP 101L                       /* maximum number of pairs per base */
#define BOND_UPPER_LIMIT 2.5          /* for function torsion */
#define HTWIST0 0.05                  /* minimum helical twist */
#define BOND_FACTOR 1.15              /* bond distance criterion */
#define NBOND_FNUM 2.0                /* estimated # of bond from # of atoms */
#define NON_WC_IDX 6                  /* non-Watson-Crick base index */
#define AXIS_LENGTH 3.5               /* reference axis length */
#define NUM_BASE_ATOMS BUF512         /* max. no. of base atoms in a residue */
#define NUM_RESIDUE_ATOMS BUF512      /* max. no. of atoms in a residue */
#define NUM_DINUCLEOTIDE_ATOMS BUFBIG /* max. no. of atoms per dinucleotide */
#define EMPTY_NUMBER -9999.99
#define EMPTY_CRITERION -9999
#define MAXBASE 30000        /* maximum number of bases in regular/fiber */
#define NELE 12              /* 12 elements */
#define O3P_UPPER 2.5        /* upper limit for O3-P connection */
#define RTNNUM 37            /* number of returned values from check_pair */
#define PSTNUM 29            /* number of parameters kept in pair_stat */
#define OLCRT 1.2            /* criterion for overlaps */
#define MBASES 50            /* maximum bases per unit */
#define MAXCH 100            /* maximum # of chains in a biological unit */
#define END_STACK_XANG 125.0 /* find_pair; bdl070 (105 deg) */
#define MAXCLEN 52           /* Maximum length of a color name */

#define RA_LIST                                                                \
  " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
#define WC_LIST "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
#define CX_LIST "ACGITUX"
#define CB_LIST "ACGITU"
#define NT_LIST                                                                \
  "  A", "  C", "  G", "  I", "  T", "  U", "ADE", "CYT", "GUA", "INO", "THY", \
      "URA", " +A", " +C", " +G", " +I", " +T", " +U"

#define WATER_LIST "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP"

#define NUM_SAA 20 /* number of standard amino acids */
#define NUM_ATM 12 /* number of side chain atoms, excluding Ca */
#define AA_LIST                                                                \
  "ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU", /* hydrophobic */           \
      "ASP", "GLU", "LYS", "ARG",                  /* charged */               \
      "SER", "THR", "TYR", "HIS", "CYS", "ASN", "GLN", "TRP", /* polar */      \
      "GLY"                                                   /* smallest */

#define OVERLAP 0.01

#define PBLKALC_FILE "pblk_lkg.alc" /* peptide block linkage */
#define SNAP_PEP_PDB "snap_pep.pdb"
#define SNAP_PEP_ALC "snap_pep.alc"
#define TRSP_RMS 0.25
#define DUMCHAR '@'

/* personal_begin */
/* for snap */
#define LBIG 1000000
#define DNA_BASE "ACGT"
#define SNAP_AAA "snap_aaa.pdb"
#define ABLKALC_FILE                                                           \
  "ablk_lkg.alc" /* amino acid planar moiety block linkage                     \
                  */
#define SNAP_NTS "snap_nts.pdb"
#define SNAP_OPTS "snap_options"
#define WITH_BASE 1 /* base or side-chain */
#define WITH_BKBN 2 /* backbone */
#define WITH_BOTH 3

#define PDBX "PDBx:"
#define O2_STACK "o2_stack.dat"
/* personal_end */

#endif /* _X3DNA_H */

/// @brief - x3dna_fnc.h

void populate_nt_info(long num_residue, long **seidx, char **ResName,
                      char *ChainID, long *ResSeq, char **Miscs, char *bseq,
                      char **nt_info);

void populate_nt_list(long num_residue, long **seidx, long *RY, char *bseq,
                      char **AtomName, double **xyz, long **nt_list);

long **read_input(char *inpfile, char *pdbfile, char *outfile, long *ds,
                  long *num_bp, long *ip, long *hetatm);

void print_header(long ds, long num_bp, long num, char *pdbfile, FILE *fp);

void output_Borg_P_C1_C4(long num_residue, double **org, double **xyz,
                         long **nt_list, char **nt_info);

void atom_list(long ds, long num_bp, long **pair_num, long **seidx, long *RY,
               char **bp_seq, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, long **phos, long **c6_c8,
               long **sugar, long **chi);

void get_nt_torsion(long num_residue, double **org, double **xyz,
                    long **nt_list, double **nt_torsion);

void get_ss_Zp_Dp(long num_residue, double **org, double **orien, double **xyz,
                  long **nt_list, double **ss_Zp_Dp);

void output_nt_torsion(long num_residue, char **nt_info, long **nt_list,
                       double **nt_torsion, double **ss_Zp_Dp, FILE *fp);

void get_nt_bb_torsion(double **nt_bb_torsion, long num_residue, long **seidx,
                       long *RY, char **AtomName, char **ResName, char *ChainID,
                       long *ResSeq, char **Miscs, double **xyz);

void backbone_torsion(long ds, long num_bp, long **pair_num, char **bp_seq,
                      long **sugar, long **chi, double **xyz,
                      double **nt_bb_torsion, FILE *fp);

void p_c1_dist(long ds, long num_bp, char **bp_seq, long **phos, long **chi,
               double **xyz, long *bphlx, FILE *fp);

void lambda_d3(long num_bp, char **bp_seq, long **chi, long **c6_c8,
               double **xyz, FILE *fp);

void print_axyz(long num_bp, char **bp_seq, long **aidx, char *aname,
                double **xyz);

void groove_width(long parallel, long num_bp, char **bp_seq, long **phos,
                  double **xyz, long *bphlx, FILE *fp);

void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch,
                          double opening);

void set_chain_nmarkers019_to_symbols(long num, long *nmarkers, char *cmarkers);

void get_bp_3char_symbols(long bp_type, char zdir, char *bp_sym);

void ref_frames(long ds, long num_bp, long **pair_num, char **bp_seq,
                long **seidx, long *RY, char **AtomName, char **ResName,
                char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                FILE *fp, double **orien, double **org, long *WC_info,
                long *str_type, long irna, long **o3p_brk);

void bpstep_par(double **rot1, double *org1, double **rot2, double *org2,
                double *pars, double **mst_orien, double *mst_org);

void helical_par(double **rot1, double *org1, double **rot2, double *org2,
                 double *pars, double **mst_orien, double *mst_org);

void print_par(char **bp_seq, long num_bp, long ich, long ishel, double **param,
               FILE *fp);

void output_ave_std(long num, double **parcln, int dnum, char *fmt, FILE *fp);

void prt_stepstr(char **step_str, long num_step, long *bphlx, long ishel,
                 double **param, FILE *fp);

void prt_step_par(char **bp_seq, long num_bp, long *bphlx, long ishel,
                  double **param, FILE *fp);

void bz_check(double **r1, double *o1, double **r2, double *o2, long bz,
              long *bz_junction, long *z_step);

void get_mtwist(long nbpm1, long *bphlx, long *WC_info, double **twist_rise,
                double *twist_p, double *twist_n);

void get_parameters(long ds, long num_bp, char **bp_seq, double **orien,
                    double **org, long *WC_info, FILE *fp, double **twist_rise,
                    double *mst_orien, double *mst_org, double *mst_orienH,
                    double *mst_orgH, long *bphlx, long istart, long istep,
                    long bz, long *str_type, long **pair_num, char **nt_info);

void parvec2mtx(double *parvec, long num, double **parmtx);

void print_ss_rebuild_pars(double **pars, long num_bp, char *str, char **bp_seq,
                           FILE *fp);

void print_ds_rebuild_pars(double **bp_par, double **step_par, long num_bp,
                           char *str, char **bp_seq, FILE *fp);

void print_ref(char **bp_seq, long num_item, long ich, double *org,
               double *orien, FILE *fp);

void write_mst(long ds, long num_bp, long **pair_num, char **bp_seq,
               double *mst_orien, double *mst_org, long **seidx,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               double **xyz, char **Miscs, long **htm_water,
               double **twist_rise, char *strfile);

void print_xyzP(long parallel, long nbpm1, char **bp_seq, long **phos,
                double *mst_orien, double *mst_org, double **xyz, FILE *fp,
                char *title_str, double **aveP, long p_offset);

void print_PP(long parallel, double **twist_rise, long num_bp, char **bp_seq,
              long **phos, double *mst_orien, double *mst_org,
              double *mst_orienH, double *mst_orgH, double **xyz, long *WC_info,
              long *bphlx, long abi, long **chi, FILE *fp);

void str_classify(double twist_p, double twist_n, long str_type, long parallel,
                  long num_bp, FILE *fp);

double a_hlxdist(long idx, double **xyz, double *hlx_axis, double *hlx_pos);

void print_radius(char **bp_seq, long nbpm1, long ich, double **p_radius,
                  double **o4_radius, double **c1_radius, long *bphlx,
                  FILE *fp);

void helix_radius(long ds, long num_bp, char **bp_seq, double **orien,
                  double **org, long **phos, long **chi, double **xyz,
                  long *bphlx, FILE *fp);

void print_shlx(char **bp_seq, long nbpm1, long ich, double *shlx_orien,
                double *shlx_org, FILE *fp);

void get_helix_axis(long ds, long num_bp, char **bp_seq, double **orien,
                    double **org, long *bphlx, FILE *fp);

void get_axis(long nvec, long **idx, long num, double **xyz, long nb, long *C1b,
              long *C1e, double *std_rise, double *hrise, double *haxis,
              double *hstart, double *hend);

void print_poc_r3d(double *rave, double *hstart, double *hend);

void global_analysis(long ds, long num_bp, long num, char **bp_seq, long **chi,
                     long **phos, double **xyz, FILE *fp);

void base_overlap(long ds, long num_bp, long num, long num_residue,
                  long **pair_num, long *bRY, char **bp_seq, long **seidx,
                  char **AtomName, double **xyz, long *idx, double **orien,
                  double **org, FILE *fp);

long ratom_xyz(long *ratom_list, long only_ring, double **xyz, double *oave,
               double **oxyz);

void get_zoave(long istep, long ds, double **orien, double **org, double *oave,
               double *zave);

void get_bp_zoave(long ia, long ib, double **orien, double **org, double *oave,
                  double *zave);

void ring_oidx(long num, long num_residue, long *RY, long **seidx,
               char **AtomName, double **xyz, long *idx, long **ring_atom);

void get_cntatom(long *ringlist, long **connect, long *idx);

double get_oarea(long r1, long r2, long **ring_atom, double *oave, double *zave,
                 double **xyz, long only_ring);

void verify_oarea(void);

void cehs_base_atoms(char **AtomName, long ib, long ie, long *num_batom,
                     long *batom);

void cehs_bppar(double **rot1, double *org1, double **rot2, double *org2,
                double *pars, double **mst_orien, double *mst_org);

void cehs_pars(long num_bp, long istart, long istep, long **pair_num,
               char **bp_seq, long **seidx, long **c6_c8, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, double *bp_orien, double *bp_org,
               FILE *fp);

void schnaap_global(long num_bp, long num, char **bp_seq, long **chi,
                    double **xyz, double *bp_orien, double *bp_org, FILE *fp);

void out_cehs(long num_bp, char **bp_seq, long *bphlx, double **orien,
              double **org, FILE *fp);

void compdna(double **rot1, double *org1, double **rot2, double *org2,
             double *pars, double **mst_orien, double *mst_org);

void out_compdna(long num_bp, char **bp_seq, long *bphlx, double **orien,
                 double **org, FILE *fp);

void my_curves(double **rot1, double *org1, double **rot2, double *org2,
               double *pars);

void curves_mbt(long ibp, double **orien, double **org, double **cvr,
                double *cvo);

void out_curves(long num_bp, char **bp_seq, long *bphlx, double **orien,
                double **org, FILE *fp);

void freehelix(double **rot1, double *org1, double **rot2, double *org2,
               double *pars, double **mst_orien, double *mst_org);

void out_freehelix(long num_bp, char **bp_seq, long *bphlx, double **orien,
                   double **org, FILE *fp);

void sgl_helix(double **rot1, double **rot2, double *rot_ang, double *rot_hlx);

void ngeom(double **rot1, double *org1, double **rot2, double *org2,
           double *pars, double **mst_orien, double *mst_org);

void out_ngeom(long num_bp, char **bp_seq, long *bphlx, double **orien,
               double **org, FILE *fp);

void nuparm(double **rot1, double *org1, double **rot2, double *org2,
            double *pars, double **mst_orien, double *mst_org, double *hpars,
            long get_hpar);

void out_nuparm(long num_bp, char **bp_seq, long *bphlx, double **orien,
                double **org, FILE *fp);

void rna(double **rot1, double *org1, double *pvt1, double **rot2, double *org2,
         double *pvt2, double *pars, double **mst_orien, double *mst_org);

void pvt_dxdy(double **rot1, double *org1, double *pvt1, double *pars,
              double **mst_orien, double *mst_org);

void out_rna(long ds, long num_bp, char **bp_seq, long *bphlx, double **orien,
             double **org, FILE *fp);

void other_pars(long num_bp, char **bp_seq, long *bphlx, double **orien,
                double **org);

/* analyze.c */

/* anyhelix.c */

/* app_fncs.c */

long string_contains_only_those_characters(char *str, char *chars_set);

void bpid_wc_str(long bpid, double zdir, char *wc);

int case_strcmp(const char *s1, const char *s2);

int case_strncmp(const char *s1, const char *s2, long n);

char *case_strstr(const char *haystack, const char *needle);

char *case_strchr(const char *s, int c);

void kbd_input(char *msg);

long get_line_number(char *filename, long skips);

long is_empty_string(const char *str);

long is_equal_string(const char *str1, const char *str2);

long is_equal_case_string(const char *str1, const char *str2);

void null_line_comment(char *str);

long is_comment_line(char *line);

long is_empty_line(char *line);

long is_skip_line(char *line);

void bname_ext(char *src, char *ext, char *dst);

double get_point2line_perp_distance(double *pnt, double *line_p1,
                                    double *line_p2);

void get_tag_string_pair(char *prefix, char *tag, char *btag, char *etag);

void get_xml_tag(FILE *fpxml, char *prefix, char *line, char *connector,
                 char *tag, char *tag_str);

void print_xml_tag(FILE *fpxml, char *prefix, char *line, char *tag, char *otag,
                   FILE *fp);

void get_xml_tag_long(FILE *fpxml, char *prefix, char *line, char *tag,
                      long *lval);

void get_xml_tag_double(FILE *fpxml, char *prefix, char *line, char *tag,
                        double *dval);

long tag_match(char *prefix, char *line, char *tag);

void extract_attribute(char *line, char *attr, char *attr_val, long to_lower);

void print_ptable();

long set_switch_default_true(char *option);

void check_required_file(char *filename, char *invalid_str, char *msg);

void define_frame_by_3atoms(long num, double **xyz, double **refmat);

long get_xmlArgNumber(char *str, char *pstr);

void base_frame(long num_residue, char *bseq, long **seidx, long *res_type,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, char *BDIR, double **orien,
                double **org);

void peptide_info(long num_residue, long **seidx, char **AtomName,
                  char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                  double **xyz, long *res_type, long *cidx, long *mchain);

void base_blks(long num_residue, long *res_type, double **orien, double **org,
               char *bseq, char *BDIR, char *alcfile);

void set_default_misc_pars(miscPars *misc_pars);

void lsplane_xyz(double **xyz, long num_plane, long *atom_plane, double **nxyz,
                 double *z);

long read_PairInfo(char *inpfile, long **pair_info);

long is_linked_by_gap(long i, long j, double **o3_p);

long is_linked(long i, long j, double **o3_p);

double distance_ab(double **o3_p, long ia, long ib, long ipa, long ipb);

void write_rotmat(double **rotmat);

void read_rotmat(char *rotfile, double **rotmat);

void reverse_y_z_columns(double **R);

long get_num_nt(long num_residue, long *RY);

void peptide_frame(long num_residue, char *BDIR, long *res_type, long *mchain,
                   double **xyz, double **orien, double **org);

void peptide_blks(long num_residue, char *BDIR, long *cidx, double **orien,
                  double **org, char *alcfile);

void set_my_globals(char *pgname);

void clear_my_globals(void);

long check_global_options(char *option);

void get_AA_frames(long num_residue, long **seidx, long *res_type,
                   char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                   char **Miscs, double **xyz, char *BDIR, double **orien,
                   double **org);

void verify_Cb_coordinates(long nmatch, double **ePxyz, double **fitted_xyz);

void residue_chain_resnum(char chain_id, long res_seq, char *misc, char *idmsg);

void residue_strid(char chain_id, long res_seq, char *misc, char *rname,
                   char *idmsg);

void convert_resNameSpace(char *resName, char replacement, char *newName);

void get_planarFrame(char *aa, long pnum, char *planar_atoms[], char *BDIR,
                     char *idmsg, long ib, long ie, char **AtomName,
                     double **xyz, double *orien_i, double *org_i);

void get_argFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_pheFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_tyrFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_trpFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_hisFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_lysFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_asnFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_glnFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_aspFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_gluFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);

void get_planarAA_frames(long num_residue, long **seidx, char **AtomName,
                         char **ResName, char *ChainID, long *ResSeq,
                         char **Miscs, double **xyz, long *res_type,
                         long *cidx);

void planarAA_blks(long num_residue, long **seidx, char **ResName, char *BDIR,
                   long *cidx, double **paa_orien, double **paa_org,
                   char *alcfile);

void print_resid(long num_residue, long **seidx, char **ResName, char *ChainID,
                 long *ResSeq, char **Miscs, long *res_type);

void snap_atype(char **AtomName, long num_residue, long **seidx, long *res_type,
                long **atom_cidx);

long number_of_aa(long num_residue, long *res_type);

long number_of_nt(long num_residue, long *res_type);

void get_snap_par(double **rot1, double *org1, double **rot2, double *org2,
                  char *direction, double *trs_dist, double *rot_dist,
                  double *pars, double *orgP, double **rotP);

void write_snap_par(FILE *fp, char *direction, double dist, double rot_ang,
                    double *pars, double *orgP, double **rotP);

void write_atom_xyz(FILE *fp, char *fmt, double *xyz, double dft);

long set2frame(long inum, long *ivec, long **seidx, double **xyz, double *morg,
               double **mst, long *serial, double **xyz_pair);

long set2frameCa(long inum, long *ivec, long **seidx, long *res_type,
                 char **AtomName, double **xyz, double *morg, double **mst,
                 long *serial, double **xyz_pair);

long get_nextModelNumber(char *pdbfile);

void output_naa_str(char *filename, char *fmode, char *idmsg, long num,
                    long *serial, char **AtomName, char **ResName,
                    char *ChainID, long *ResSeq, double **xyz_pair,
                    char **Miscs, long out_org);

void cleanup_files(long renew, long cleanup);

/* cehs.c */

/* cmn_fncs.c */
long set_3letter_base_pdb(char *res_name, char *spdb);

void set_std_base_pdb(char *bdir, long irna, char bname, char *spdb);

void set_std_base_pdb00(char *bdir, long irna, char bname, char *spdb);

void print_used_time(time_t time0);

void parcat(char *str, double par, char *format, char *bstr);

void print_bp_crit(miscPars *misc_pars, FILE *fp);

char *my_getline(FILE *fp);

long csplit(char *str, char *item[], long itemsize, char sepc);

char *trim(char *a);

char *ltrim(char *a);

char *rtrim(char *a);

long itemize(char *str, char *item[], long itemsize);

long item_list(char *str, char *item[], long itemsize, char *sep_chars);

void refs_right_left(long bnum, double **orien, double **org, double **r1,
                     double *o1, double **r2, double *o2);

void refs_i_j(long b1, long b2, double *bp_orien, double *bp_org, double **r1,
              double *o1, double **r2, double *o2);

void ref_frame_i(long bnum, double *bp_orien, double *bp_org, double **r,
                 double *o);

void mst2orien(double *orien_vec, long ioffset, double **mst);

void orien2mst(double *orien_vec, long ioffset, double **mst);

void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx);

void mtx_2_x_y_z(double **mtx, double *x, double *y, double *z);

void cehs_average(long inum_base, long *ivec, double **orien, double **org,
                  double **mst, double *morg);

void geom_average(long inum_base, long *ivec, double **orien, double **org,
                  double **mst, double *morg);

void pair2mst(long inum_base, long *ivec, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, char **Miscs, double **xyz,
              double **orien, double **org, long **seidx, double *mst_orien,
              double *mst_org, long **htm_water, miscPars *misc_pars, FILE *fp);

void get_chi_angle(long num_residue, long *RY, char *bseq, long **seidx,
                   double **xyz, char **AtomName, char **ResName, char *ChainID,
                   long *ResSeq, char **Miscs, double *chi, long **idxCN);

FILE *open_tmpfile(void);

FILE *open_file(char *filename, char *filemode);

long close_file(FILE *fp);

long exist_file(char *filename);

void remove_file(char *filename);

void rename_file(char *src, char *dst);

void copy_file_pointer(FILE *fpi, FILE *fpo, char *msg);

void cpcat_file(char *src, char *dst, char *method);

long upperstr(char *a);

long lowerstr(char *a);

char *my_strdup(const char *src);

void print_sep(FILE *fp, char x, long n);

void check_slash(char *BDIR);

void delete_end_slash(char *str);

long lround(double d);

void del_extension(char *fullname, char *okname);

void bname_noext(char *src, char *dst);

void fatal(char *fmt, ...);

void print_pdb_title(char *pdbfile, char *chain_list, FILE *fp);

long number_of_atoms(char *pdbfile, long hetatm, char *ALT_LIST);

long read_pdb(char *pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs,
              long hetatm, char *ALT_LIST);

void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs);

void reset_xyz(long num, double **xyz, char *fmt);

void deduce_misc(char **Miscs, char **AtomName, long i, char *str);

long is_dna_with_backbone(long ib, long ie, char **AtomName);

void normalize_resName_atomName(long is_dna, const char *rname0,
                                const char *aname0, char *rname, char *aname);

void pdb_record(long ib, long ie, long *inum, long idx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, double **xyz,
                char **Miscs, FILE *fp);

void write_pdb(long num, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, double **xyz, char **Miscs, char *pdbfile);

void write_mmcif(long num, char **AtomName, char **ResName, char *ChainID,
                 long *ResSeq, double **xyz, char *pdbfile);

void write_pdbcnt(long num, char **AtomName, char **ResName, char *ChainID,
                  long *ResSeq, double **xyz, long **connect, char *pdbfile);

void move_position(double **d, long nr, long nc, double *mpos);

long **residue_idx(long num, long *ResSeq, char **Miscs, char *ChainID,
                   char **ResName, long *num_residue);

long frag_contain_metal(long ib, long ie, long *is_metal);

void atom_metal(long num_atoms, char **AtomName, long *is_metal);

void residue_wtype(long num_residue, long **seidx, char **ResName,
                   char **AtomName, double **xyz, char **Miscs, long *res_type,
                   long only_ntaa);

long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib,
                   long ie);

void normalize_atom_symbol(char *asym);

void get_atomlist(char **atomlist, long *num_sa);

long has_atom_name(long ib, long ie, char **AtomName, char *aname);

void get_baselist(char **baselist, long *num_sb);

void get_seq(long num_residue, long **seidx, char **AtomName, char **ResName,
             char *ChainID, long *ResSeq, char **Miscs, double **xyz,
             char *bseq, long *RY);

void get_bpseq(long ds, long num_bp, long **pair_num, long **seidx,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, char **bp_seq, long *RY);

long strmatch_idx(char *str, char **strmat, long nb, long ne);

long num_strmatch(char *str, char **strmat, long nb, long ne);

void get_idmsg(char *rname, char cid, long snum, char icode, char *idmsg);

long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg);

double torsion(double **d);

double torsion2(double **d);

void get_BDIR(char *BDIR, char *filename);

void align2zaxis(long num, double *haxis, double **rotmat, double **xyz,
                 double **xyzH);

void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx);

double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz,
                  double **R, double *orgi);

void ls_plane(double **bxyz, long n, double *pnormal, double *ppos,
              double *odist, double *adist);

void arb_rotation(double *va, double ang_deg, double **rot_mtx);

double vec_ang(double *va, double *vb, double *vref);

void get_vector(double *va, double *vref, double deg_ang, double *vo);

void rotate(double **a, long i, long j, long k, long l, double *g, double *h,
            double s, double tau);

void eigsrt(double *d, double **v, long n);

void jacobi(double **a, long n, double *d, double **v);

void dludcmp(double **a, long n, long *indx, double *d);

void dlubksb(double **a, long n, long *indx, double *b);

void dinverse(double **a, long n, double **y);

void rotx(double ang_deg, double **rotmat);

void roty(double ang_deg, double **rotmat);

void rotz(double ang_deg, double **rotmat);

void get_alc_nums(char *alcname, long *num, long *nbond);

void read_alc(char *alcname, long *num, long *nbond, char **AtomName,
              double **xyz, long *ibase, long **linkage);

void write_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
               long **linkage, char *alcfile);

void free_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
              long zero_1, long **linkage);

void cnct_org(long num_bp, long ia, long ib, char **tAtomName, double **txyz,
              long *tibase, long **tlinkage, double **org_xyz);

void dsort(long n, double *a, long *idx);

void lsort(long n, long *a, long *idx);

void lreverse(long ia, long n, long *lvec);

void fig_title(FILE *fp);

void ps_title_cmds(FILE *fp, char *imgfile, long *bbox);

void get_fig_xy(long num, double **xyz, long nO, double **oxyz, long *urxy,
                long frame_box, FILE *fp);

void get_pxy(double *xy1, double *xy2, double r, double *px, double *py);

void alc2fig(long nobj, long *idx, long *depth, long **allobj, double **blkxyz,
             double **oxyz, long *ibase, long faces[][5], long *opts, FILE *fp);

void get_ps_xy(char *imgfile, long *urxy, long frame_box, FILE *fp);

void alc2ps(long nobj, long *idx, long **allobj, double **blkxyz, double **oxyz,
            long *ibase, long faces[][5], long *opts, FILE *fp);

void bring_atoms(long ib, long ie, long ra_num, char **AtomName, long *nmatch,
                 long *batom);

void all_bring_atoms(long num_residue, long *RY, long **seidx, char **AtomName,
                     long *num_ring, long **ring_atom);

void base_idx(long num, char *bseq, long *ibase, long single);

long basepair_idx(char *bpi);

void plane_xyz(long num, double **xyz, double *ppos, double *nml,
               double **nxyz);

void prj2plane(long num, long ra_num, char **AtomName, double **xyz, double z0,
               double **nxyz);

void adjust_xy(long num, double **xyz, long nO, double **oxyz,
               double scale_factor, long default_size, long *urxy);

void get_depth(long nobj, long *zval, long *depth);

void raster3d_header(long num, double **xyz, double scale_factor,
                     long no_header, long frame_box, FILE *fp);

void get_r3dpars(double **base_col, double *hb_col, double *width3,
                 double **atom_col, char *label_style);

void r3d_rod(long itype, double *xyz1, double *xyz2, double rad, double *rgbv,
             FILE *fp);

void r3d_dash(double *xyz1, double *xyz2, double hb_width, double *hb_col,
              FILE *fp);

void r3d_sphere(double *xyz1, double rad, double *rgbv, FILE *fp);

void cpk_model(long num, long *idx, double **xyz, double ballrad,
               double **colrgb, FILE *fp);

void r3d_tripln(long itype, double *xyz1, double *xyz2, double *xyz3,
                double *rgbv, FILE *fp);

void r3d_block_edge(double *rgbv, long ioffset8, double **blkxyz, double w1,
                    FILE *fp);

void base_label(double **rxyz, char *label_style, double *rgbv, char *bname_num,
                FILE *fp);

void fill_base_ring(long num_residue, long num_ring, long **ring_atom,
                    double **xyz, long *ibase, char *bseq, double **base_col,
                    char *label_style, long label_ring, long *ResSeq, FILE *fp);

void process_alc(char *alcfile, char *imgfile, double scale_factor, long *opts);

void alc_3images(long *opts, long nobj, long num_blk, long num_blk8,
                 long nO_lkg, long nO, double **oxyz, double **blkxyz,
                 long *blkibase, long **linkage, double scale_factor,
                 char *imgfile);

void get_alc_objs(long num_blk, double **blkxyz, long nO, double **oxyz,
                  long nO_lkg, long **linkage, long faces[][5], long **allobj);

void get_fig_pars(double *dot_sep, long *dlcol, long *dwidth, long *bp1width,
                  long *bp2width, long **bc_idx, double *msat, double *Msat,
                  long *o_sides, long *line_width, long *join_style,
                  long *cap_style, long *mfcol);

void frame_xyz(long side_view, double *morg, double **mst, long num,
               double **xyz);

void change_xyz(long side_view, double *morg, double **mst, long num,
                double **xyz);

void get_side_view(long ib, long ie, double **xyz);

void get_CNidx(long ds, long num_bp, long **chi, long **idx, long *nvec,
               long *C1b, long *C1e);

void add_3axes(long *num, char **AtomName, long *ibase, double **xyz,
               long *nbond, long **linkage, long side_view, double axis_len);

char *get_sequence(char *Wbase, long *num_bp);

char **single2double(long nbp, char *bseq, char *Wbase, char *Cbase);

long is_valid_base(char c, char *valid_bases);

long repeat_num(void);

char *read_sequence(char *seqfile, char *valid_bases, long *nbp);

char *read_repeat(char *crepeat, long fixed, char *valid_bases, long *nbp);

void combine_pstnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                    char **tAtomName, char **tResName, char *tChainID,
                    long *tResSeq, double **txyz, char **tAtomName2,
                    double **txyz2);

void reverse_stnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                   char **tAtomName, char **tResName, char *tChainID,
                   long *tResSeq, double **txyz, char **tAtomName2,
                   double **txyz2, long basep);

void pair_checking(long ip, long ds, long num_residue, char *pdbfile,
                   long *num_bp, long **pair_num);

void double_print_msg(char *msg, FILE *fp);

void drct_checking(long ds, long num_bp, long **pair_num, long **seidx,
                   char **AtomName, double **xyz, long *parallel, long *bbexist,
                   long **o3p_brk, FILE *fp);

void residue_idstr(char chain_id, long res_seq, char *rname, char *idmsg);

void base_str(char chain_id, long res_seq, char *misc, char *rname, char bcode,
              long stnd, char *idmsg);

void write_lkglist(long nbond, long **linkage, char **AtomName, char **ResName,
                   char *ChainID, long *ResSeq, char **Miscs);

void hbond_info(long **pair_num, char *bseq, long **seidx, long *idx,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, long *RY, long *num_hbond,
                long **hb_linkage);

void hbond_pdb(long num, long num_residue, char *bseq, long **seidx, long *idx,
               long *RY, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, double **xyz, long *num_hbond,
               long **hb_linkage, long pwise);

void hbond_list(long i, long j, char **AtomName, char **ResName, char *ChainID,
                long *ResSeq, double **xyz, char **Miscs, char *bseq,
                long **seidx, long *idx, long **hb_linkage, miscPars *misc_pars,
                long **num_list, long *num_hbond, long ilayer, FILE *fph);

void hb_numlist(long i, long j, char basei, char basej, long **seidx, long *idx,
                char **AtomName, double **xyz, miscPars *misc_pars,
                long *num_hb, long **num_list);

void hb_information(long num_bp, long **pair_num, char **bp_seq, long **seidx,
                    long *idx, char **AtomName, double **xyz, long *WC_info,
                    FILE *fp);

long good_hbatoms(miscPars *misc_pars, char *atom1, char *atom2, long idx1,
                  long idx2);

void read_lkginfo(char *lkgfile, long num, long *nbond, long **linkage);

void read_hbinfo(char *hbfile, long num, long *num_hbond, long **hb_linkage);

void update_hb_idx(long idx, double *dtmp, long *ddidx, double *hb_dist,
                   long cur_idx);

void hb_atompair(long num_hbonds, char **hb_atom1, char **hb_atom2,
                 double *hb_dist, long *lkg_type, miscPars *misc_pars);

long validate_hbonds(long num_hbonds, double *hb_dist, long *lkg_type,
                     char *hb_type, char basei, char basej, char **hb_atom1,
                     char **hb_atom2);

void get_hbond_ij(long i, long j, char basei, char basej, miscPars *misc_pars,
                  long **seidx, long *idx, char **AtomName, double **xyz,
                  char *hb_info);

char donor_acceptor(char basei, char basej, char *hb_atom1, char *hb_atom2);

long asym_idx(char *asym, char atoms_list[NELE][3], long dft_lval);

void atom_info(long idx, char atoms_list[NELE][3], double *covalence_radii,
               double *vdw_radii);

void aname2asym(const char *aname0, char *my_asym, long num_sa,
                char **atomlist);

void atom_idx(long num, char **AtomName, char **Miscs, long *idx);

void get_bonds(long num, char **AtomName, double **xyz, long num_residue,
               long *RY, long **seidx, long **connect);

void atom_linkage(long ib, long ie, long *idx, double **xyz, char **Miscs,
                  char *ChainID, long nbond_estimated, long *nbond,
                  long **linkage);

void lkg2connect(char **AtomName, long ib, long ie, long nbond, long **linkage,
                 long **connect);

void init_htm_water(long waters, long num, long num_residue, long *idx,
                    long **htm_water);

void identify_htw(long num_residue, long **seidx, long *RY, char **AtomName,
                  char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                  double **xyz, long **htm_water);

long attached_residues(long inum_base, long *ivec, long *ivec2, long **seidx,
                       double **xyz, long **htm_water, miscPars *misc_pars);

void print_pairinfo(long i, long j, char basei, char basej, double *rtn_val,
                    double *chi, miscPars *misc_pars, long **seidx, long *idx,
                    char **AtomName, double **xyz, char *bseq, long detailed,
                    FILE *fp);

void check_pair(long i, long j, char *bseq, long **seidx, double **xyz,
                double **NC1xyz, double **orien, double **org, long *idx,
                char **AtomName, miscPars *misc_pars, double *rtn_val,
                long *bpid, long **ring_atom, long network);

void o3_p_xyz(long ib, long ie, char *aname, char **AtomName, double **xyz,
              double *o3_or_p, long idx);

long is_baseatom(char *atomname);

void base_info(long num_residue, char *bseq, long **seidx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, double **orien, double **org,
               double **NC1xyz, double **o3_p);

void help3dna_usage(char *program_name);

void help3dna(char *program_name);

void delH_pdbfile(char *inpfile, char *outfile);

void contact_msg(long prt_msg);

long str_pmatch(char *str, char *sstr);

long case_str_pmatch(char *str, char *sstr);

long is_numeric(char *str);

double cvt2double(char *str);

long cvt2long(char *str);

long equalsign_pos(char *str);

long get_lvalue(char *str, long vmin, long vmax);

double get_dvalue(char *str, double vmin, double vmax);

void get_strvalue(char *str, char *dst, long expand_tilde);

void reverse_string(char *str);

void cvtstr_set1toc2(char *str, char *set1, char c2);

void cvtstr_c1toc2(char *str, char c1, char c2);

double z1_z2_angle_in_0_to_90(double *z1, double *z2);

void do_nothing(void);

void skip_lines(long num, FILE *fp);

void check_havefile(char *filename, char *msg);

void print_frame(FILE *fp, double *O, double **R);

/* comb_str.c */

/* ex_str.c */

/* fiber.c */

/* find_pair.c */

/* fncs_slre.c */

const char *slre_match(enum slre_option options, const char *re,
                       const char *buf, int buf_len, ...);

int lux_match(enum slre_option options, const char *re, const char *buf);

int lux_bcmatch(const char *buf, const char *re);

int lux_ncmatch(const char *buf, const char *re);

/* frame_mol.c */

/* get_part.c */

/* mutate_bases.c */

/* nrutil.c */

void nrerror(char *error_text);

void vector_boundary_check(long nl, long nh, char *fun_name);

void matrix_boundary_check(long nrl, long nrh, long ncl, long nch,
                           char *fun_name);

char *cvector(long nl, long nh);

char *cvector_nr(long nl, long nh);

double *dvector(long nl, long nh);

double *dvector_nr(long nl, long nh);

long *lvector(long nl, long nh);

long *lvector_nr(long nl, long nh);

char **cmatrix(long nrl, long nrh, long ncl, long nch);

char **cmatrix_nr(long nrl, long nrh, long ncl, long nch);

double **dmatrix(long nrl, long nrh, long ncl, long nch);

double **dmatrix_nr(long nrl, long nrh, long ncl, long nch);

long **lmatrix(long nrl, long nrh, long ncl, long nch);

long **lmatrix_nr(long nrl, long nrh, long ncl, long nch);

void free_cvector(char *v, long nl, long nh);

void free_dvector(double *v, long nl, long nh);

void free_lvector(long *v, long nl, long nh);

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);

double dval_sqr(double dval);

void dval_swap(double *pa, double *pb);

void lval_swap(long *pa, long *pb);

void cval_swap(char *pa, char *pb);

double dval_max(double a, double b);

double dval_min(double a, double b);

long lval_max(long a, long b);

long lval_min(long a, long b);

double abs_dval_diff(double a, double b);

long lval_in_set(long lval, long ib, long ie, long *s);

long dval_in_range(double dval, double dlow, double dhigh);

long lval_in_range(long lval, long llow, long lhigh);

void max_dmatrix(double **d, long nr, long nc, double *maxdm);

void min_dmatrix(double **d, long nr, long nc, double *mindm);

void ave_dmatrix(double **d, long nr, long nc, double *avedm);

void std_dmatrix(double **d, long nr, long nc, double *stddm);

double max_dvector(double *d, long nl, long nh);

double min_dvector(double *d, long nl, long nh);

double ave_dvector(double *d, long n);

double std_dvector(double *d, long n);

void init_cmatrix(char **cmtx, long nrl, long nrh, long ncl, long nch,
                  char init_val);

void init_dmatrix(double **dmtx, long nrl, long nrh, long ncl, long nch,
                  double init_val);

void init_lmatrix(long **lmtx, long nrl, long nrh, long ncl, long nch,
                  long init_val);

void init_cvector(char *cvec, long ib, long ie, char init_val);

void init_cvector_all(char *cvec, long ib, long ie, char init_val);

void init_dvector(double *dvec, long ib, long ie, double init_val);

void init_lvector(long *lvec, long ib, long ie, long init_val);

void copy_dvector(double *d, double *s, long nl, long nh);

void copy_lvector(long *d, long *s, long nl, long nh);

int dval_compare(const void *v1, const void *v2);

int lval_compare(const void *v1, const void *v2);

int cstr_compare(const void *v1, const void *v2);

void negate_xyz(double *xyz1);

double p1p2_dist(double *xyz1, double *xyz2);

void p1p2_ave(double *xyz1, double *xyz2, double *ave);

long within_limits(double *xyz1, double *xyz2, double dlow, double dhigh);

void sumxyz(double *xyz1, double *xyz2, double *sxyz);

void avexyz(double *xyz1, double *xyz2, double *mxyz);

void ddxyz(double *xyz1, double *xyz2, double *dxyz);

void cpxyz(double *xyz1, double *xyz2);

void vec_orth(double *va, double *vref);

double dot(double *va, double *vb);

void cross(double *va, double *vb, double *vc);

long sign_control(double *va, double *vb, double *vref);

double veclen(double *va);

void vec_norm(double *va);

double dot2ang(double dotval);

double magang(double *va, double *vb);

double rad2deg(double ang);

double deg2rad(double ang);

void copy_dmatrix(double **a, long nr, long nc, double **o);

void copy_lmatrix(long **a, long nr, long nc, long **o);

void multi_matrix(double **a, long nra, long nca, double **b, long nrb,
                  long ncb, double **o);

void multi_vec_matrix(double *a, long n, double **b, long nr, long nc,
                      double *o);

void multi_vec_Tmatrix(double *a, long n, double **b, long nr, long nc,
                       double *o);

void transpose_matrix(double **a, long nr, long nc, double **o);

void identity_matrix(double **d, long n);

/* o1p_o2p.c */

/* pdb2img.c */

/* r3d_atom.c */

/* reb_fncs.c */

void link_o3_p(long num_residue, long **seidx, char **AtomName, double **xyz,
               char *ChainID, long **connect);

void atom_lkg(long num, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, double **xyz, char *outfile);

void atomic_pdb1(long num_bp, long num_atoms, long num_max_per_residue,
                 long is_helical, long xdir, char **bp_seq, double **step_par,
                 char *BDIR, char *outfile);

void atomic_pdb2(long parallel, long num_bp, long num_atoms,
                 long num_max_per_residue, long is_helical, long xdir,
                 char **bp_seq, double **bp_par, double **step_par, char *BDIR,
                 char *outfile);

void base_c1_atoms(char **AtomName, long ib, long ie, long *num_batom,
                   long *batom);

void extract_base_atoms(long num, char **AtomName, double **xyz, long *bnum,
                        char **bAtomName, double **bxyz);

void atomic_base_p(long parallel, long num_bp, long num_atoms,
                   long num_max_per_residue, long is_helical, long xdir,
                   char **bp_seq, double **bp_par, double **step_par,
                   char *BDIR, long *pidx, char *outfile);

void num_PDB_atoms(long num_bp, long is_single, char **bp_seq, char *BDIR,
                   long *num_atoms, long *num_max_per_residue);

void set_bp_pdb(long num_max_per_residue, char *fname1, char *fname2,
                double *param, long xdir, long *num1, long *num2,
                char **AtomName1, char **AtomName2, double **xyz1,
                double **xyz2, char ap);

void base_fname(char bname, char *BDIR, char *fname);

void block_alc1(long num_bp, long is_single, long is_helical, long xdir,
                char **bp_seq, double **step_par, char *BDIR, char *outfile);

void block_alc2(long num_bp, long is_helical, long xdir, char **bp_seq,
                double **bp_par, double **step_par, char *BDIR, char *outfile);

void set_bp_alc(char *fname1, char *fname2, double *param, long xdir, long num,
                long nbond, double **bp_xyz, char ap);

void xbpfunc(double *param, double **orien, double **mst, double *pos,
             double *mpos);

void xhelfunc(double *param, double **orien, double **mst, double *pos,
              double *mpos);

void print_ref_frames(FILE *fp, double *pos_next, double **orien_next);

/* rebuild.c */

/* regular_dna.c */

/* rotate_mol.c */

/* stack2img.c */

/* std_base.c */

/* step_hel.c */
