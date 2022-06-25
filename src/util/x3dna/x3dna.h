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
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <ctime>

// RNAMake Headers
#include <base/paths.hpp>
#include <base/string.hpp>
#include <math/matrix_3x3.hpp>
#include <math/numerical.hpp>
#include <math/vector_3.hpp>
#include <util/io/pdb_parser.hpp>

namespace util::x3dna {

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

#define DUMMY (-1L)

#define BUF32 32
#define BUF512 512
#define BUF1K 1024
#define BUFBIG 8192

#define UNUSED_PARAMETER(x) (void)(x)
#define NELEMS(x) ((sizeof(x)) / (sizeof((x)[0])))

/* ********* SLRE: http://code.google.com/p/slre/ */
enum slre_option { SLRE_CASE_SENSITIVE = 0, SLRE_CASE_INSENSITIVE = 1 };


struct miscPars{
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
};


struct struct_Gvars {
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
};

/* global variables declaration */
extern struct_Gvars Gvars;

#define DEBUG_LEVEL 6

#define SKIPS "#\0"      /* lines to be skipped */

#define UNKATM "XX"
#define SFACTOR 2L /* scale factor for increasing memory */
#define PI 3.141592653589793
#define XEPS 1.0e-7
#define XBIG 1.0e+18
#define XBIG_CUTOFF 1.0e+16
#define MFACTOR 10000.0
#define NMISC 34    /* number of characters of miscellaneous items */

#define FIG_BOUND 166    /* boundary offset: 10/72*1200 */

#define PAR_FILE "misc_3dna.par" /* miscellaneous parameters */
#define BASE_FILE "baselist.dat" /* 3-letter to 1-letter base residue */
#define ATOM_FILE "atomlist.dat" /* 2-letter to atomic symbol */

#define REF_FILE "ref_frames.dat" /* reference frames */
#define BESTP_FILE "bestpairs.pdb" /* best base-pairs */
#define BPORDER_FILE "bp_order.dat"   /* base-pair ordering in <find_pair> */
#define SEVEN_FILE "cf_7methods.par" /* compare seven methods in <analyze> */

#define BOND_UPPER_LIMIT 2.5          /* for function torsion */
#define HTWIST0 0.05                  /* minimum helical twist */
#define BOND_FACTOR 1.15              /* bond distance criterion */
#define NBOND_FNUM 2.0                /* estimated # of bond from # of atoms */
#define NUM_RESIDUE_ATOMS BUF512      /* max. no. of atoms in a residue */
#define EMPTY_NUMBER (-9999.99)
#define NELE 12              /* 12 elements */
#define O3P_UPPER 2.5        /* upper limit for O3-P connection */
#define RTNNUM 37            /* number of returned values from check_pair */
#define PSTNUM 29            /* number of parameters kept in pair_stat */
#define OLCRT 1.2            /* criterion for overlaps */
#define END_STACK_XANG 125.0 /* find_pair; bdl070 (105 deg) */

#define RA_LIST                                                                \
  " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
#define WC_LIST "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
#define CB_LIST "ACGITU"
#define NT_LIST                                                                \
  "  A", "  C", "  G", "  I", "  T", "  U", "ADE", "CYT", "GUA", "INO", "THY", \
      "URA", " +A", " +C", " +G", " +I", " +T", " +U"

#define WATER_LIST "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP"

#define OVERLAP 0.01

/* personal_begin */
/* for snap */
#define SNAP_OPTS "snap_options"
#define WITH_BASE 1 /* base or side-chain */
#define WITH_BKBN 2 /* backbone */

/* personal_end */

#endif /* _X3DNA_H */

/// @brief - x3dna_fnc.h

void populate_nt_info(long num_residue, long **seidx, char **ResName,
                      char *ChainID, long *ResSeq, char **Miscs, char *bseq,
                      char **nt_info);

long **read_input(char *inpfile, char *pdbfile, char *outfile, long *ds,
                  long *num_bp, long *ip, long *hetatm);

void print_header(long ds, long num_bp, long num, char *pdbfile, FILE *fp);

void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch,
                          double opening);

void get_bp_3char_symbols(long bp_type, char zdir, char *bp_sym);

void bpstep_par(double **rot1, double *org1, double **rot2, double *org2,
                double *pars, double **mst_orien, double *mst_org);

void helical_par(double **rot1, double *org1, double **rot2, double *org2,
                 double *pars, double **mst_orien, double *mst_org);

void print_par(char **bp_seq, long num_bp, long ich, long ishel, double **param,
               FILE *fp);

void output_ave_std(long num, double **parcln, int dnum, char *fmt, FILE *fp);


void write_mst(long ds, long num_bp, long **pair_num, char **bp_seq,
               double *mst_orien, double *mst_org, long **seidx,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               double **xyz, char **Miscs, long **htm_water,
               double **twist_rise, char *strfile);

long ratom_xyz(long *ratom_list, long only_ring, double **xyz, double *oave,
               double **oxyz);

void get_bp_zoave(long ia, long ib, double **orien, double **org, double *oave,
                  double *zave);

void ring_oidx(long num, long num_residue, long *RY, long **seidx,
               char **AtomName, double **xyz, long *idx, long **ring_atom);

void get_cntatom(long *ringlist, long **connect, long *idx);

double get_oarea(long r1, long r2, long **ring_atom, double *oave, double *zave,
                 double **xyz, long only_ring);

int case_strcmp(const char *s1, const char *s2);

int case_strncmp(const char *s1, const char *s2, long n);

char *case_strchr(const char *s, int c);

long is_equal_string(const char *str1, const char *str2);

long is_skip_line(char *line);

void bname_ext(char *src, char *ext, char *dst);

void get_tag_string_pair(char *prefix, char *tag, char *btag, char *etag);

void get_xml_tag(FILE *fpxml, char *prefix, char *line, char *connector,
                 char *tag, char *tag_str);

void get_xml_tag_long(FILE *fpxml, char *prefix, char *line, char *tag,
                      long *lval);

void get_xml_tag_double(FILE *fpxml, char *prefix, char *line, char *tag,
                        double *dval);

long tag_match(char *prefix, char *line, char *tag);

void base_frame(long num_residue, char *bseq, long **seidx, long *res_type,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, char *BDIR, double **orien,
                double **org);

void set_default_misc_pars(miscPars *misc_pars);

long read_PairInfo(char *inpfile, long **pair_info);

long is_linked(long i, long j, double **o3_p);

double distance_ab(double **o3_p, long ia, long ib, long ipa, long ipb);

void reverse_y_z_columns(double **R);

void set_my_globals(char *pgname);

void clear_my_globals();

void residue_strid(char chain_id, long res_seq, char *misc, char *rname,
                   char *idmsg);

void convert_resNameSpace(char *resName, char replacement, char *newName);

void snap_atype(char **AtomName, long num_residue, long **seidx, long *res_type,
                long **atom_cidx);

void cleanup_files(long renew, long cleanup);

void set_std_base_pdb(char *bdir, long irna, char bname, char *spdb);

char *my_getline(FILE *fp);

char *trim(char *a);

long itemize(char *str, char *item[], long itemsize);

long item_list(char *str, char *item[], long itemsize, char *sep_chars);

void refs_right_left(long bnum, double **orien, double **org, double **r1,
                     double *o1, double **r2, double *o2);

void mst2orien(double *orien_vec, long ioffset, double **mst);

void orien2mst(double *orien_vec, long ioffset, double **mst);

void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx);

void mtx_2_x_y_z(double **mtx, double *x, double *y, double *z);

void cehs_average(long inum_base, long *ivec, double **orien, double **org,
                  double **mst, double *morg);

void pair2mst(long inum_base, long *ivec, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, char **Miscs, double **xyz,
              double **orien, double **org, long **seidx, double *mst_orien,
              double *mst_org, long **htm_water, miscPars *misc_pars, FILE *fp);

FILE *open_file(char *filename, char *filemode);

long close_file(FILE *fp);

long exist_file(char *filename);

void remove_file(char *filename);

void rename_file(char *src, char *dst);

long upperstr(char *a);

long lowerstr(char *a);

void print_sep(FILE *fp, char x, long n);

void check_slash(char *BDIR);

void del_extension(char *fullname, char *okname);

void fatal(char *fmt, ...);

long number_of_atoms(const String &, const miscPars &);

long read_pdb(String pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs,
              long hetatm, char *ALT_LIST);

void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs);

void deduce_misc(char **Miscs, char **AtomName, long i, char *str);

long is_dna_with_backbone(long ib, long ie, char **AtomName);

void normalize_resName_atomName(long is_dna, const char *rname0,
                                const char *aname0, char *rname, char *aname);

void pdb_record(long ib, long ie, long *inum, long idx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, double **xyz,
                char **Miscs, FILE *fp);

void move_position(double **d, long nr, long nc, double *mpos);

long **residue_idx(long num, long *ResSeq, char **Miscs, char *ChainID,
                   char **ResName, long *num_residue);

long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib,
                   long ie);

void normalize_atom_symbol(char *asym);

void get_atomlist(char **atomlist, long *num_sa);

long has_atom_name(long ib, long ie, char **AtomName, char *aname);

void get_baselist(char **baselist, long *num_sb);

void get_seq(long num_residue, long **seidx, char **AtomName, char **ResName,
             char *ChainID, long *ResSeq, char **Miscs, double **xyz,
             char *bseq, long *RY);

long num_strmatch(char *str, char **strmat, long nb, long ne);

void get_idmsg(char *rname, char cid, long snum, char icode, char *idmsg);

long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg);

void get_BDIR(char *BDIR, char *filename);

void align2zaxis(long num, double *haxis, double **rotmat, double **xyz,
                 double **xyzH);

void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx);

double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz,
                  double **R, double *orgi);

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

void lsort(long n, long *a, long *idx);

void lreverse(long ia, long n, long *lvec);

void bring_atoms(long ib, long ie, long ra_num, char **AtomName, long *nmatch,
                 long *batom);

void all_bring_atoms(long num_residue, long *RY, long **seidx, char **AtomName,
                     long *num_ring, long **ring_atom);

void change_xyz(long side_view, double *morg, double **mst, long num,
                double **xyz);

void get_side_view(long ib, long ie, double **xyz);

long is_valid_base(char c, char *valid_bases);

long repeat_num();

char *read_repeat(char *crepeat, long fixed, char *valid_bases, long *nbp);

void pair_checking(long ip, long ds, long num_residue, char *pdbfile,
                   long *num_bp, long **pair_num);

void base_str(char chain_id, long res_seq, char *misc, char *rname, char bcode,
              long stnd, char *idmsg);

void hb_numlist(long i, long j, char basei, char basej, long **seidx, long *idx,
                char **AtomName, double **xyz, miscPars *misc_pars,
                long *num_hb, long **num_list);

long good_hbatoms(miscPars *misc_pars, char *atom1, char *atom2, long idx1,
                  long idx2);

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

long str_pmatch(char *str, char *sstr);

long case_str_pmatch(char *str, char *sstr);

double cvt2double(char *str);

long cvt2long(char *str);

long equalsign_pos(char *str);

long get_lvalue(char *str, long vmin, long vmax);

double get_dvalue(char *str, double vmin, double vmax);

void get_strvalue(char *str, char *dst, long expand_tilde);

double z1_z2_angle_in_0_to_90(double *z1, double *z2);

void print_frame(FILE *fp, double *O, double **R);

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

void dval_swap(double *pa, double *pb);

void lval_swap(long *pa, long *pb);

double dval_min(double a, double b);

long lval_min(long a, long b);

long lval_in_set(long lval, long ib, long ie, long *s);

long dval_in_range(double dval, double dlow, double dhigh);

long lval_in_range(long lval, long llow, long lhigh);

void ave_dmatrix(double **d, long nr, long nc, double *avedm);

void std_dmatrix(double **d, long nr, long nc, double *stddm);

void init_cmatrix(char **cmtx, long nrl, long nrh, long ncl, long nch,
                  char init_val);

void init_dmatrix(double **dmtx, long nrl, long nrh, long ncl, long nch,
                  double init_val);

void init_lmatrix(long **lmtx, long nrl, long nrh, long ncl, long nch,
                  long init_val);

void init_cvector(char *cvec, long ib, long ie, char init_val);

void init_dvector(double *dvec, long ib, long ie, double init_val);

void init_lvector(long *lvec, long ib, long ie, long init_val);

int dval_compare(const void *v1, const void *v2);


double p1p2_dist(double *xyz1, double *xyz2);

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

void multi_matrix(double **a, long nra, long nca, double **b, long nrb,
                  long ncb, double **o);

void transpose_matrix(double **a, long nr, long nc, double **o);

void identity_matrix(double **d, long n);

void xhelfunc(double *param, double **orien, double **mst, double *pos,
              double *mpos);
