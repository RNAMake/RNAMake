//
// Created by Joe Yesselman on 6/23/22.
//

#ifndef RNAMAKE_SRC_UTIL_IO_PDB_PARSER_HPP_
#define RNAMAKE_SRC_UTIL_IO_PDB_PARSER_HPP_

#include <fstream>

#include <base/types.hpp>
#include <math/vector_3.hpp>

namespace util::io {

enum class PDBLineType{
  ATOM, HETATOM, MODEL, ENDMDL, END, TER, UNKNOWN
};

struct PDBAtomLineData {
public:
  int atom_id;
  String atom_name;
  char atom_loc;
  String res_name;
  int res_num;
  char i_code;
  double x;
  double y;
  double z;
  double occ;
  double b_factor;
};

/*
 * PDBS are column based with fixed indices for columns
 * 1-6:   "ATOM  ";
 * 7-11: ATOM ID
 * 13-16: Atom Name;
 * 17: Location indicator
 * 18-20: resname;
 * 22: chain identifier
 * 23-26: resnum;
 * 27: icode
 * 31-38: X,real(8.3)
 * 39-46: Y,real(8.3)
 * 47-54: Z,real(8.3)
 * 55-60: Occupancy
 * 61-66: TempFactor
 * 73-76: segID
 * 77-78: element
 * 79-80: Charge
 */
class PDBParser {
public:
  ~PDBParser() {
    if(_in.is_open()) { _in.close(); }
  }
public:
  void parse(const String &);

  void next();

  PDBLineType get_line_type();

  const PDBAtomLineData & get_line_atom_data();

  bool done();

private:

  PDBAtomLineData _line_data;
  std::ifstream _in;
  String _line;
};

}

#endif // RNAMAKE_SRC_UTIL_IO_PDB_PARSER_HPP_
