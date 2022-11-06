//
// Created by Joe Yesselman on 6/23/22.
//

#include <base/log.hpp>
#include <util/io/pdb_parser.hpp>

namespace util::io {

void PDBParser::parse(const String & path) {
  if(!std::filesystem::exists(path)) {
    String msg = "file " + path + " does not exist";
    base::log_and_throw<base::InputException>(msg);
  }
  _in.open(path);
}

void PDBParser::next() {
  getline(_in, _line);
}

PDBLineType PDBParser::get_line_type() {
  if (_line.size() < 6) {
    return PDBLineType::UNKNOWN;
  }
  String startswith = _line.substr(0, 6);
  if (startswith == ("ATOM  ")) {
    return PDBLineType::ATOM;
  }
  else if (startswith == ("HETATM")) {
    return PDBLineType::HETATOM;
  }
  else if(startswith == "TER  ") {
    return PDBLineType::TER;
  }
  else if(startswith == "MODEL ") {
    return PDBLineType::MODEL;
  }
  return PDBLineType::UNKNOWN;

}

const PDBAtomLineData & PDBParser::get_line_atom_data() {
  try{
    _line_data.atom_id = std::stoi(_line.substr(6, 5));
    _line_data.atom_name = base::string::trim(_line.substr(12, 4));
    _line_data.atom_loc = _line[16];
    _line_data.res_name = base::string::trim(_line.substr(17, 4));
    _line_data.res_num = std::stoi(_line.substr(22, 4));
    _line_data.i_code = _line[26];
    _line_data.x = std::stod(_line.substr(30, 8));
    _line_data.y = std::stod(_line.substr(38, 8));
    _line_data.z = std::stod(_line.substr(46, 8));
    _line_data.occ = std::stod(_line.substr(54, 6));
    _line_data.b_factor = std::stod(_line.substr(60, 6));
  }
  catch(...) {
    throw base::InputException("invalid pdb line");
  }

  return _line_data;
}


bool PDBParser::done() {
  return !_in.good();
}

}
