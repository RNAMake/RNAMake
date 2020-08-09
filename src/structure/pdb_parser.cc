//
//  pdb_parser.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

//RNAMake Headers
#include "base/string.h"
#include "base/file_io.h"
#include "math/xyz_vector.h"
#include "structure/pdb_parser.h"

namespace structure {

ResidueOPs const &
PDBParser::parse(
        String const & pdb_file,
        int protein,
        int rna,
        int others) {

    residues_ = ResidueOPs();

    Strings lines = base::get_lines_from_file(pdb_file);
    String startswith;
    String atomname, resname, resnum, chid, alt;
    math::Point coords;
    String sx, sy, sz;
    double x, y, z;

    Strings atomnames, resnames, chainids, icodes, resnums;
    math::Points coordinates;

    for (auto const & line : lines) {
        startswith = line.substr(0, 6);
        if (startswith.compare("ATOM  ") == 0 ||
            startswith.compare("HETATM") == 0) {
            atomname = line.substr(12, 4);
            atomname = base::trim(atomname);
            resname = line.substr(17, 4);
            resname = base::trim(resname);
            chid = line.substr(21, 1);
            alt = line.substr(16, 1);

            sx = line.substr(30, 8);
            sx = base::trim(sx);
            sy = line.substr(38, 8);
            sy = base::trim(sy);
            sz = line.substr(46, 8);
            sz = base::trim(sz);


            x = std::stod(sx);
            y = std::stod(sy);
            z = std::stod(sz);
            coords = math::Point(x, y, z);

            resnum = line.substr(22, 4);
            resnum = base::trim(resnum);
            atomnames.push_back(atomname);
            resnames.push_back(resname);
            chainids.push_back(chid);
            resnums.push_back(resnum);
            icodes.push_back(line.substr(26, 1));
            coordinates.push_back(coords);


        } else if (startswith.compare("ENDMDL") == 0 || startswith.substr(0, 3).compare("END") == 0) {
            break;
        }
    }

    String key;
    std::map<String, AtomOPs> residue_atoms;
    int already_has = 0;

    for (int i = 0; i < atomnames.size(); i++) {
        if (resnames[i].compare("HOH") == 0) { continue; }
        key = resnames[i] + " " + resnums[i] + " " + chainids[i] + " " + icodes[i];
        if (residue_atoms.find(key) == residue_atoms.end()) {
            residue_atoms[key] = AtomOPs();
        }
        already_has = 0;
        for (auto const & a : residue_atoms[key]) {
            if (a->name().compare(atomnames[i]) == 0) {
                //already_has = 1;
                //break;
            }
        }
        if (already_has) { continue; }
        residue_atoms[key].push_back(AtomOP(new Atom(atomnames[i], coordinates[i])));
    }
    ResidueType rtype;
    Strings spl;
    String icode = "";
    ResidueOP r;
    for (auto & kv : residue_atoms) {
        if (kv.second.size() < 6) { continue; }
        spl = base::split_str_by_delimiter(kv.first, " ");
        if (!rts_.contains_rtype(spl[0]) && !others) { continue; }
        if (!others) {
            rtype = rts_.get_rtype_by_resname(spl[0]);
        } else {
            auto atom_names = Strings();
            for (auto const & a : kv.second) { atom_names.push_back(a->name()); }
            rtype = _get_new_residue_type(spl[0], atom_names);
        }

        if (!protein && rtype.set_type() == SetType::PROTEIN) { continue; }
        if (!rna && rtype.set_type() == SetType::RNA) { continue; }


        icode = "";
        if (spl.size() > 3) { icode = spl[3]; }
        r = ResidueOP(new Residue(rtype, spl[0], std::stoi(spl[1]), spl[2], icode));
        r->setup_atoms(kv.second);
        residues_.push_back(r);

    }

    return residues_;
}

//adapted from new code to allow for small molecule residue types
ResidueType
PDBParser::_get_new_residue_type(
        String const & res_name,
        Strings const & atom_names) {

    auto atom_name_map = StringIntMap();
    int i = 0;
    for (auto const & name : atom_names) {
        atom_name_map[name] = i;
        i++;
    }
    auto extra_alt_names = Strings();

    return ResidueType(res_name, atom_name_map, SetType::UNKNOWN);
}

}
























