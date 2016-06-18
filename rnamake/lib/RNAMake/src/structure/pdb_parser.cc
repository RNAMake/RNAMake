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

ResidueOPs const &
PDBParser::parse(
    String const & pdb_file,
    int protein,
    int rna) {
    
    residues_ = ResidueOPs();
    
    Strings lines = get_lines_from_file(pdb_file);
    String startswith;
    String atomname, resname, resnum, chid, alt;
    Point coords;
    String sx, sy, sz;
    double x, y, z;
    
    Strings atomnames, resnames, chainids, icodes, resnums;
    Points coordinates;
    
    for(auto const & line : lines) {
        startswith = line.substr(0,6);
        if(startswith.compare("ATOM  ") == 0 ||
           startswith.compare("HETATM") == 0) {
            atomname = line.substr(12, 4);
            atomname = trim(atomname);
            resname  = line.substr(17, 4);
            resname  = trim(resname);
            chid     = line.substr(21, 1);
            alt      = line.substr(16, 1);
            
            sx       = line.substr(30, 8);
            sx       = trim(sx);
            sy       = line.substr(38, 8);
            sy       = trim(sy);
            sz       = line.substr(46, 8);
            sz       = trim(sz);
            
            
            x        = std::stod(sx);
            y        = std::stod(sy);
            z        = std::stod(sz);
            coords   = Point(x, y, z);
            
            resnum  = line.substr(22, 4);
            resnum  = trim(resnum);
            
            atomnames.push_back(atomname);
            resnames.push_back(resname);
            chainids.push_back(chid);
            resnums.push_back(resnum);
            icodes.push_back(line.substr(26, 1));
            coordinates.push_back(coords);
            
            
        }
        
        else if(startswith.compare("ENDMDL") == 0 || startswith.substr(0,3).compare("END") == 0) {
            break;
        }
    }
    
    String key;
    std::map<String, AtomOPs> residue_atoms;
    int already_has = 0;
    
    for(int i = 0; i < atomnames.size(); i++) {
        if(resnames[i].compare("HOH") == 0) { continue; }
        key = resnames[i] + " " + resnums[i] + " " + chainids[i] + " " + icodes[i];
        if(residue_atoms.find(key) == residue_atoms.end()) {
            residue_atoms[key] = AtomOPs();
        }
        already_has = 0;
        for(auto const & a : residue_atoms[key]) {
            if(a->name().compare(atomnames[i]) == 0) {
                //already_has = 1;
                //break;
            }
        }
        if(already_has) { continue;}
        residue_atoms[key].push_back(AtomOP(new Atom(atomnames[i], coordinates[i])));
    }
    
    
    ResidueType rtype;
    Strings spl;
    String icode = "";
    ResidueOP r;
    for(auto & kv : residue_atoms) {
        if(kv.second.size() < 6) { continue; }
        spl = split_str_by_delimiter(kv.first, " ");
        if(!rts_.contains_rtype(spl[0])) { continue; }
        rtype = rts_.get_rtype_by_resname(spl[0]);
        
        if(!protein && rtype.set_type() == SetType::PROTEIN) { continue; }
        if(!rna && rtype.set_type() == SetType::RNA) { continue; }
 
        icode = "";
        if(spl.size() > 3) { icode = spl[3]; }
        r = ResidueOP(new Residue(rtype, spl[0], std::stoi(spl[1]), spl[2], icode));
        r->setup_atoms(kv.second);
        residues_.push_back(r);
        
    }
    
    return residues_;
}



























