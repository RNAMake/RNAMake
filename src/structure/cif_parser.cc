//
//  pdb_parser.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include <sstream>
#include <set>

//RNAMake Headers
#include "base/string.h"
#include "base/file_io.h"
#include "math/xyz_vector.h"
#include "structure/cif_parser.h"

//    Residue(
//            ResidueType const & rtype,
//            String const & name,
//            int const & num,
//            String const & chain_id,
//            String const & i_code) :

namespace structure {

std::unordered_set<int>
CIFParser::_get_entity_ids(
                String const & header_contents
                    ) {
    
    // helper method that finds the entity id's that contain RNA
    auto ss= std::stringstream(header_contents); 
    auto line = String{};
    auto headers = std::map<String,int>(); 
    auto entity_contents = String{};
    auto index(0);
    //first we need to map the header names to the index
    while(std::getline(ss,line,'\n')) {
        if(line[0] == '_' && line.find(".") != std::string::npos) {
            const auto key = base::trim(base::split_str_by_delimiter(line,".")[1]);
            headers[key] = index++;
             
        } else if(line.find("loop_") == std::string::npos){
            entity_contents += line;
        }
    }
    // sometimes there will be ";"'s separating tokens.. replace it with whitespace
    std::replace(entity_contents.begin(),entity_contents.end(),';',' ');
    
    // tokenize the contents and see how long they are
    auto tks = base::tokenize_line(entity_contents);
    const auto num_entries = tks.size()/headers.size() ;
    
    // check that the sizing is correct
    assert(num_entries*tks.size() == entity_contents.size());
   
    // find the indicies we need
    const auto res_type_index = headers.at("type");
    const auto id_index = headers.at("entity_id");
   
    // loop through and add the indices if they are RNA, add the index
    auto allowed_ids = std::unordered_set<int>{};
    for(auto index = 0; index < num_entries; ++ index) {
        if( tks[index*headers.size()+res_type_index].find("polyribonucleotide") != std::string::npos) {
            allowed_ids.insert( 
                            std::stoi(tks[index*headers.size()+id_index])
                                );
        }
    }

    return allowed_ids;
}

ResidueOPs const &
CIFParser::parse(
        String const & cif_file,
        int protein,
        int rna,
        int others) {
    
    auto infile = std::ifstream(cif_file);
    auto line = String{}; 
    auto raw_contents = String{}; 
    
    auto entity_info = String{}; 

    while(std::getline(infile,raw_contents,'#')) {
        if(raw_contents.find("polyribonucleotide") != std::string::npos) {
           entity_info = raw_contents;  
        } else if(raw_contents.find("ATOM") == std::string::npos || raw_contents.find("HETATM") == std::string::npos) {
            continue;
        } else {
            break;
        }
    }
    
    auto allowed_chains = _get_entity_ids(entity_info);
    
    auto ss = std::stringstream(raw_contents);
    auto atom_contents = Strings{}; 
    auto headers = std::map<String,int>{};    

    auto index(0);
    while(std::getline(ss,line,'\n')) {
        if(line[0] == '_' && line.find(".") != std::string::npos) {
            const auto key = base::trim(base::split_str_by_delimiter(line,".")[1]);
            headers[key] = index++;
             
        } else if (line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos) {
            atom_contents.push_back(line);
        }
    }
    
    const auto x_index = headers.at("Cartn_x");
    const auto y_index = headers.at("Cartn_y");
    const auto z_index = headers.at("Cartn_z");
    const auto atom_name_index = headers.at("auth_atom_id");
    const auto resname_index = headers.at("auth_comp_id");
    const auto chain_id_index = headers.at("label_entity_id"); 
    const auto residue_id_index = headers.at("label_seq_id"); 
    auto i(0);
    
    std::map<String, AtomOPs> residue_atoms;
    
    for(auto& entry : atom_contents) {
        
        const auto tks = base::tokenize_line(entry);
        auto resname = tks[resname_index];
        resname = base::trim(resname);
        auto id = std::stoi(tks[chain_id_index]); 
        if (allowed_chains.find(id) == allowed_chains.end()) { 
            continue;
        } else {
            const auto key =  tks[resname_index] + "|" + tks[residue_id_index]+"|" + tks[chain_id_index];
            if(residue_atoms.find(key) == residue_atoms.end()) {
                residue_atoms[key] = AtomOPs{};
            }
            auto name = tks[atom_name_index];
            residue_atoms[key].push_back(
                std::make_shared<structure::Atom>(
                                base::trim(name),
                                math::Point(
                                       std::stod(tks[x_index]),
                                       std::stod(tks[y_index]),
                                       std::stod(tks[z_index])
                                    )
                                    )
                    );
        }
        
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
CIFParser::_get_new_residue_type(
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
























