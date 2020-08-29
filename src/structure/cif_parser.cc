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
#include <regex>

//RNAMake Headers
#include "base/string.h"
#include "base/file_io.h"
#include "math/xyz_vector.h"
#include "structure/cif_parser.h"

// td  Residue(
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
    auto tks = base::tokenize_line(entity_contents);
    const auto num_entries = tks.size()/headers.size() ;
    // check that the sizing is correct
    if(num_entries*headers.size() != tks.size()) {
        const auto msg = String{"header contents are mis-sized, mal-formed or the parser does not work"};
        throw CIFParseError{msg};
    }
    // find the indicies we need
    const auto res_type_index = headers.at("pdbx_description");
    const auto id_index = headers.at("id");
    // loop through and add the indices if they are RNA, add the index
    auto allowed_ids = std::unordered_set<int>{};
    for(auto index = 0; index < num_entries; ++ index) {
        const auto& chain_name = tks[index*headers.size()+res_type_index] ;
        if( chain_name.find("polyribonucleotide") != std::string::npos ||
            chain_name.find("RNA") != std::string::npos    ) {
            allowed_ids.insert( 
                            std::stoi(tks[index*headers.size() +id_index])
                                );
        }
    }
    return allowed_ids;
}


std::map<String, AtomOPs>
CIFParser::_get_atoms(
            String const & atom_info,
            std::unordered_set<int> const & allowed_chains
                        ) {
    // helper method that gets the atoms out of an input area 
    auto ss = std::stringstream(atom_info);
    auto atom_contents = Strings{}; 
    auto headers = std::map<String,int>{};    
    auto line = String{};
    
    auto index(0);
    while(std::getline(ss,line,'\n')) {
        if(line.find("_atom_site") != std::string::npos) {
            const auto key = base::trim(base::split_str_by_delimiter(line,".")[1]);
            headers[key] = index++;
             
        } else if (line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos) {
            atom_contents.push_back(line);
        }
    }
    // pretty dumb looking but need to identify each of the indices needed up front 
    const auto x_index = headers.at("Cartn_x");
    const auto y_index = headers.at("Cartn_y");
    const auto z_index = headers.at("Cartn_z");
    const auto atom_name_index = headers.at("auth_atom_id");
    const auto resname_index = headers.at("auth_comp_id");
    const auto chain_id_numeric_index = headers.at("label_entity_id"); 
    const auto chain_id_index = headers.at("label_asym_id"); 
    const auto residue_id_index = headers.at("label_seq_id"); 
    const auto ins_code_index = headers.at("pdbx_PDB_ins_code"); 
     
    auto residue_atoms = std::map<String, AtomOPs>{} ;
    for(auto& entry : atom_contents) {
        
        const auto tks = base::tokenize_line(entry);
        auto resname = tks[resname_index];
        resname = base::trim(resname);
        auto id = std::stoi(tks[chain_id_numeric_index]); 
        
        if (allowed_chains.find(id) == allowed_chains.end()) { 
            continue;
        } else {
            auto insertion_code = tks[ins_code_index];
             
            if(insertion_code == "." || insertion_code == "?") {
                insertion_code = "NA"; 
            }
            const auto key =  tks[resname_index] + "|" + tks[residue_id_index]+"|" 
                    + tks[chain_id_index] + "|" + insertion_code ;
            if(residue_atoms.find(key) == residue_atoms.end()) {
                residue_atoms[key] = AtomOPs{};
            }
            
            auto name = tks[atom_name_index];
            std::replace(name.begin(),name.end(),'\"',' '); 
            name = base::trim(name);
            residue_atoms[key].emplace_back(
                std::make_shared<structure::Atom>(
                                name,
                                math::Point(std::stod(tks[x_index]),std::stod(tks[y_index]),std::stod(tks[z_index])
                                    )
                                    )
                    );
        }
        
    }

    return residue_atoms;
}



ResidueOPs const &
CIFParser::parse(
        String const & cif_file,
        int protein,
        int rna,
        int others) {
    // This is the driver method for the class. It performs the following actions:
    // 1. Clear residues (if not empty)
    // 2. Read in selected raw contents from the file in question 
    //      - throw if file is ill-formed
    // 3. Used "_entity" information to determine which atoms should be used
    // 4. Organize atom objects for the atoms we would like to keep
    // 5. Contstruct Residues from the atom groupings... TODO this part of the method needs some work.
    if(!base::file_exists(cif_file)) {
        throw std::runtime_error("doesn't exist");
    }
    // in case you use the same parser for multiple files  
    if(!residues_.empty()) {
        residues_.clear();
    }

    auto infile = std::ifstream(cif_file);
    auto line = String{}; 
    auto raw_contents = String{}; 
   
    if(!infile.is_open()) {
        const auto msg = String{"Error opening file: " + cif_file + ". Are you sure it exists?"};
        throw CIFParseError(msg);
    }

    auto entity_info = String{}; 
    auto atom_info = String{};
    auto contents = std::map<String,String>{
                                    {"_entity" , ""},
                                    {"_atom_site" , ""},
                                    //{"_pdbx_poly_seq_scheme" , ""},
                                    //{"_ndb_struct_na_base_pair" , ""}
                                 };

    auto var_start = std::regex("\n_[a-zA-Z_]*");
    auto matches = std::smatch{};
    auto section = String{};
    
    while(std::getline(infile,section,'#')) {
        if(section.find("loop_\n") == std::string::npos){
            continue;
        }
        std::regex_search(section,matches,var_start); 
        auto key = String{matches[0]};
        key = base::trim(key);
        
        if(contents.find(key) != contents.end()) {
            contents[key] = section;
        }
    }
    // always close!
    infile.close(); 
    
    // checking to make sure that all the required content is present in the file
    for(auto& kv : contents) {
        if(kv.second.empty()) { 
            auto msg = String{"Required section " + kv.first + " not found in CIF file " + cif_file + ". This file is likely ill formed and/or a non-ASCII file."};
            throw CIFParseError(msg);
        }
    }
    
    // deal with entity information
    auto allowed_ids = _get_entity_ids(contents["_entity"]);
    // group the atoms 
    auto residue_atoms = _get_atoms(contents["_atom_site"],allowed_ids);
    
    for (auto & kv : residue_atoms) {
        const auto res_tokens = base::split_str_by_delimiter(kv.first,"|"); 
        ResidueType rtype;
        if (!rts_.contains_rtype(res_tokens[0]) && !others) {
            continue; 
        }
        if (!others) {
            rtype = rts_.get_rtype_by_resname(res_tokens[0]);
        } else {
            auto atom_names = Strings();
            for (auto const & a : kv.second) { atom_names.push_back(a->name()); }
            rtype = _get_new_residue_type(res_tokens[0], atom_names);
        }

        if (!protein && rtype.set_type() == SetType::PROTEIN && res_tokens[0].size() > 1) { continue; }
        if (!rna && rtype.set_type() == SetType::RNA) { continue; }


        residues_.emplace_back(
                std::make_shared<Residue>(rtype, res_tokens[0], std::stoi(res_tokens[1]), res_tokens[2], res_tokens[3])
                );
        (*residues_.rbegin())->setup_atoms(kv.second); // that's uglyyyy, but it works :). TODO take atoms at ctor

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
























