//
//  residue_type_set.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <dirent.h>
#include <fstream>

//RNAMake Headers
#include <structure/residue_type_set.h>
#include <base/paths.h>

namespace structure {

    ResidueTypeSet::ResidueTypeSet():
            residue_types_(ResidueTypeOPs()) {
        auto base_path = base::resources_path() + "/residue_types/";

        _read_rtypes_from_dir(base_path + "RNA/", SetType::RNA);
        _read_rtypes_from_dir(base_path + "PROTEIN/", SetType::PROTEIN);
    }


    void
    ResidueTypeSet::_read_rtypes_from_dir(
            String const & path,
            SetType const & set_type) {

        DIR *pDIR;
        struct dirent *entry;
        pDIR = opendir(path.c_str());
        while ((entry = readdir(pDIR)) != NULL) {
            auto fname = String(entry->d_name);
            if (fname.length() < 4) { continue; }

            auto name = _get_rtype_name(path + fname);
            auto atom_map = _get_atom_map_from_file(path + fname);
            auto alt_names = _get_extra_resnames_for_specific_res(name);
            auto rtype = std::make_shared<ResidueType>(name, atom_map, set_type, alt_names);
            residue_types_.push_back(rtype);
        }
        closedir(pDIR);
        delete entry;
    }

    String
    ResidueTypeSet::_get_rtype_name(
            String const & fname) {
        auto name_spl = base::split_str_by_delimiter(fname, "/");
        auto type_file_name = name_spl.back();
        auto spl_2 = base::split_str_by_delimiter(type_file_name, ".");
        return spl_2[0];

    }

    StringIntMap
    ResidueTypeSet::_get_atom_map_from_file(
            String const & fname) {

        auto line = String();
        std::ifstream input;
        input.open(fname);
        getline(input, line);
        input.close();
        auto atom_names = base::split_str_by_delimiter(line, " ");
        auto i = 0;
        auto atom_map = StringIntMap();
        for (auto const & name : atom_names) {
            atom_map[name] = i;
            i++;
        }
        return atom_map;
    }

    ResidueTypeCOP
    ResidueTypeSet::get_residue_type(
            String const & resname) const {

        for (auto restype : residue_types_) {
            if (restype->is_valid_residue_name(resname)) { return restype; }
        }

        throw ResidueTypeException("cannot find residue with name :" + resname);
    }


    bool
    ResidueTypeSet::contains_residue_type(
            String const & resname) const {

        for (auto const & restype : residue_types_) {
            if (restype->is_valid_residue_name(resname)) { return true; }
        }
        return false;
    }

    Strings
    ResidueTypeSet::_get_extra_resnames_for_specific_res(
            String const & res_name) {
        auto alt_names = Strings();
        if        (res_name == "GUA") {
            alt_names = base::split_str_by_delimiter("MIA GDP GTP M2G 1MG 7MG G7M QUO I YG", " ");
        } else if (res_name == "ADE") {
            alt_names = base::split_str_by_delimiter("A23 3DA 1MA 12A AET 2MA", " ");
        } else if (res_name == "URA") {
            alt_names = base::split_str_by_delimiter("PSU H2U 5MU 4SU 5BU 5MC U3H 2MU 70U BRU DT", " ");
        } else if (res_name == "CYT") {
            alt_names = base::split_str_by_delimiter("CBR CCC", " ");
        }
        return alt_names;
    }

}