//
//  rna_structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

//RNAMake Headers
#include "structure/rna_structure.h"

RNAStructure::RNAStructure(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends,
        SimpleStringOPs const & end_ids,
        SimpleStringOP const & name,
        int block_end_add,
        SimpleStringOP const & dot_bracket) :
        primitives::RNAStructure<Basepair, Structure, Chain, Residue>(
                structure, basepairs, ends, end_ids, name),
        block_end_add_(block_end_add),
        dot_bracket_(dot_bracket),
        protein_beads_(Beads()) {}

RNAStructure::RNAStructure(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends,
        SimpleStringOPs const & end_ids,
        SimpleStringOP const & name,
        int block_end_add,
        SimpleStringOP const & dot_bracket,
        Beads const & protein_beads ):
        primitives::RNAStructure<Basepair, Structure, Chain, Residue>(
                structure, basepairs, ends, end_ids, name),
        block_end_add_(block_end_add),
        dot_bracket_(dot_bracket),
        protein_beads_(protein_beads) {}

RNAStructure::RNAStructure(
        RNAStructure const & rs,
        int new_uuid):
        primitives::RNAStructure<Basepair, Structure, Chain, Residue>(),
        block_end_add_(rs.block_end_add_),
        dot_bracket_(rs.dot_bracket_),
        protein_beads_(rs.protein_beads_) {

    structure_ = std::make_shared<Structure>(*rs.structure_, new_uuid);
    end_ids_   = rs.end_ids_;
    name_      = rs.name_;

    for(auto const & bp : rs.basepairs_) {
        if(new_uuid) {
            auto bp_res = rs.get_bp_res(*bp);
            auto res1 = structure_->get_residue(bp_res[0]->num(), bp_res[0]->chain_id(), bp_res[0]->i_code());
            auto res2 = structure_->get_residue(bp_res[1]->num(), bp_res[1]->chain_id(), bp_res[1]->i_code());
            auto bp_new = std::make_shared<Basepair>(*bp, res1->uuid(), res2->uuid(), Uuid());
            basepairs_.push_back(bp_new);

        }
        else {
            basepairs_.push_back(std::make_shared<Basepair>(*bp));
        }
    }

    for(auto const & end : rs.ends_) {
        if(new_uuid) {
            auto bp_res = rs.get_bp_res(*end);
            auto res1 = structure_->get_residue(bp_res[0]->num(), bp_res[0]->chain_id(), bp_res[0]->i_code());
            auto res2 = structure_->get_residue(bp_res[1]->num(), bp_res[1]->chain_id(), bp_res[1]->i_code());
            auto end_new = std::make_shared<Basepair>(*end, res1->uuid(), res2->uuid(), Uuid());
            ends_.push_back(end_new);

        }
        else {
            ends_.push_back(std::make_shared<Basepair>(*end));
        }
    }

}

RNAStructure::RNAStructure(
        String const & s,
        ResidueTypeSet const & rts):
        primitives::RNAStructure<Basepair, Structure, Chain, Residue>() {

    auto spl = split_str_by_delimiter(s, "&");
    name_ = std::make_shared<SimpleString>(spl[0]);
    block_end_add_ = std::stoi(spl[1]);
    structure_ = std::make_shared<Structure>(spl[2], rts);
    auto bp_strs = split_str_by_delimiter(spl[3], "@");
    for(int i = 0; i < bp_strs.size()-1; i++) {
        //basepairs_.push_back()
    }
}

String
RNAStructure::to_str() {
    auto s = name_->to_str() + "&" + std::to_string(block_end_add_) + "&" + structure_->to_str() + "&";
    for(auto const & bp : basepairs_) {
        auto bp_res = get_bp_res(*bp);
        s += bp->to_str() + ";";
        s += std::to_string(bp_res[0]->num()) + "|" + bp_res[0]->chain_id() + "|" + bp_res[0]->i_code() + ";";
        s += std::to_string(bp_res[1]->num()) + "|" + bp_res[1]->chain_id() + "|" + bp_res[1]->i_code() + "@";
    }
    s += "&";
    for(auto const & end : ends_) {
        auto bp_res = get_bp_res(*end);
        s += end->to_str() + ";";
        s += std::to_string(bp_res[0]->num()) + "|" + bp_res[0]->chain_id() + "|" + bp_res[0]->i_code() + ";";
        s += std::to_string(bp_res[1]->num()) + "|" + bp_res[1]->chain_id() + "|" + bp_res[1]->i_code() + "@";
    }
    s += "&";
    for(auto const & end_id : end_ids_) { s += end_id->to_str() + " "; }
    s += "&";
    for(auto const & b : protein_beads_) { s += b.to_str() + ";"; }
    s += "&" + dot_bracket_->to_str() + "&";
    return s;
}

/*


// get end infomation functions ////////////////////////////////////////////////////////////////////


int
RNAStructure::get_end_index(BasepairOP const & end) {
    if(std::find(ends_.begin(), ends_.end(), end) == ends_.end()) {
        throw RNAStructureException("cannot find end: " + end->name() + " in ends ");
    }
    
    int pos = (int)(std::find(ends_.begin(), ends_.end(), end) - ends_.begin());
    return pos;
}

int
RNAStructure::get_end_index(String const & end_name) {
    for(int i = 0; i < ends_.size(); i++) {
        if(ends_[i]->name() == end_name) { return i; }
    }
    
    throw RNAStructureException("could not find a end basepair with name: " + end_name);
    
    return -1;
}


// output functions ////////////////////////////////////////////////////////////////////////////////


String const
RNAStructure::to_pdb_str(
    int renumber) {
    return structure_->to_pdb_str(renumber);
}

void
RNAStructure::to_pdb(
    String const fname,
    int renumber) {
    return structure_->to_pdb(fname, renumber);
}


// non member functions ////////////////////////////////////////////////////////////////////////////


std::shared_ptr<BasepairOPs>
end_from_basepairs(
    StructureOP const & s,
    BasepairOPs const & bps) {
    
    auto chain_ends = ResidueOPs();
    for(auto const & c : s->chains()) {
        chain_ends.push_back(c->first());
        if(c->length() > 1) {
            chain_ends.push_back(c->last());
        }
    }
    
    auto res_map = std::map<Uuid, ResidueOP, UuidCompare>();
    for(auto const & r : chain_ends ) { res_map[r->uuid()] = r;}

    auto ends = std::make_shared<BasepairOPs>();
    for(auto const & bp : bps) {
        if(bp->bp_type() != "cW-W") { continue; }
        if(res_map.find(bp->res1()->uuid()) != res_map.end() &&
           res_map.find(bp->res2()->uuid()) != res_map.end() ) {
            ends->push_back(bp);
        }
    }
    
    return ends;
    
}

std::shared_ptr<BasepairOPs>
subselect_basepairs_with_res(
    ResidueOPs const & res,
    BasepairOPs const & all_bps) {
    
    auto res_map = std::map<Uuid, ResidueOP, UuidCompare>();
    auto bps = std::make_shared<BasepairOPs>();
    
    for(auto const & r : res ) { res_map[r->uuid()] = r;}

    for(auto const & bp : all_bps) {
        if(res_map.find(bp->res1()->uuid()) != res_map.end() &&
           res_map.find(bp->res2()->uuid()) != res_map.end() )  {
            bps->push_back(bp);
        }
    }

    return bps;
}

 */

RNAStructureOP
rna_structure_from_pdb(
        String const & path,
        ResidueTypeSet const & rts) {
    auto s = structure_from_pdb(path, rts);
    auto bps = basepairs_from_x3dna(path, *s);
    auto ends = ends_from_basepairs(*s, bps);
    auto end_ids = SimpleStringOPs();
    auto name = std::make_shared<SimpleString>("");

    for(auto const & end : ends) {
        bps.erase(std::remove(bps.begin(), bps.end(), end), bps.end());
    }

    auto first_ss_id = String("");
    auto i = 0;
    for(auto const & end : ends) {
        auto ss_id = assign_end_id(*s, bps, ends, end);
        end_ids.push_back(std::make_shared<SimpleString>(ss_id));
        if(i == 0) { first_ss_id = ss_id; }
        i++;
    }

    auto seq_and_ss = primitives::end_id_to_seq_and_db(first_ss_id);
    auto dot_bracket = std::make_shared<SimpleString>(seq_and_ss[1]);
    return std::make_shared<RNAStructure>(s, bps, ends, end_ids, name, -1, dot_bracket);
}

BasepairOPs
basepairs_from_x3dna(
        String const & path,
        Structure const & s) {

    auto basepairs = BasepairOPs();
    auto x3dna_parser = X3dna();
    auto x_basepairs = x3dna_parser.get_basepairs(path);
    ResidueOP res1, res2;
    BasepairOP bp;
    for(auto const & xbp : x_basepairs) {
        res1 = s.get_residue(xbp.res1.num, xbp.res1.chain_id, xbp.res1.i_code);
        res2 = s.get_residue(xbp.res2.num, xbp.res2.chain_id, xbp.res2.i_code);
        if (res1 == nullptr || res2 == nullptr) { continue; }
        try {
            res1->get_atom("C1'");
            res2->get_atom("C1'");
        }
        catch(ResidueException) { continue; }

        auto res = ResidueOPs{res1, res2};
        auto bp_name_str = primitives::calc_bp_name(res);
        auto bp_type = primitives::get_bp_type<Residue>(res, xbp.bp_type);
        auto bp_name = std::make_shared<SimpleString>(bp_name_str);
        auto center = _calc_center(res);
        bp = std::make_shared<Basepair>(res1->uuid(), res2->uuid(), xbp.r, center,
                                        res1->get_atom("C1'")->coords(), res2->get_atom("C1'")->coords(),
                                        bp_name, xbp.bp_type, bp_type, Uuid());
        basepairs.push_back(bp);
    }
    return basepairs;
}

bool
are_rna_strucs_equal(
        RNAStructure const & rs1,
        RNAStructure const & rs2,
        bool check_uuid) {

    if(rs1.num_residues() != rs2.num_residues()) { return false; }
    bool result;
    for(int i = 0; i < rs1.num_residues(); i++) {
        result = are_residues_equal(rs1.get_residue(i), rs2.get_residue(i), check_uuid);
        if(!result) { return false; }
    }

    if(rs1.num_basepairs() != rs2.num_basepairs()) { return false; }
    for(int i = 0; i < rs1.num_basepairs(); i++) {
        result = are_basepairs_equal(rs1.get_basepair(i), rs2.get_basepair(i), check_uuid);
        if(!result) { return false; }
    }

    if(rs1.num_ends() != rs2.num_ends()) { return false; }
    for(int i = 0; i < rs1.num_ends(); i++) {
        result = are_basepairs_equal(rs1.get_end(i), rs2.get_end(i), check_uuid);
        if(!result) { return false; }
    }
    return true;

}








