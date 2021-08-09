//
// Created by Joseph Yesselman on 12/15/17.
//

#include <base/log.h>
#include <structure/pose.h>

namespace structure {

PoseOP
get_pose_from_pdb(
        String const & path,
        PDBParser & pdb_parser) {

    LOGI << "generating a pose from pdb file: " + path;

    auto filename = base::filename(path);
    auto name = std::make_shared<base::SimpleString>(filename.substr(0, filename.length() - 4));

    auto residues = pdb_parser.parse(path);
    auto rna_structure = std::make_shared<Structure>(Residues(), Cutpoints());
    auto protein_structure = std::make_shared<Structure>(Residues(), Cutpoints());
    auto small_molecules = std::make_shared<Structure>(Residues(), Cutpoints());

    if(residues->has_protein()) {
        protein_structure = get_structure_from_residues(residues->protein_residues);
    }
    if(residues->has_small_molecules()) {
        small_molecules = get_structure_from_residues(residues->small_molecule_residues);
    }

    auto bps = Basepairs();
    auto end_indexes = Indexes();
    auto end_ids = base::SimpleStringCOPs();
    auto dot_bracket = base::SimpleStringCOP(nullptr);

    if(residues->has_RNA()) {
        rna_structure = get_structure_from_residues(residues->RNA_residues);
        auto x = util::X3dna();
        auto x3dna_basepairs = x.get_basepairs(path);
        bps = get_basepairs_from_x3dna(x3dna_basepairs, *rna_structure)->get_data();
        end_indexes = get_end_indexes_from_basepairs(*rna_structure, bps)->get_data();
        for(auto const & ei : end_indexes) {
            auto end_id = generate_end_id(*rna_structure, bps, bps[ei]);
            end_ids.push_back(std::make_shared<base::SimpleString>(end_id));
        }

        if(end_indexes.size() == 0) {
            dot_bracket = std::make_shared<base::SimpleString>(generate_secondary_structure(*rna_structure, bps));
        }
        else {
            dot_bracket = std::make_shared<base::SimpleString>(
                    primitives::get_dot_bracket_from_end_id(end_ids[0]->get_str()));

        }
    }

    LOGI << filename + " has " + std::to_string(protein_structure->get_num_residues()) +
            " protein residue(s) in " + std::to_string(protein_structure->get_num_chains()) + " chain(s)";
    LOGI << filename + " has " + std::to_string(small_molecules->get_num_residues()) + " small_molecules";
    LOGI << filename + " has " + std::to_string(rna_structure->get_num_residues()) +
            " RNA residue(s) in " + std::to_string(rna_structure->get_num_chains()) + " chain(s)";
    LOGI << filename + " has " + std::to_string(bps.size() + end_indexes.size()) + " basepairs with " +
            std::to_string(end_indexes.size()) + " basepair ends, sites to connect other structure too";

    return std::make_shared<Pose>(*rna_structure, bps, end_indexes, end_ids, name,
                                  *protein_structure, *small_molecules, dot_bracket);
}

}