//
// Created by Joseph Yesselman on 12/3/17.
//

#include <iostream>
#include "../common.hpp"

#include <structure/residue_type_set.h>
#include <structure/segment_factory.h>
#include <structure/segment.h>

TEST_CASE( "Test all atom segment") {

    //TODO segment unittest uses segment factory maybe reweork it into two files?

    SUBCASE("test loading segments from pdbs") {
        auto rts = structure::ResidueTypeSet();
        auto seg_factory = structure::SegmentFactory(rts);
        auto path = base::unittest_resource_dir() + "all_atom/HELIX.IDEAL.2/HELIX.IDEAL.2.pdb";
        auto seg = seg_factory.segment_from_pdb(path, util::SegmentType::HELIX, false);

        auto segs = seg_factory.all_segments_from_pdb(path, util::SegmentType::HELIX, false);

        int i = -1;
        for (auto & seg: segs) {
            i++;
            seg_factory.align_segment_to_ref_frame(*seg);
            seg->write_pdb("test." + std::to_string(i) + ".pdb");
        }



        //TODO change the get_json to to a string function
//        auto j = seg->get_json();
//        auto seg2 = structure::Segment(j, rts);
//        CHECK(seg2.is_equal(*seg, false));

//        auto p = math::Point(2, 2, 2);
//        seg2.move(p);
//        CHECK(!seg2.is_equal(*seg, false));
    }

    SUBCASE("test generating segments from components") {
        auto rts = structure::ResidueTypeSet();
        auto seg_factory = structure::SegmentFactory(rts);
        auto path = base::unittest_resource_dir() + "/all_atom/HELIX.IDEAL.2/HELIX.IDEAL.2.pdb";

        auto parser = structure::PDBParser(rts);
        auto residues = parser.parse(path);

        auto rna_structure = structure::get_structure_from_residues(residues->RNA_residues);

        auto x = util::X3dna();
        x.set_rebuild_files(false);
        auto x3dna_bps = x.get_basepairs(path);
        auto bps = get_basepairs_from_x3dna(x3dna_bps, *rna_structure)->get_data();

        auto proteins = structure::Structure(structure::Residues(), Cutpoints());
        auto small_molecules = structure::Structure(structure::Residues(), Cutpoints());

        auto seg1 = seg_factory.segment_from_components("HELIX.IDEAL.2", *rna_structure, bps, proteins,
                                                        small_molecules, util::SegmentType::HELIX);

        auto seg2 = seg_factory.segment_from_pdb(path, util::SegmentType::HELIX, false);

        CHECK(seg1->is_equal(*seg2, false));


    }



}
