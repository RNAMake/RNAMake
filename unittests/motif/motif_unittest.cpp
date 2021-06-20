

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "util/basic_io.hpp"
#include "util/steric_lookup.hpp"
#include "structure/is_equal.h"
#include "motif/motif.h"
#include "motif/motif_factory.h"

TEST_CASE( "Test Motifs the core of everything!" ) {
    
    auto path = base::motif_dirs() + "base.motif";
    auto m = motif::file_to_motif(path);
    path = base::motif_dirs() + "ref.motif";
    auto ref_m = motif::file_to_motif(path);
    auto rts = structure::ResidueTypeSet();
    
    CHECK(m->ends().size() == 2);
    CHECK(m->basepairs().size() == 5);
    
    SUBCASE("test copy constuctor for motif") {
        auto m_copy = motif::Motif(*m);
        
        CHECK(m_copy.ends().size() == 2);
        CHECK(m_copy.basepairs().size() == 5);

        for(auto const & r : m->residues()) {
            CHECK(m_copy.get_residue(r->uuid()) != nullptr);
        }
        
        for(auto const & bp : m->basepairs()) {
            CHECK(m_copy.get_basepair(bp->uuid()).size() == 1);
        }
        
        
        CHECK(m_copy.atoms().size() == m->atoms().size());
        auto atoms1 = m->atoms();
        auto atoms2 = m_copy.atoms();
        
        CHECK(are_atom_vectors_equal(atoms1, atoms2));

    }
    
    SUBCASE("test stringifying motif") {
        auto s = m->to_str();
        auto m_copy = motif::Motif(s, rts);
        
        CHECK(m_copy.ends().size() == 2);
        CHECK(m_copy.basepairs().size() == 5);

        for(auto const & r : m->residues()) {
            CHECK(m_copy.get_residue(r->num(), r->chain_id(), r->i_code()) != nullptr);
        }
        
        CHECK(m_copy.atoms().size() == m->atoms().size());
        auto atoms1 = m->atoms();
        auto atoms2 = m_copy.atoms();
        
        CHECK(are_atom_vectors_equal(atoms1, atoms2));
    }

    SUBCASE("test aligning motifs") {
        m->move(math::Point(10, 10, 10));
        auto m_aligned = get_aligned_motif(ref_m->ends()[0], m->ends()[0], m);
        
        auto dist = ref_m->ends()[0]->d().distance(m_aligned->ends()[0]->d());
        auto r_dist = ref_m->ends()[0]->r().difference(m_aligned->ends()[0]->r());

        CHECK(dist < 1.0);
        CHECK(r_dist < 0.0001);
        
    }
    
    SUBCASE("Aligning motifs directly from pdbs") {
        auto mf = motif::MotifFactory{};
        auto motifs = mf.motif_from_file(base::unittest_resource_dir() + "/pdbs/124D.pdb");
    }

    SUBCASE("test that repeat aligning does not cause error") {
        m->move(math::Point(10, 10, 10));
        auto ref_bp = ref_m->ends()[0]->state();
        
        auto dist = 0.0f, r_dist = 0.0f;
        
        for(int i = 0; i < 100; i++) {
            align_motif(ref_bp, m->ends()[0], m);
            
            dist = ref_bp->d().distance(m->ends()[0]->d());
            r_dist = ref_bp->r().difference(m->ends()[0]->r());
            
        }
        
        CHECK(dist < 1.0);
        CHECK(r_dist < 0.0001);

    }
    
    SUBCASE("test getting new unique indentifers") {
        auto m_copy = motif::Motif(*m);
        m_copy.new_res_uuids();
        
        SUBCASE("all residues and basepairs have new indentifers cant search with originals anymore") {
        
            for(auto const & r : m->residues()) {
                CHECK(m_copy.get_residue(r->uuid()) == nullptr);
            }
        
            for(auto const & bp : m->basepairs()) {
                CHECK(m_copy.get_basepair(bp->uuid()).size() == 0);
            }
        }
        
        SUBCASE("should still be able to find them using names to find them") {
            for(auto const & r : m->residues()) {
                CHECK(m_copy.get_residue(r->num(), r->chain_id(), r->i_code()) != nullptr);
            }
        }

        
    }
    
    SUBCASE("test copying uuids from one motif to another") {
        auto path = base::motif_dirs() + "base.motif";
        auto m1 = motif::file_to_motif(path);
        auto m2 = motif::file_to_motif(path);

        CHECK(m1->id() != m2->id());
        m1->copy_uuids_from_motif(*m2);
        
        CHECK(m1->id() == m2->id());
        
        int i = 0;
        for(auto const & bp : m1->basepairs()) {
            CHECK(bp->uuid() == m2->basepairs()[i]->uuid());
            i++;
        }
        
        i = 0;
        for(auto const & r : m1->residues()) {
            CHECK(r->uuid() == m2->residues()[i]->uuid());
            i++;
        }
        
    }

    SUBCASE("test steric look based on beads and atoms") {
        auto path = base::motif_dirs() + "ref.motif";
        auto m1 = motif::file_to_motif(path);
        auto lookup = util::StericLookupNew();
        auto points = math::Points();
        for(auto const & b : m1->get_beads()) {
            points.push_back(b.center());
        }
        lookup.add_points(points);
    }
}







