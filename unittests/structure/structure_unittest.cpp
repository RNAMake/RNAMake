

#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/structure.h"
#include "structure/is_equal.h"


TEST_CASE( "Test Structure" ) {
    auto path = base::unittest_resource_dir() + "/structure/test_str_to_structure.dat";
    auto lines = base::get_lines_from_file(path);
    auto rts = structure::ResidueTypeSet();
    auto s = std::make_shared<structure::Structure>(lines[0], rts);
    
    SUBCASE("Test ability to get specific residues from structure") {
        auto r = s->get_residue(107, "A", "");
        
        CHECK(r != nullptr);
        
        auto r2 = s->get_residue(r->uuid());
        
        CHECK(r == r2);
        
        SUBCASE("should get nullptr for residues that dont exist") {
            
            CHECK(s->get_residue(1000, "A", "") == nullptr);
            CHECK(s->get_residue(106, "B", "") == nullptr);
            
            auto r3 = s->get_residue(107, "A", "");
            auto r3_copy = std::make_shared<structure::Residue>(*r3);
            r3_copy->new_uuid();
            
            CHECK(s->get_residue(r3_copy->uuid()) == nullptr);
            
        }
    }
    
    SUBCASE("Should contain the correct number of residues") {
        CHECK(s->residues().size() == 157);
    }
    
    SUBCASE("Should contain the correct number of atoms") {
        CHECK(s->atoms().size() == 3357);
    }
    
    SUBCASE("Test ability to stringify structure") {
        auto str = s->to_str();
        auto s2 = std::make_shared<structure::Structure>(str, rts);
        
        CHECK(s->residues().size() == s2->residues().size());
        CHECK(are_structures_equal(s, s2, 0));
    }
    
    SUBCASE("Test ability to copy structure") {
        auto s2 = std::make_shared<structure::Structure>(*s);
        
        CHECK(are_structures_equal(s, s2));
    }
    
    SUBCASE("Collecting beads for all residues") {
        CHECK(s->get_beads().size() == 470);
        
        SUBCASE("excluding a residue should remove 3 beads") {
            auto r = s->get_residue(106, "A", "");
            CHECK(s->get_beads(structure::ResidueOPs{r}).size() == 467);
        }
            
    }
    
    SUBCASE("Test moving structure") {
        auto p = math::Point(10, 0, 0);
        auto s2 = std::make_shared<structure::Structure>(*s);
        s2->move(p);
        
        auto new_atoms = s2->atoms();
        auto org_atoms = s->atoms();
        
        float dist;
        int flag = 0;
        for(int i = 0; i < org_atoms.size(); i++) {
            if(org_atoms[i].get() == NULL) { continue; }
            dist = 10 - org_atoms[i]->coords().distance(new_atoms[i]->coords());
            if(dist > 0.001) { flag = 1; }
        }
        
        CHECK(flag == 0);
        
        
    }
    
    SUBCASE("Test applying a tranform to coordinates") {
        auto path = base::unittest_resource_dir() + "/structure/test_transform.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto r = math::matrix_from_str(lines[0]);
        auto trans = math::vector_from_str(lines[1]);
        auto t = math::Transform(r, trans);
        auto s2 = std::make_shared<structure::Structure>(*s);
        
        s2->transform(t);
        
        auto new_atoms = s2->atoms();
        auto org_atoms = s->atoms();
        
        float dist;
        int flag = 0;
        for(int i = 0; i < org_atoms.size(); i++) {
            if(org_atoms[i].get() == NULL) { continue; }
            dist = 10 - org_atoms[i]->coords().distance(new_atoms[i]->coords());
            if(dist > 0.001) { flag = 1; }
        }
        
        CHECK(flag == 0);
        
    }
    
}


















