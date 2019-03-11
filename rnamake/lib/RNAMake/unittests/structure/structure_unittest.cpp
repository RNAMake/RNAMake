
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/structure.h"
#include "structure/is_equal.hpp"


TEST_CASE( "Test Structure", "[Structure]" ) {
    auto path = base::unittest_resource_dir() + "/structure/test_str_to_structure.dat";
    auto lines = base::get_lines_from_file(path);
    auto rts = ResidueTypeSet();
    auto s = std::make_shared<Structure>(lines[0], rts);
    
    SECTION("Test ability to get specific residues from structure") {
        auto r = s->get_residue(107, "A", "");
        
        REQUIRE(r != nullptr);
        
        auto r2 = s->get_residue(r->uuid());
        
        REQUIRE(r == r2);
        
        SECTION("should get nullptr for residues that dont exist") {
            
            REQUIRE(s->get_residue(1000, "A", "") == nullptr);
            REQUIRE(s->get_residue(106, "B", "") == nullptr);
            
            auto r3 = s->get_residue(107, "A", "");
            auto r3_copy = std::make_shared<Residue>(*r3);
            r3_copy->new_uuid();
            
            REQUIRE(s->get_residue(r3_copy->uuid()) == nullptr);
            
        }
    }
    
    SECTION("Should contain the correct number of residues") {
        REQUIRE(s->residues().size() == 157);
    }
    
    SECTION("Should contain the correct number of atoms") {
        REQUIRE(s->atoms().size() == 3357);
    }
    
    SECTION("Test ability to stringify structure") {
        auto str = s->to_str();
        auto s2 = std::make_shared<Structure>(str, rts);
        
        REQUIRE(s->residues().size() == s2->residues().size());
        REQUIRE(are_structures_equal(s, s2, 0));
    }
    
    SECTION("Test ability to copy structure") {
        auto s2 = std::make_shared<Structure>(*s);
        
        REQUIRE(are_structures_equal(s, s2));
    }
    
    SECTION("Collecting beads for all residues") {
        REQUIRE(s->get_beads().size() == 470);
        
        SECTION("excluding a residue should remove 3 beads") {
            auto r = s->get_residue(106, "A", "");
            REQUIRE(s->get_beads(ResidueOPs{r}).size() == 467);
        }
            
    }
    
    SECTION("Test moving structure") {
        auto p = math::Point(10, 0, 0);
        auto s2 = std::make_shared<Structure>(*s);
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
        
        REQUIRE(flag == 0);
        
        
    }
    
    SECTION("Test applying a tranform to coordinates") {
        auto path = base::unittest_resource_dir() + "/structure/test_transform.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto r = math::matrix_from_str(lines[0]);
        auto trans = math::vector_from_str(lines[1]);
        auto t = math::Transform(r, trans);
        auto s2 = std::make_shared<Structure>(*s);
        
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
        
        REQUIRE(flag == 0);
        
    }
    
}


















