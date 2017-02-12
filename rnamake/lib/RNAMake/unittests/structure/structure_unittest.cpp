
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "math/numerical.h"
#include "structure/structure.h"

TEST_CASE( "Test Structure", "[Structure]" ) {
    auto rts = ResidueTypeSet();
    auto path = py_unittest_path() + "/resources/p4p6.pdb";
    auto s = structure_from_pdb(path, rts);


    SECTION("Test ability to get specific residues from structure") {
        for(auto const & r : *s) {
        //    std::cout << r->num() << "|" << r->chain_id() << "|" << r->i_code() << "|" << std::endl;
        }

        auto r = s->get_residue(107, 'A', ' ');
        
        REQUIRE(r != nullptr);
        
        auto r2 = s->get_residue(r->uuid());
        
        REQUIRE(r == r2);
        
        SECTION("should get nullptr for residues that dont exist") {
            
            REQUIRE(s->get_residue(1000, 'A', ' ') == nullptr);
            REQUIRE(s->get_residue(106, 'B', ' ') == nullptr);
            
            auto r3 = s->get_residue(107, 'A', ' ');
            auto r3_copy = std::make_shared<Residue>(*r3, 1);
            REQUIRE(s->get_residue(r3_copy->uuid()) == nullptr);
            
        }
    }

    SECTION("Test ability to stringify structure") {
        auto str = s->to_str();
        auto s2 = std::make_shared<Structure>(str, rts);
        REQUIRE(are_structures_equal(s, s2, 0));
    }

    SECTION("Test ability to copy structure") {
        auto s2 = std::make_shared<Structure>(*s);
        REQUIRE(are_structures_equal(s, s2));

        s2 = std::make_shared<Structure>(*s, 1);
        REQUIRE(are_structures_equal(s, s2) == 0);
        REQUIRE(are_structures_equal(s, s2, 0) == 1);

    }

    SECTION("Test moving structure") {
        auto p = Point(10, 0, 0);
        auto s2 = std::make_shared<Structure>(*s);
        s2->move(p);
        
        auto r1 = s->get_residue(0);
        auto r2 = s2->get_residue(0);

        auto dist = r1->center().distance(r2->center());
        REQUIRE(are_floats_equal(dist, 10.0f));
    }

    SECTION("Test applying a tranform to coordinates") {
        auto path = unittest_resource_dir() + "/structure/structure_transformations.dat";
        auto lines = get_lines_from_file(path);

        for(auto i = 0; i < lines.size()-2; i+=3) {

            auto s = structure_from_pdb(resources_path()+"/start/start.pdb", rts);
            auto r = matrix_from_str(lines[i]);
            auto trans = vector_from_str(lines[i+1]);
            auto t = Transform(r, trans);
            auto s2 = std::make_shared<Structure>(lines[i+2], rts);

            s->transform(t);

            for(auto const & r : *s) {
                auto r2 = s2->get_residue(r->num(), r->chain_id(), r->i_code());
                auto dist = r->center().distance(r2->center());
                REQUIRE(dist < 0.001);
            }
        }
    }

    SECTION("Test applying a tranform to coordinates using fast transforms") {
        auto path = unittest_resource_dir() + "/structure/structure_transformations.dat";
        auto lines = get_lines_from_file(path);

        for(auto i = 0; i < lines.size()-2; i+=3) {

            auto s = structure_from_pdb(resources_path()+"/start/start.pdb", rts);
            auto r = matrix_from_str(lines[i]);
            auto trans = vector_from_str(lines[i+1]);
            auto t = Transform(r, trans);
            auto s2 = std::make_shared<Structure>(lines[i+2], rts);

            s->fast_transform(r.transpose(), trans);

            for(auto const & r : *s) {
                auto r2 = s2->get_residue(r->num(), r->chain_id(), r->i_code());
                auto dist = r->center().distance(r2->center());
                REQUIRE(dist < 0.001);
            }
        }
    }

    SECTION("Test getting motif state of structure") {
        auto ss = s->get_state();
        REQUIRE(s->num_residues() == ss->num_residues());

        auto ss_copy = std::make_shared<state::Structure>(*ss);
        REQUIRE(ss_copy->num_residues() == ss->num_residues());

        auto str = ss->to_str();
        ss_copy = std::make_shared<state::Structure>(str);
        REQUIRE(ss_copy->num_residues() == ss->num_residues());

    }
    
}


















