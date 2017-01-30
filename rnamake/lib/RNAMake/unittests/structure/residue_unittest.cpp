


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/file_io.h"
#include "base/settings.h"
#include "math/numerical.h"
#include "structure/residue.h"
#include "structure/is_equal.hpp"


TEST_CASE( "Test Residues for Structure", "[Residue]" ) {
    auto rts = ResidueTypeSet();
    auto path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    auto lines = get_lines_from_file(path);
    auto residues  = ResidueOPs();
    for(auto const & l : lines) {
        if(l.size() < 10) { break; } // end of file
        auto r = std::make_shared<Residue>(l, rts);
        residues.push_back(r);
    }

    SECTION("test getting atoms by name") {
        auto gtype = rts.get_type("GUA");
        auto atoms = AtomOPs();
        atoms.push_back(std::make_shared<Atom>("P", Point(0, 1, 2)));
        auto res = std::make_shared<Residue>(atoms, gtype, "GUA", 1, "A", "");

        // can find a real atom
        auto p_atom = res->get_atom("P");
        REQUIRE(p_atom->name() == "P");
        REQUIRE(p_atom == res->get_atom(0));

        REQUIRE_THROWS_AS(res->get_atom(999), ResidueException);
        REQUIRE_THROWS_AS(res->get_atom("OP1"), ResidueException);
        REQUIRE_THROWS_AS(res->get_atom("fake"), ResidueException);
    }

    SECTION("test has atoms") {
        auto gtype = rts.get_type("GUA");
        auto atoms = AtomOPs();
        atoms.push_back(std::make_shared<Atom>("P", Point(0, 1, 2)));
        auto res = std::make_shared<Residue>(atoms, gtype, "GUA", 1, "A", "");

        REQUIRE(res->has_atom("P") == true);
        REQUIRE(res->has_atom("P1") == false);
        REQUIRE(res->has_atom("OP1") == false);
        REQUIRE(res->has_atom(0) == true);
        REQUIRE(res->has_atom(999) == false);

        res = residues[0];
        // first residue does not have a P atom
        REQUIRE(res->has_atom("P") == false);
        REQUIRE(res->has_atom("C1'") == true);
    }

    SECTION("are residues detecting connections properly") {
        auto r1 = residues[0];
        auto r2 = residues[1];
        auto r3 = residues[2];
        
        SECTION("connecting from 5' to 3'") { REQUIRE(r1->connected_to(*r2) == 1);  }
        SECTION("should not be connected")  { REQUIRE(r1->connected_to(*r3) == 0);  }
        SECTION("connecting from 3' to 5'") { REQUIRE(r3->connected_to(*r2) == -1); }
        
    }

    SECTION("are residues stringifing correctly") {
        auto r = residues[0];
        auto s = r->to_str();
        auto r2 = std::make_shared<Residue>(s, rts);

        SECTION("residues should be the same but not have the same id") {
            REQUIRE(are_residues_equal(r, r2, 0));
        }
    }

    SECTION("are residues copying correctly") {

        auto r = residues[0];
        auto r2 = std::make_shared<Residue>(*r);

        REQUIRE(are_residues_equal(r, r2));

        r2 = std::make_shared<Residue>(*r, 1);

        REQUIRE(are_residues_equal(r, r2) == 0);
        REQUIRE(are_residues_equal(r, r2, 0) == 1);

    }

    SECTION("can move residues properly") {
        auto r = residues[0];
        auto c1 = r->center();
        auto p = Point(1, 0, 0);

        r->move(p);
        auto c2 = r->center();
        auto dist = c1.distance(c2);
        REQUIRE(are_floats_equal(dist, 1.0));
    }

    SECTION("are residues generating steric beads properly") {

        auto r = residues[1];

        REQUIRE(r->num_beads() == 0);
        r->build_beads();
        auto beads = r->beads();

        SECTION("produced the right number of beads") { REQUIRE(beads.size() == 3); }
        
        REQUIRE(beads[0].btype() == BeadType::PHOS);
        REQUIRE(beads[1].btype() == BeadType::SUGAR);
        REQUIRE(beads[2].btype() == BeadType::BASE);

    }

    SECTION("test rotations to be consistent with python") {
        auto path = unittest_resource_dir()+"/residue/residue_transformations.dat";
        auto lines = get_lines_from_file(path);
        for(int i  = 0; i < lines.size()-2; i+=3) {
            auto r = std::make_shared<Residue>(*residues[0]);
            auto rot = matrix_from_str(lines[i]);
            auto trans = vector_from_str(lines[i+1]);
            auto t = Transform(rot, trans);
            auto r_rot = std::make_shared<Residue>(lines[i+2], rts);

            r->transform(t);
            auto dist = r->center().distance(r_rot->center());
            REQUIRE(dist < 0.001);
        }
    }

    SECTION("test fast rotations to be consistent with python") {
        auto path = unittest_resource_dir()+"/residue/residue_transformations.dat";
        auto lines = get_lines_from_file(path);
        for(int i  = 0; i < lines.size()-2; i+=3) {
            auto r = std::make_shared<Residue>(*residues[0]);
            auto rot = matrix_from_str(lines[i]);
            auto trans = vector_from_str(lines[i+1]);
            auto t = Transform(rot, trans);
            auto r_rot = std::make_shared<Residue>(lines[i+2], rts);

            r->fast_transform(rot.transpose(), trans);
            auto dist = r->center().distance(r_rot->center());
            REQUIRE(dist < 0.001);
        }
    }

    SECTION("get residue state") {
        auto r = residues[1];
        r->build_beads();
        auto rs = r->get_state();

        // got all the beads
        REQUIRE(rs->num_beads() == r->num_beads());

        SECTION("are beads properly copied over and not linked") {
            auto p = Point(1, 0, 0);
            rs->move(p);

            auto beads1 = rs->beads();
            auto beads2 = r->beads();

            auto dist = beads1[0].distance(beads2[0]);
            REQUIRE(are_floats_equal(dist, 1.0));
        }

        auto s = rs->to_str();
        auto rs_copy = std::make_shared<state::Residue>(s);

        REQUIRE(rs->num_beads() == rs_copy->num_beads());

        rs_copy = std::make_shared<state::Residue>(*rs);
        REQUIRE(rs->num_beads() == rs_copy->num_beads());
        REQUIRE(*rs == *rs_copy);

    }

}
























