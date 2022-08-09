#include "../common.hpp"

#include <base/paths.hpp>
#include <structure/all_atom/residue_type.h>
#include <structure/all_atom/residue_type_set.h>

TEST_CASE("Test ResidueType for Structure") {

  using namespace structure::all_atom;

  SUBCASE("Can generate new residue type") {
    auto name = String("GUA");
    auto atom_map = StringIntMap();
    atom_map["P"] = 0;
    auto strings = Strings();
    auto rt = ResidueType(name, atom_map, SetType::RNA, strings);

    CHECK(rt.get_short_name() == 'G');
    CHECK(rt.get_atom_index("P") == 0);
  }
  SUBCASE("test residue type set") {
    ResidueTypeSet rts;
    SUBCASE("Normal residues should be found by 3-letter name") {
      CHECK(rts.contains_residue_type("GUA"));
    }
    SUBCASE("short hand names and alt names also should be found") {
      CHECK(rts.contains_residue_type("A"));
      CHECK(rts.contains_residue_type("rC"));
    }
    CHECK(rts.contains_residue_type("FAKE") == 0);
    CHECK(rts.contains_residue_type("rCC") == 0);
    CHECK(rts.contains_residue_type("AA") == 0);
  }
  SUBCASE("Getting correct residue type by name from set") {
    ResidueTypeSet rts;
    ResidueType rt = rts.get_residue_type("GUA");

    CHECK(rt.get_short_name() == 'G');
    REQUIRE_THROWS_AS(rts.get_residue_type("FAKE"), ResidueTypeException);
  }
  SUBCASE("Getting correct residue type by name from set") {
    ResidueTypeSet rts;
    ResidueType rt = rts.get_residue_type("GUA");

    CHECK(rt.get_short_name() == 'G');
    REQUIRE_THROWS_AS(rts.get_residue_type("FAKE"), ResidueTypeException);
  }
  SUBCASE("can load amino acid residue types") {
    ResidueTypeSet rts;
    auto rt = rts.get_residue_type("ARG");
    CHECK(rt.get_short_name() == 'A');
    CHECK(rt.get_atom_index("N") == 0);
  }
}
