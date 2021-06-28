

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/resource_manager.h"

TEST_CASE ("Test Resource Manager ") {

  SUBCASE("test construction of singleton") {
      REQUIRE_NOTHROW(resources::Manager::instance());
  }

  SUBCASE("test storing instance") {
    auto & rm = resources::Manager::instance();
    auto m = rm.motif("HELIX.IDEAL");
  }

  SUBCASE("test individual queries") {
    auto & rm = resources::Manager::instance();
    REQUIRE_NOTHROW(rm.motif("HELIX.IDEAL"));
    REQUIRE_NOTHROW(rm.motif("", "CC_LL_GG_RR"));
    REQUIRE_NOTHROW(rm.motif("", "", "A5-B7"));
    REQUIRE_NOTHROW(rm.motif("HELIX.IDEAL", "CC_LL_GG_RR"));
    REQUIRE_NOTHROW(rm.motif("HELIX.IDEAL", "", "A5-B7"));
    REQUIRE_NOTHROW(rm.motif("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7"));

    REQUIRE_NOTHROW(rm.motif_state("HELIX.IDEAL"));
    REQUIRE_NOTHROW(rm.motif_state("", "CC_LL_GG_RR"));
    REQUIRE_NOTHROW(rm.motif_state("", "", "A5-B7"));
    REQUIRE_NOTHROW(rm.motif_state("HELIX.IDEAL", "CC_LL_GG_RR"));
    REQUIRE_NOTHROW(rm.motif_state("HELIX.IDEAL", "", "A5-B7"));
    REQUIRE_NOTHROW(rm.motif_state("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7"));

    REQUIRE_NOTHROW(rm.motif_state_ensemble("CC_LL_GG_RR"));

    REQUIRE_THROWS_AS(rm.motif("TEST"), resources::ResourceManagerException);
    REQUIRE_THROWS_AS(rm.motif("", "TEST"), resources::ResourceManagerException);
    REQUIRE_THROWS_AS(rm.motif("", "", "TEST"), resources::ResourceManagerException);
  }

  SUBCASE("test adding motif_state_ensembles") {
    auto & rm = resources::Manager::instance();
    auto path = base::unittest_resource_dir() + "resources/test.dat";
    REQUIRE_NOTHROW(resources::Manager::instance().register_extra_motif_ensembles(path));

    // should exist
    CHECK(rm.has_supplied_motif_ensemble("TWOWAY.3TD0.1", "A7-B42") == 1);

    // does not exist wrong end_name
    CHECK(rm.has_supplied_motif_ensemble("TWOWAY.3TD0.1", "A8-B42") == 0);

    auto me = rm.get_supplied_motif_ensemble("TWOWAY.3TD0.1", "A7-B42");
      CHECK(me->size() == 1);
  }

  SUBCASE("test adding motif ensemble") {
    auto & rm = resources::Manager::instance();
    auto csv_path = base::unittest_resource_dir() + "/motif_ensemble/ires.csv";
    auto mes = std::vector<motif::MotifEnsembleOP>();
    motif::motif_ensemble_from_csv_file(csv_path, mes);
    for(auto const & me : mes) {
      auto name = me->members()[0]->motif->name();
      auto end_name = me->members()[0]->motif->end_name(0);
      rm.register_motif_ensemble(name, end_name, me);
      auto new_me = rm.get_supplied_motif_ensemble(name, end_name);
      CHECK(new_me->size() == 1);
      CHECK(new_me->members()[0]->motif->name() == "ires");
      CHECK(new_me->members()[0]->energy == 1);
    }

  }

  SUBCASE("test adding new motifs to resource manager") {

    REQUIRE_NOTHROW(resources::Manager::instance().add_motif(
      base::unittest_resource_dir() + "motif/GAAA_tetraloop"));

    auto m1 = resources::Manager::instance().motif("GAAA_tetraloop");
    CHECK(m1->name() == "GAAA_tetraloop");
  }

}