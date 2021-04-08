

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/added_motif_library.h"

TEST_CASE( "Test Added Motif Library" ) {
    
    auto path = base::motif_dirs() + "base.motif";
    auto m = motif::file_to_motif(path);
    
    auto mlib = resources::AddedMotifLibrary();
    mlib.add_motif(m);
    
    SUBCASE("try getting back motif that was stored in library") {
        REQUIRE_NOTHROW(mlib.get(m->name()));
        REQUIRE_NOTHROW(mlib.get("", m->end_ids()[0]));
        REQUIRE_NOTHROW(mlib.get("", "", m->ends()[0]->name()));
        REQUIRE_NOTHROW(mlib.get(m->name(), m->end_ids()[0], m->ends()[0]->name()));
        
        REQUIRE_THROWS_AS(mlib.get(), resources::AddedMotifLibraryException);
        REQUIRE_THROWS_AS(mlib.get("TEST"), resources::AddedMotifLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "TEST"), resources::AddedMotifLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "", "TEST"), resources::AddedMotifLibraryException);

    }
    
    SUBCASE("trying to catch adding the same motif twice") {
        REQUIRE_THROWS_AS(mlib.add_motif(m), resources::AddedMotifLibraryException);
    }
    
    SUBCASE("test detecting whether a motif exists in the library") {
        CHECK(mlib.contains(m->name()) == 1);
        CHECK(mlib.contains("", m->end_ids()[0]) == 1);

        CHECK(mlib.contains("TEST") == 0);
        
    }
    
}