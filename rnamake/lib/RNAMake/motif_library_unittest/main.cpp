//
//  main.cpp
//  motif_library_unittest
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "sqlite3_connection.h"
#include "motif_library.h"
#include "motif_type.h"
#include "motif.h"
#include "residue_type_set.h"

int
test_sqlite3_connection() {
    Sqlite3Connection conn(HELIX);
    conn.query("SELECT * FROM motifs");
    Strings values = conn.next();
    ResidueTypeSet rts;
    Motif m ( values[0], rts );
    m.to_pdb("test.pdb");
    // HELIX.4KQ0.0
    return 1;
}

int
test_load_all() {
    MotifLibrary mlib(HELIX);
    mlib.load_all(10);
    return 1;
}

int
test_get_motif() {
    MotifLibrary mlib(HELIX);
    Motif m = mlib.get_motif("HELIX.IDEAL");
    m.to_pdb("test2.pdb");
    return 1;
}

int
test_align_motif() {
    MotifLibrary mlib(HELIX);
    Motif m = mlib.get_motif("HELIX.IDEAL");
    Motif m2 = m.copy();
    align_motif(m.ends()[0], m2.ends()[1], m2);
    m.to_pdb("node.0.pdb");
    m2.to_pdb("node.1.pdb");
    return 1;
    
}




int main(int argc, const char * argv[]) {
    if (test_sqlite3_connection() == 0)    { std::cout << "test_sqlite3_connection failed" << std::endl; }
    if (test_load_all() == 0)              { std::cout << "test_load_all failed" << std::endl; }
    if (test_get_motif() == 0)             { std::cout << "test_get_motif failed" << std::endl; }
    if (test_align_motif() == 0)           { std::cout << "test_get_motif failed" << std::endl; }

    return 0;
}
