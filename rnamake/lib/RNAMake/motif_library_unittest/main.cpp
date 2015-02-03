//
//  main.cpp
//  motif_library_unittest
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "sqlite3_connection.h"
#include "motif_type.h"
#include "motif.h"
#include "residue_type_set.h"

int test_sqlite3_connection() {
    Sqlite3Connection conn(HELIX);
    conn.query("SELECT * FROM motifs");
    conn.next();
    ResidueTypeSet rts;
    Strings values = conn.values();
    Motif m ( values[0], rts );
    std::cout << m.name() << std::endl;
    m.to_pdb("test.pdb");
    return 1;
}

int main(int argc, const char * argv[]) {
    if (test_sqlite3_connection() == 0)     { std::cout << "test_sqlite3_connection failed" << std::endl; }
    return 0;
}
