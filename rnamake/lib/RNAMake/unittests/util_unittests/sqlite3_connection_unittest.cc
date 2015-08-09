//
//  sqlite3_connection_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "sqlite3_connection_unittest.h"
#include "util/sqlite3_connection.h"
#include "util/settings.h"

int
Sqlite3ConnectionUnittest::test_creation() {
    String path = resources_path()+"/motif_libraries_new/bp_steps.db";
    Sqlite3Connection conn(path);
    
    return 1;
}

int
Sqlite3ConnectionUnittest::test_count() {
    String path = resources_path()+"/motif_libraries_new/bp_steps.db";
    Sqlite3Connection conn(path);
    if(conn.count() != 529) { return 0; }
    
    return 1;
}

int
Sqlite3ConnectionUnittest::run() {
    if (test_creation() == 0)    {  std::cout << "test_creation failed" << std::endl; }
    if (test_count() == 0)       {  std::cout << "test_count failed" << std::endl; }
    return 0;
}