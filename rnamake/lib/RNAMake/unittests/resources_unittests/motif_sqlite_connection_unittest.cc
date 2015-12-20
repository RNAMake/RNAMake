//
//  motif_sqlite_connection_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_sqlite_connection_unittest.h"
#include "util/sqlite3_connection.h"
#include "resources/motif_sqlite_connection.h"
#include "util/random_number_generator.h"

int
MotifSqliteConnectionUnittest::test_creation() {
    String path = resources_path()+"/motif_libraries_new/bp_steps.db";
    MotifSqliteConnection conn(path);
    
    return 1;
}

int
MotifSqliteConnectionUnittest::test_next() {
    String path = resources_path()+"/motif_libraries_new/bp_steps.db";
    MotifSqliteConnection conn(path);
    conn.query("SELECT * from data_table");
    auto data = conn.next();
    int count = 0;
    while(data->data.size() != 0) {
        count++;
        data = conn.next();
    }
    std::cout << count << std::endl;
   

    return 1;

}

void
MotifSqliteConnectionUnittest::test_memory() {
    String path = resources_path()+"/motif_libraries_new/ideal_helices.db";
    MotifSqliteConnection conn(path);
    //auto total = conn.count();
    auto names = Strings{"HELIX.IDEAL", "HELIX.IDEAL.2", "HELIX.IDEAL.3", "HELIX.IDEAL.4"};
    auto rng = RandomNumberGenerator();
    
    int count = 0;
    for(int i = 0; i < 1000000; i++) {
        auto name = names[rng.randrange(names.size())];
        conn.query("SELECT * from data_table WHERE name='"+name+"'");
        auto data = conn.next();
        conn.clear();
        count += data->data.length();
        if(i % 100 == 0) { count = 0; }
    }
}


int
MotifSqliteConnectionUnittest::run() {
    //if (test_creation() == 0)             { std::cout << "test_creation failed" << std::endl; }
    //if (test_next() == 0)                 { std::cout << "test_next failed" << std::endl; }
    test_memory();
    return 0;
}