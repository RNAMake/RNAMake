//
// Created by Erik Whiting on 12/10/22.
//

#include <filesystem>
#include <sqlite3.h>

#include "persistence.h"

using namespace segment_data_structure;
using namespace std;

namespace persistence {
  void Persistence::save_to_database(SegmentGraphAllAtom sg) {
    Persistence::save_to_database(sg, "user_database");
  }

  void Persistence::save_to_database(SegmentGraphAllAtom sg, String name) {
    // Create directory if it doesn't exist
    filesystem::path current_dir = filesystem::current_path();
    String directory_name = name + "_dir";
    if (!filesystem::is_directory(directory_name)) {
      filesystem::create_directory(directory_name);
    }
    // Create database
    sqlite3 *db;
    String db_path = directory_name + "/" + name + ".db";
    bool db_exists = filesystem::exists(db_path);;
    int conn_failure = sqlite3_open(db_path.c_str(), &db);
    if (conn_failure) {
      throw PersistenceException("Can't connect to database");
    }
    if (!db_exists) {
      // Create the tables
      create_table(segment_table_sql(), db);
      create_table(segment_graph_table_sql(), db);
      create_table(segment_map_table_sql(), db);
    }
    // Upsert motif data

    sqlite3_close(db);
  }

  SegmentGraphAllAtom Persistence::retrieve_from_database(String name) {
    std::cout << "Retrieve stuff\n";
  }

  SegmentGraphAllAtom Persistence::retrieve_from_database(String name, String path) {
    std::cout << "Retrieve specific stuff\n";
  }

  void Persistence::create_table(String sql, sqlite3 *db) {
    char *error_msg = 0;
    sqlite3_exec(
      db,
      sql.c_str(),
      nullptr,
      0,
      &error_msg
    );
  }

  String Persistence::segment_table_sql() {
    String segment_table_sql = "CREATE TABLE segments (" \
                                  "id       INT PRIMARY KEY NOT NULL," \
                                  "name     TEXT," \
                                  "context  TEXT," \
                                  "end_id   TEXT," \
                                  "end_name TEXT," \
                                  "data     TEXT" \
                                ");";
    return segment_table_sql;
  }

  String Persistence::segment_graph_table_sql() {
    String segment_graph_table_sql = "CREATE TABLE segment_graphs (" \
                                        "id       INT PRIMARY KEY NOT NULL," \
                                        "name     TEXT NOT NULL," \
                                        "created  INT" \
                                      ");";
    return segment_graph_table_sql;
  }

  string Persistence::segment_map_table_sql() {
    String segment_map_table_sql = "CREATE TABLE segment_maps (" \
                                      "id               INT PRIMARY KEY NOT NULL," \
                                      "segment_id       INT NOT NULL," \
                                      "segment_graph_id INT NOT NULL," \
                                      "coord_x       FLOAT,"
                                      "coord_y       FLOAT,"
                                      "coord_z       FLOAT"
                                    ");";
                                    // Matrix data for rotation 3x3
    return segment_map_table_sql;
  }

  bool Persistence::record_exists(sqlite3 *db, String context, String name) {
    String sql = "SELECT id FROM segments WHERE" \
                  "context = " + context + " AND " \
                  "name = " + name + ";";
    vector<vector<String>> results = query_results(db, sql);
    return (results[0].size() > 0);
  }

  vector<vector<String>> Persistence::query_results(sqlite3 *db, String query) {
    vector<vector<String>> result;
    sqlite3_stmt *stmt;
    int rc = sqlite3_prepare(db, query.c_str(), -1, &stmt, 0);
    if (rc != SQLITE_OK) {
      throw PersistenceException(query.c_str());
    }

    rc = sqlite3_step(stmt);
    int cols = sqlite3_column_count(stmt);

    for (int i = 0; i < cols; i++) {
      result.push_back(vector<String>());
    }

    for (int i = 0; i < cols; i++) {
      result[i].push_back(string((char*)sqlite3_column_text(stmt, i)));
      sqlite3_step(stmt);
    }

    return result;
  }
}
