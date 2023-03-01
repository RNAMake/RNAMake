//
// Created by Erik Whiting on 12/10/22.
//

#include <filesystem>
#include <sqlite3.h>

#include "persistence.h"

using namespace segment_data_structure;
using namespace std;

int selectCallback(void*, int, char**, char**);

namespace persistence {
  // For ensuring we unwrap databases with the same
  // algorithm they were persisted with:
  const String PERSISTENCE_VERSION = "1";

  Persistence::Persistence() {};

  void Persistence::save_to_database(SegmentGraphAllAtom &sg) const {
    Persistence::save_to_database(sg, "user_database", "user_segment");
  }

  void Persistence::save_to_database(SegmentGraphAllAtom &sg, const String db_name) const {
    Persistence::save_to_database(sg, db_name, "user_segment");
  }

  void Persistence::save_to_database(SegmentGraphAllAtom &sg, const String db_name, const String sg_name) const {
    // Create directory if it doesn't exist
    filesystem::path current_dir = filesystem::current_path();
    String directory_name = db_name + "_dir";
    if (!filesystem::is_directory(directory_name)) {
      filesystem::create_directory(directory_name);
    }
    // Create database
    sqlite3 *db;
    String db_path = directory_name + "/" + db_name + ".db";
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
    // Insert graph data
    int sg_id = save_segment_graph_to_database(sg_name, directory_name, db_name);
    // Upsert motif data
    auto graph = sg.get_graph(); // When sg is a constant, this method doesn't work
    // Use get_node_data to get the node
    // Write a function for that in SegmentGraphAllAtom (whatever sg is)
    for (auto & i : sg) {
      auto segment = graph.get_node_data(i);
      save_segment_to_database(segment, directory_name, db_name, sg_id);
    }
    sqlite3_close(db);
  }

  const int Persistence::save_segment_graph_to_database(
    const String name, const String directory_name, const String database_name
  ) const {
    sqlite3 *db;
    String db_path = directory_name + "/" + database_name + ".db";
    bool db_exists = filesystem::exists(db_path);
    int conn_failure = sqlite3_open(db_path.c_str(), &db);
    if (conn_failure) {
      throw PersistenceException("Can't connect to database");
    }
    String sql = insert_segment_graph_sql(name);
    char *error_msg = 0;
    sqlite3_exec(
      db,
      sql.c_str(),
      nullptr,
      0,
      &error_msg
    );
    auto last_id = sqlite3_last_insert_rowid(db);
    return last_id;
  }

  void Persistence::save_segment_to_database(
    const SegmentOP &segment, String directory_name, String database_name, int sg_id
    ) const {
    sqlite3 *db;
    String db_path = directory_name + "/" + database_name + ".db";
    bool db_exists = filesystem::exists(db_path);
    int conn_failure = sqlite3_open(db_path.c_str(), &db);
    if (conn_failure) {
      throw PersistenceException("Can't connect to database");
    }
    // Write segment to database, but not if it's already in there
    String segment_data = segment->to_str();
    if (!record_exists(db, segment_data, PERSISTENCE_VERSION)) {
      String sql = insert_segment_sql(segment_data, segment->get_name());
      char *error_msg = 0;
      sqlite3_exec(
        db,
        sql.c_str(),
        nullptr,
        0,
        &error_msg
      );
      math::Vector3s segment_coordinates;
      for (auto const &r : *segment) {
        for (auto const &atom : r) {
          segment_coordinates.push_back(atom.get_coords());
        }
      }
      String map_sql = insert_segment_map_sql(sg_id, sqlite3_last_insert_rowid(db), segment_coordinates);
      save_segment_map_to_database(map_sql, db);
    }
  }

  void Persistence::save_segment_map_to_database(String sql, sqlite3 *db) const {
    char *error_msg = 0;
    sqlite3_exec(
      db,
      sql.c_str(),
      nullptr,
      0,
      &error_msg
    );
  }

  SegmentGraphAllAtom Persistence::retrieve_from_database(String name) const {
    std::cout << "Retrieve stuff\n";
  }

  SegmentGraphAllAtom Persistence::retrieve_from_database(String name, String path) const {
    std::cout << "Retrieve specific stuff\n";
  }

  void Persistence::create_table(const String sql, sqlite3 *db) const {
    char *error_msg = 0;
    sqlite3_exec(
      db,
      sql.c_str(),
      nullptr,
      0,
      &error_msg
    );
    if (error_msg != 0) {
      auto message = "Table creation error for: " + sql;
      throw PersistenceException(message.c_str());
    }
  }

  const String Persistence::segment_table_sql() const {
    String segment_table_sql = "CREATE TABLE segments (" \
                                  "id                   INTEGER PRIMARY KEY AUTOINCREMENT," \
                                  "name                 TEXT," \
                                  "data                 TEXT," \
                                  "persistence_version  TEXT" \
                                ");";
    return segment_table_sql;
  }

  const String Persistence::segment_graph_table_sql() const {
    String segment_graph_table_sql = "CREATE TABLE segment_graphs (" \
                                        "id       INTEGER PRIMARY KEY AUTOINCREMENT," \
                                        "name     TEXT NOT NULL" \
                                      ");";
    return segment_graph_table_sql;
  }

  const string Persistence::segment_map_table_sql() const {
    String segment_map_table_sql = "CREATE TABLE segment_maps (" \
                                      "id               INTEGER PRIMARY KEY AUTOINCREMENT," \
                                      "segment_id       INT NOT NULL," \
                                      "segment_graph_id INT NOT NULL," \
                                      "coord_x       FLOAT,"
                                      "coord_y       FLOAT,"
                                      "coord_z       FLOAT"
                                    ");";
                                    // matrix code
                                    // Matrix data for rotation 3x3
                                    // Do we even need coordinates? They are recorded in
                                    // the segments table data (get_str output) column.
    return segment_map_table_sql;
  }

  const String Persistence::insert_segment_sql(
    const String segment_data, const String name
  ) const {
    // Check existence of segment here, skip if it does
    String sql = "INSERT INTO segments (name, data, persistence_version) VALUES (";
    sql += "\"" + name + "\"" + ", ";
    sql += "\"" + segment_data + "\"" + ", ";
    sql += "\"" + PERSISTENCE_VERSION + "\"" + ");";
    return sql;
  }

  const String Persistence::insert_segment_graph_sql(const String name) const {
    String sql = "INSERT INTO segment_graphs (name) VALUES (";
    sql += "\"" + name + "\"" + ");";
    return sql;
  }

  const String Persistence::insert_segment_map_sql(const int segment_id, const int segment_graph_id, math::Vector3s& coordinates) const {
    String sql = "INSERT INTO segment_maps (segment_id, segment_graph_id) VALUES (";
    sql += std::to_string(segment_id) + ", " + std::to_string(segment_graph_id) + ");";
    // Ask about coordinates here
    return sql;
  }

  bool Persistence::record_exists(
    sqlite3 *db, const String segment_data, const String persistence_version
  ) const {
    String sql = "SELECT id FROM segments WHERE " \
                  "data = \"" + segment_data + "\" AND " \
                  "persistence_version = \"" + persistence_version + "\";";
    auto results = SQLRecords(db, sql);
    if (results.size() > 0) {
      std::cout << "Segment with this data already exists\n";
      return true;
    } else {
      return false;
    }
  }

  vector<vector<String>> Persistence::SQLRecords(sqlite3 *db, String sql) const {
    vector<vector<String>> records;
    char *errmsg;
    int ret = sqlite3_exec(db, sql.c_str(), selectCallback, &records, &errmsg);
    if (ret != SQLITE_OK) {
      std::cerr << "Error in select statement " << sql << "[" << errmsg << "]\n";
    } else {
      std::cerr << records.size() << " records returned.\n";
    }
    return records;
  }

}

int selectCallback(void *p_data, int num_fields, char **p_fields, char **p_col_names) {
  vector<vector<String>>* records = static_cast<vector<vector<String>>*>(p_data);
  try {
    records->emplace_back(p_fields, p_fields + num_fields);
  } catch (...) {
    // abort select on failure, don't let exception propogate thru sqlite3 call-stack
    return 1;
  }
  return 0;
}
