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
    bool db_exists = filesystem::exists(db_path);
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
    int index = 1;
    for (auto &i : sg) {
      save_segment_to_database(sg[i], directory_name, db_name, sg_id);
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
    const Segment &segment, String directory_name, String database_name, int sg_id
    ) const {
    sqlite3 *db;
    String db_path = directory_name + "/" + database_name + ".db";
    bool db_exists = filesystem::exists(db_path);
    int conn_failure = sqlite3_open(db_path.c_str(), &db);
    if (conn_failure) {
      throw PersistenceException("Can't connect to database");
    }
    // Write segment to database, but not if it's already in there
    String segment_data = segment.to_str();
    if (!record_exists(db, segment_data, PERSISTENCE_VERSION)) {
      String sql = insert_segment_sql(segment_data, segment.get_name());
      char *error_msg = 0;
      sqlite3_exec(
        db,
        sql.c_str(),
        nullptr,
        0,
        &error_msg
      );
      // Once loaded back in, apply the rotation, then apply rotation
      auto coords_string = segment.get_end_center(0).get_str();
      auto rf_string = segment.get_end_ref_frame(0).get_str();
      String map_sql = insert_segment_map_sql(
        sqlite3_last_insert_rowid(db),
        sg_id,
        coords_string,
        rf_string
      );
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

  SegmentOP Persistence::retrieve_segment_from_database(String name, String db_path) const {
    sqlite3 *db;
    bool db_exists = filesystem::exists(db_path);
    int conn_failure = sqlite3_open(db_path.c_str(), &db);
    if (conn_failure) {
      throw PersistenceException("Can't connect to database");
    }
    String sql = "SELECT * FROM segments WHERE name = \"" + name + "\";";
    auto segment_data = SQLRecords(db, sql);
    auto seg = get_segment_from_str(segment_data[0][2]);
    return make_shared<Segment>(seg);
  }

  const SegmentGraphAllAtom Persistence::retrieve_segment_graph_from_database(String name, String db_path) const {
    sqlite3 *db;
    bool db_exists = filesystem::exists(db_path);
    int conn_failure = sqlite3_open(db_path.c_str(), &db);
    if (conn_failure) {
      throw PersistenceException("Can't connect to database");
    }
    String sql = "SELECT segments.data FROM segment_graphs "\
                 "JOIN segment_maps ON segment_maps.segment_id = segment_graphs.id " \
                 "JOIN segments ON segment_maps.segment_id = segments.id " \
                 "WHERE segment_graphs.name = \"" + name + "\";";
    auto segment_strings = SQLRecords(db, sql);
    std::vector<SegmentOP> segments;
    for (auto results : segment_strings) {
      for (auto result : results) {
        // Get segment from string
        auto seg = structure::all_atom::get_segment_from_str(result);
        segments.push_back(make_shared<Segment>(seg));
      }
    }
    segment_data_structure::SegmentGraphAllAtom sg;
    sg.add(segments[0]);
    for (int i = 1; i < segments.size(); i++) {
      sg.add(segments[i], 0, sg.get_end_name(0, 1));
      // Transverse the graph
    }
    return sg;
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
    // "Save the connectivity"
    String segment_map_table_sql = "CREATE TABLE segment_maps (" \
                                      "id               INTEGER PRIMARY KEY AUTOINCREMENT," \
                                      "segment_id       INT NOT NULL," \
                                      "segment_graph_id INT NOT NULL," \
                                      "is_root          INT," \
                                      "coord_data       TEXT," \
                                      "ref_frame        TEXT" \
                                    ");";
                                    // matrix code
                                    // Matrix data for rotation 3x3
                                    // Do we even need coordinates? They are recorded in
                                    // the segments table data (get_str output) column.
                                    // Need coordinates, translate rotation as well
                                    // root is first node without a parent (is not aligned)
                                    // look for get_roots
                                    // end[0] contains rotation and translation
                                    // align_segment when you get the data and are ready to apply it
                                    // "I think it's in segment"
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

  const String Persistence::insert_segment_map_sql(
    const int segment_id, const int segment_graph_id,
    const string coord_data, const string ref_frame
  ) const {
    String sql = "INSERT INTO segment_maps (segment_id, segment_graph_id, coord_data, ref_frame) VALUES (";
    sql += std::to_string(segment_id) + ", " + std::to_string(segment_graph_id);
    sql += ", \"" + coord_data + "\", \"" + ref_frame + "\");";
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
