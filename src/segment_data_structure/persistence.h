//
// Created by Erik Whiting on 12/10/22.
//

#ifndef RNAMAKE_GRAPH_DATA_PERSISTENCE
#define RNAMAKE_GRAPH_DATA_PERSISTENCE

#include "segment_graph.h"

using namespace std;
using namespace segment_data_structure;
using namespace structure::all_atom;

namespace persistence {
  class Persistence {
    public:
      Persistence();
      // SEE NOTE IN persistence.cpp FOR WHY SegmentGraphAllAtom ISN'T CONST
      void save_to_database(SegmentGraphAllAtom&) const;
      void save_to_database(SegmentGraphAllAtom&, const String) const;
      void save_to_database(SegmentGraphAllAtom&, const String, const String) const;
      void save_segment_to_database(const SegmentOP&, String, String, int) const;
      void save_segment_map_to_database(String, sqlite3*) const;
      const int save_segment_graph_to_database(const String, const String, const String) const;
      SegmentGraphAllAtom retrieve_from_database(String) const;
      SegmentOP retrieve_segment_from_database(String, String) const;
      SegmentGraphAllAtom retrieve_segment_graph_from_database(String, String) const;
      vector<vector<String>> SQLRecords(sqlite3*, String) const;
    private:
      bool record_exists(sqlite3*, const String, const String) const;
      const String segment_table_sql() const;
      const String segment_graph_table_sql() const;
      const String segment_map_table_sql() const;
      void create_table(const String, sqlite3*) const;
      const String insert_segment_sql(const String, const String) const;
      const String insert_segment_graph_sql(const String) const;
      const String insert_segment_map_sql(const int, const int, math::Vector3s&) const;
  };

  class PersistenceException : public std::exception {
    private:
      char *message;

    public:
      PersistenceException(char *msg) : message(msg) {};
      PersistenceException(const char *msg) : message((char*)msg) {};
      char *what() { return message; }
  };
}

#endif // RNAMAKE_GRAPH_DATA_PERSISTENCE
