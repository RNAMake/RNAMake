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
      static void save_to_database(const SegmentGraphAllAtom);
      static void save_to_database(const SegmentGraphAllAtom, String);
      static SegmentGraphAllAtom retrieve_from_database(String);
      static SegmentGraphAllAtom retrieve_from_database(String, String);
    private:
      static bool record_exists(sqlite3*, String, String);
      static vector<vector<String>> query_results(sqlite3*, String);
      static String segment_table_sql();
      static String segment_graph_table_sql();
      static String segment_map_table_sql();
      static void create_table(String, sqlite3*);
      static Strings insert_statement(const SegmentGraphAllAtom&);
      static String pdb_data(const Segment&);
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
