//
// Created by Joseph Yesselman on 12/30/17.
//

#include <util/sqlite/database.h>

namespace util::sqlite {

int Database::close() {
  int err = sqlite3_close(_db);
  _open = false;
  _db = nullptr;
  return err;
}

// open (connect) the database
void Database::_connect(const String &name, int flags) {
  _name = name;
  if (!std::filesystem::exists(name)) { _created = true; }
  int err = sqlite3_open_v2(_name.c_str(), &_db, flags, nullptr);
  _open = !err;
}
}// namespace util::sqlite