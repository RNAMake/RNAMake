//
// Created by Joseph Yesselman on 12/31/17.
//

#ifndef RNAMAKE_NEW_TABLE_DETAILS_H
#define RNAMAKE_NEW_TABLE_DETAILS_H

#include <base/types.hpp>
#include <utility>

#include <util/sqlite/field.hpp>

namespace util::sqlite {

class TableDetails {
public:
  struct ColumnDetails {
    String name;
    String type;
    bool is_primary;
  };

public:
  inline TableDetails(const String &name)
      : name_(name), columns_(std::vector<ColumnDetails>()) {}

  ~TableDetails() = default;

public: // iterators
  typedef std::vector<ColumnDetails>::const_iterator const_iterator;

  [[nodiscard]] const_iterator begin() const { return columns_.begin(); }
  [[nodiscard]] const_iterator end() const { return columns_.end(); }

public: // getters
  [[nodiscard]] inline const ColumnDetails &get_column(int i) const {
    return columns_[i];
  }

  [[nodiscard]] inline size_t size() const { return columns_.size(); }

  void add_column(String const &name, String const &type,
                  bool is_primary = false) {
    if (!_is_valid_sqlite_type(type)) {
      throw SqliteException("not a valid sqlite3 type: " + type);
    }
    if (_name_exists(name)) {
      throw SqliteException("col name already exists in table cannot repeat");
    }

    columns_.push_back(ColumnDetails{name, type, is_primary});
  }

  [[nodiscard]] bool has_primary_key() const {
    if (std::any_of(columns_.begin(), columns_.end(),
                    [](const ColumnDetails &c) { return c.is_primary; })) {
      return true;
    } else {
      return false;
    }
  }

public:
  [[nodiscard]] inline const String &name() const { return name_; }

  inline ColumnDetails const &operator[](Index i) const { return columns_[i]; }

private:
  bool _is_valid_sqlite_type(String const &type) {
    if (type == "TEXT" || type == "BLOB" || type == "INTEGER" ||
        type == "INT" || type == "REAL") {
      return true;
    } else {
      return false;
    }
  }

  bool _name_exists(String const &name) {
    for (auto const &col : columns_) {
      if (col.name == name) {
        return true;
      }
    }
    return false;
  }

private:
  String name_;
  std::vector<ColumnDetails> columns_;
};

typedef std::shared_ptr<TableDetails> TableDetailsOP;

} // namespace util::sqlite

#endif // RNAMAKE_NEW_TABLE_DETAILS_H
