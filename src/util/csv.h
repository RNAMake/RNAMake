//
// Created by Joseph Yesselman on 2/16/20.
//

#ifndef RNAMAKE_NEW_CSV_H
#define RNAMAKE_NEW_CSV_H

#include <map>
#include <fstream>

#include <base/types.h>
#include <base/log.h>
#include <base/string.h>
#include <base/file_io.h>

namespace util {
namespace csv {

class CSVException : public std::runtime_error {
public:
    CSVException(): std::runtime_error("") {}

    /**
     * Standard constructor for CSVException
     * @param   message   Error message for csv parsing
     */
    CSVException(String const & message):
    std::runtime_error(message) {}


};

struct Data {
    int i;
    float f;
    String s;
};

class RowTopology {
public:
    inline
    RowTopology(
            DataTypes & data_types,
            std::map<String, int> & col_names):
            data_types_(std::move(data_types)),
            col_names_(std::move(col_names)) {}

public:
    int
    get_col_pos(
            String const & col_name) const {
        return _get_col_pos(col_name);
    }

    DataType
    get_col_data_type(
            String const & col_name) const {
        auto pos = _get_col_pos(col_name);
        return data_types_[pos];
    }

    bool
    does_col_exist(
            String const & name) const {
        return !(col_names_.find(name) == col_names_.end());
    }

private:
    int
    _get_col_pos(
            String const & col_name) const {
        if(col_names_.find(col_name) == col_names_.end()) {
            throw std::runtime_error(col_name + " col is not present in csv table");
        }
        return col_names_.at(col_name);
    }

private:

    DataTypes data_types_;
    std::map<String, int> col_names_;
};

class Row {
public:
    inline
    Row(
            RowTopology const & topology,
            std::vector<Data> & data):
            topology_(topology),
            data_(std::move(data)) {}

public: //getters
    int
    get_int_val(
            String const & col_name) const {
        auto pos = topology_.get_col_pos(col_name);
        auto data_type = topology_.get_col_data_type(col_name);
        if(data_type != DataType::INT) {
            throw std::runtime_error("attempting to access csv val: " + col_name + " but it is not a int");
        }
        return data_[pos].i;
    }

    float
    get_float_val(
            String const & col_name) const {
        auto pos = topology_.get_col_pos(col_name);
        auto data_type = topology_.get_col_data_type(col_name);
        if(data_type != DataType::FLOAT) {
            throw std::runtime_error("attempting to access csv val: " + col_name + " but it is not a float");
        }
        return data_[pos].f;
    }

    String const &
    get_string_val(
            String const & col_name) const {
        auto pos = topology_.get_col_pos(col_name);
        auto data_type = topology_.get_col_data_type(col_name);
        if(data_type != DataType::STRING) {
            throw std::runtime_error("attempting to access csv val: " + col_name + " but it is not a string");
        }
        return data_[pos].s;
    }

public: // wrapper to row topology

    bool
    does_col_exist(
            String const & name) const {
        return topology_.does_col_exist(name);
    }

private:
    RowTopology const & topology_;
    std::vector<Data> data_;
};

class Table {
public:
    explicit
    Table(
            RowTopology & topology):
            topology_(std::move(topology)) {}

public:
    typedef std::vector<Row>::iterator iterator;
    typedef std::vector<Row>::const_iterator const_iterator;

    iterator
    begin() { return rows_.begin(); }

    iterator
    end() { return rows_.end(); }

    const_iterator
    begin() const { return rows_.begin(); }

    const_iterator
    end() const { return rows_.end(); }

public:
    void
    add_row(
            std::vector<Data> & data) {
        rows_.emplace_back(topology_, data);
    }


    size_t
    num_rows() {
        return rows_.size();
    }

public: //row topology wrapper
    bool
    does_col_exist(
            String const & name) {
        return topology_.does_col_exist(name);
    }

private:
    std::vector<Row> rows_;
    RowTopology topology_;

};

typedef std::shared_ptr<Table> TableOP;

class Reader {
public:
    Reader() = default;

    ~Reader() = default;

public:
    TableOP
    read_csv(
            String const & file_name) {
        file_name_ = file_name;
        auto in = std::ifstream();

        // check to make sure file actually exists
        if(!base::file_exists(file_name)) {
            LOG_ERROR << "file does not exist: " << file_name; throw CSVException();
        }

        in.open(file_name);

        // get column names, which are always the first line
        auto col_name_map = std::map<String, int>();
        _get_col_names(in, col_name_map);
        int col_name_count = col_name_map.size();

        if(in.good()) {
            getline(in, line_);
        }
        else {
            LOG_ERROR << "there needs to be at least one data line in " << file_name; throw CSVException();
        }

        auto data_types = DataTypes(col_name_count);
        auto data = std::vector<Data>(col_name_count);

        _parse_data_col(line_, col_name_count, data_types, data);

        auto row_topology = RowTopology(data_types, col_name_map);
        auto table = std::make_shared<Table>(row_topology);
        table->add_row(data);

        while(in.good()) {
            getline(in, line_);
            if(line_.size() < 2) { break; }
            // need to reallocate since they were moved during construction of last row
            data_types = DataTypes(col_name_count);
            data = std::vector<Data>(col_name_count);
            _parse_data_col(line_, col_name_count, data_types, data);
            table->add_row(data);
        }

        return table;

    }

private: // helper functions

    void
    _get_col_names(
            std::ifstream & in,
            std::map<String, int> & col_name_map /* return */) {

        // is there are first line?
        if(in.good()) {
            getline(in, line_);
        }
        else {
            LOG_ERROR << "there are no lines in: " << file_name_; throw CSVException();
        }

        auto col_names = base::split_str_by_delimiter(line_, ",");
        if(col_names.empty()) {
            LOG_ERROR << "there are no columns in csv: " << file_name_ << ". Are you sure this is comma delimited?"; throw CSVException();
        }
        auto i = 0;
        for (auto const & name : col_names) {
            col_name_map[name] = i;
            i++;
        }
    }

    void
    _parse_data_col(
            String const & line,
            int col_name_count,
            DataTypes & data_types, /* return */
            std::vector<Data> & data /* return */) {

        auto data_cols = base::split_str_by_delimiter(line, ",");

        if(data_cols.size() != col_name_count) {
            LOG_ERROR << " there are not the same number of data columns as there are data names in: " << file_name_ << "\n";
            LOG_ERROR << "there are " << data_cols.size() << " data column and " << col_name_count << " column names";
            throw CSVException();
        }
        auto pos = 0;
        for(auto const & data_col : data_cols) {
            auto data_type = base::determine_string_data_type(data_col);
            if(data_type == DataType::STRING) {
                data[pos] = Data{-1, -1, data_col};
            }
            else if(data_type == DataType::INT) {
                data[pos] = Data{std::stoi(data_col), -1, ""};
            }
            else {
                data[pos] = Data{-1, std::stof(data_col), ""};
            }
            data_types[pos] = data_type;
            pos += 1;
        }

    }

private:
    String line_;
    String file_name_;

};


class Writer {

};

}
}


#endif //RNAMAKE_NEW_CSV_H

























