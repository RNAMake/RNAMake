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

namespace util {
namespace csv {

struct Data {
    int i;
    float f;
    String s;
};

class RowTopology {
public:
    inline
    RowTopology(
            std::vector<DataType> const & data_types,
            std::map<String, int> const & col_names):
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
        if(col_names_.find(name) == col_names_.end()) { return false;}
        else {                                          return true; }
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

    std::vector<DataType> data_types_;
    std::map<String, int> col_names_;
};

class Row {
public:
    inline
    Row(
            RowTopology const & topology,
            std::vector<Data> const & data):
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
    Table(
            RowTopology const & topology):
            topology_(topology) {}

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
            std::vector<Data> const & data) {
        rows_.push_back(Row(topology_, std::move(data)));
    }


    size_t
    num_rows() {
        return rows_.size();
    }

private:
    std::vector<Row> rows_;
    RowTopology topology_;

};

typedef std::shared_ptr<Table> TableOP;

class Reader {
public:
    Reader() {}

    ~Reader() {}

public:
    TableOP
    read_csv(
            String const & file_name) {
        auto in = std::ifstream();
        auto line = String();
        in.open(file_name);
        if(in.good()) {
            getline(in, line);
        }
        else {
            LOG_ERROR << "there are no lines in: " << file_name; exit(1);
        }

        auto col_names = base::split_str_by_delimiter(line, ",");
        auto col_name_map = std::map<String, int>();
        auto i = 0;
        for(auto const & name : col_names) {
            col_name_map[name] = i;
            i++;
        }

        if(in.good()) {
            getline(in, line);
        }
        else {
            LOG_ERROR << "there needs to be at least one data line in " << file_name; exit(1);
        }

        auto data_cols = base::split_str_by_delimiter(line, ",");

        if(data_cols.size() != col_names.size()) {
            LOG_ERROR << " there are not the same number of data columns as there are data names in: " << file_name;
            exit(1);
        }

        auto data_types = std::vector<DataType>();
        auto data = std::vector<Data>();
        for(auto const & data_col : data_cols) {
            auto data_type = base::determine_string_data_type(data_col);
            if(data_type == DataType::STRING) {
                data.push_back(Data{-1, -1, data_col});
            }
            else if(data_type == DataType::INT) {
                data.push_back(Data{std::stoi(data_col), -1, ""});
            }
            else {
                data.push_back(Data{-1, std::stof(data_col), ""});
            }
            data_types.push_back(data_type);
        }

        auto row_topology = RowTopology(data_types, col_name_map);
        auto table = std::make_shared<Table>(row_topology);
        table->add_row(data);


        while(in.good()) {
            getline(in, line);
            if(line.size() < 2) { break; }
            data = std::vector<Data>();
            for(auto const & data_col : data_cols) {
                auto data_type = base::determine_string_data_type(data_col);
                if(data_type == DataType::STRING) {
                    data.push_back(Data{-1, -1, data_col});
                }
                else if(data_type == DataType::INT) {
                    data.push_back(Data{std::stoi(data_col), -1, ""});
                }
                else {
                    data.push_back(Data{-1, std::stof(data_col), ""});
                }
            }
            table->add_row(data);
        }

        return table;

    }


};


class Writer {

};

}
}


#endif //RNAMAKE_NEW_CSV_H

























