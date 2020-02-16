//
// Created by Joseph Yesselman on 2/16/20.
//

#ifndef RNAMAKE_NEW_CSV_H
#define RNAMAKE_NEW_CSV_H

#include <map>

#include <base/types.h>


namespace util {

class CSV_Reader {
public:
    CSV_Reader(
            String const & file_name):
            file_name_(file_name) {

    }

    ~CSV_Reader() {}

private:
    enum DataType {
        INT,
        FLOAT,
        STRING
    };

    union Data {
        int i;
        float f;
        String s;
    };


    class Row {
    public:

    private:
        std::map<String, int> col_names_;
    };

public:
    void
    start() {}



private:
    String file_name_;
    std::ifstream* in_;
};


class CSV_Writer {

};

}


#endif //RNAMAKE_NEW_CSV_H
