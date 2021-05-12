//
// Created by Hassan Abdelsamad on 4/30/21.
//


#ifndef RNAMAKE_TOOLS_H

#include "util/csv.h"


bool check_logs(const std::string f1, const std::string f2) {

    io::CSVReader<1> in1("Logs_orig.csv");
    io::CSVReader<1> in2("Logs.csv");
    in1.read_header(io::ignore_extra_column, "Message");
    in2.read_header(io::ignore_extra_column, "Message");
    bool found = false;
    std::string message1; std::string message2;

    while (in1.read_row(message1)) {
        found = false;
        while(!found && in2.read_row(message2) ) {
            std::cout<< "Message 1" << message1 << std::endl;
            std::cout<< "Message 2" << message2 << std::endl;
            if(message1 == message2){
                found = true;
            }
        }
        std::cout<< "Found? " << found << std::endl;

        if (!found) {
            return false;
        }
    }
    return true;
}

#define RNAMAKE_TOOLS_H

#endif //RNAMAKE_TOOLS_H
