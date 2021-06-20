//
// Created by Hassan Abdelsamad on 4/30/21.
//


#ifndef RNAMAKE_TOOLS_H

#include "util/csv.h"
#include "common.hpp"
#include <regex>

std::string get_arguments(std::string folder_name){
    //TODO Check if RNAMAKE variable exists
    auto base = std::getenv("RNAMAKE");

    std::string base_str(base);
    std::string args_path = base_str + "/integration_tests/" + folder_name + "/ARGS.txt";
    std::ifstream args_file(args_path);
    if (!args_file){
                FAIL("File does not exist, Make sure the ARGS.txt file exists in the test folder");
    }
    std::string args((std::istreambuf_iterator<char>(args_file)),
                     std::istreambuf_iterator<char>());

    args = std::regex_replace(args, std::regex("\\{RNAMAKE\\}"), base);

    return args;
}

std::string get_expected_path(std::string folder_name){
    auto base = std::getenv("RNAMAKE");

    std::string base_str(base);
    //Get expected logs
    std::string expected_path = base_str + "/integration_tests/" + folder_name + "/EXPECTED.csv";
    std::ifstream expected_file(expected_path);
    if (!expected_file){
                FAIL("File does not exist, Make sure the EXPECTED.csv file exists in the test folder");
    }

    return expected_path;
}

void mock_main(std::string args) {
    std::set_terminate(base::print_backtrace);
    auto app = DesignRNAScaffold();

    remove("logs.csv");

    base::init_logging_with_file(base::LogLevel::DEBUG);

    app.setup_options();

    app.app_.parse(args);
    app.run();
}

bool check_logs(const std::string &f1, const std::string &f2) {

    io::CSVReader<1> in1(f1.c_str());
    io::CSVReader<1> in2(f2.c_str());
    in1.read_header(io::ignore_extra_column, "Message");
    in2.read_header(io::ignore_extra_column, "Message");
    bool found = false;
    std::string message1; std::string message2;

    while (in1.read_row(message1)) {
        found = false;
        while(!found && in2.read_row(message2) ) {

            if(message1 == message2){
                found = true;
            }
        }

        if (!found) {
            FAIL("Expected log not found: " + message1);
            return false;
        }
    }
    return true;
}

#define RNAMAKE_TOOLS_H

#endif //RNAMAKE_TOOLS_H
