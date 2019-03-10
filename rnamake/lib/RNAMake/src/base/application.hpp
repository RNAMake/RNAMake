//
//  application.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef application_hpp
#define application_hpp

#include <stdio.h>

#include "base/option.h"
#include "base/command_line_parser.hpp"

namespace base {

class Application {
public:
    inline
    Application() :
            options_(Options()),
            cl_options_(CommandLineOptions()),
            cl_parser_(CommandLineParser()) {}

    ~Application() {}

public:
    virtual
    void
    setup_options() = 0;


    virtual
    void
    parse_command_line(
            int argc,
            const char **argv) {

        cl_options_.parse_command_line(argc, argv);
        cl_parser_.assign_options(cl_options_, options_);

    }

    virtual
    void
    run() = 0;

public: //option wrappers

    template<typename T>
    void
    add_option(
            String const & name,
            T const & val,
            OptionType const & type,
            bool required = false) {

        options_.add_option(name, val, type);
        cl_options_.add_option(name, val, type, required);
    }

    void
    add_cl_options(
            Options const & opts,
            String prefix = "") {

        if (prefix.length() == 0) {
            cl_options_.add_options(opts);
            return;
        }

        for (auto const & opt : opts) {
            auto name = prefix + "." + opt->name();

            if (opt->type() == OptionType::STRING) {
                cl_options_.add_option(name, opt->get_string(), OptionType::STRING, false);
            } else if (opt->type() == OptionType::FLOAT) {
                cl_options_.add_option(name, opt->get_float(), OptionType::FLOAT, false);
            } else if (opt->type() == OptionType::INT) {
                cl_options_.add_option(name, opt->get_int(), OptionType::INT, false);
            } else if (opt->type() == OptionType::BOOL) {
                cl_options_.add_option(name, opt->get_bool(), OptionType::BOOL, false);

            }

        }
    }

    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }

    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }

    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }

    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }


protected:
    Options options_;
    CommandLineOptions cl_options_;
    CommandLineParser cl_parser_;

private:


};

}

#endif /* application_hpp */
