//
//  cl_option.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/cl_option.h"

#include "base/option.h"

namespace base {

void CommandLineOptions::add_options(Options const &opts) {
  for (auto const &opt : opts) {
    auto cl_opt = std::make_shared<CommandLineOption>(*opt);
    options_.push_back(cl_opt);
  }
}

void CommandLineOptions::_set_option(CommandLineOptionOP const &cl_opt,
                                     String const &value, bool is_bool) {
  if (cl_opt->filled()) {
    CommandLineOptionException(
        "commandline option has already has been set!: " + cl_opt->name());
  }

  cl_opt->filled(true);

  if (cl_opt->type() == OptionType::STRING) {
    try {
      /*bool fail = 0;
      for(auto const & e : value) {
          if(std::isdigit(e) || e == '.') {
              continue;
          }
          else { fail = 1; break; }
      }

      if(!fail) {
          throw CommandLineOptionException(
              cl_opt->name() + " is a STRING option but it can successfully
      convert it float");
      }*/
    } catch (std::invalid_argument) {
    }

    cl_opt->value(String(value));
  } else if (cl_opt->type() == OptionType::FLOAT) {
    try {
      auto f = std::stof(value);
      cl_opt->value(f);
    } catch (std::invalid_argument) {
      throw CommandLineOptionException(
          cl_opt->name() +
          " is a FLOAT option but cannot successfully convert it float");
    }

  } else if (cl_opt->type() == OptionType::INT) {
    try {
      auto i = std::stoi(value);
      cl_opt->value(i);
    } catch (std::invalid_argument) {
      throw CommandLineOptionException(
          cl_opt->name() +
          " is a INT option but cannot successfully convert it int");
    }
  } else if (cl_opt->type() == OptionType::BOOL) {
    if (!is_bool) {
      throw CommandLineOptionException(
          cl_opt->name() +
          " is a BOOL option must be supplied like \"-bool_option\" not " +
          "-bool_option 1");
    }

    cl_opt->value(true);
  }
}

void CommandLineOptions::parse_command_line(int const argc, char const **argv) {
  String key = "";
  CommandLineOptionOP cl_opt;

  for (int i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      continue;
    }

    key = String(argv[i]);
    key = key.substr(1);

    cl_opt = _find_option(key);

    if (argc != i + 1 && argv[i + 1][0] != '-') {
      _set_option(cl_opt, String(argv[i + 1]), false);
    } else {
      if (cl_opt->type() != OptionType::BOOL) {
        throw CommandLineOptionException(cl_opt->name() + " is a " +
                                         cl_opt->type_name() + " but " +
                                         "is set as a BOOL");
      }

      _set_option(cl_opt, "1", true);
    }
  }

  for (auto const &opt : options_) {
    if (!opt->filled() && opt->required()) {
      throw CommandLineOptionException(
          opt->name() + " is a required option and was not supplied");
    }
  }
}

}  // namespace base
