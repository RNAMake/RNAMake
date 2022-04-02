//
//  Option.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 11/9/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Option__
#define __REDESIGNC__Option__

#include <stdio.h>

#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <typeinfo>
#include <vector>

// RNAMake Headers
#include "base/types.h"

namespace base {

enum class OptionType { BOOL = 1, INT = 2, STRING = 3, FLOAT = 4 };

class OptionException : public std::runtime_error {
 public:
  OptionException(String const &message) : std::runtime_error(message) {}
};

class Option {
 public:
  Option() {}

  Option(Option const &opt)
      : type_(opt.type_),
        name_(opt.name_),
        s_val_(opt.s_val_),
        f_val_(opt.f_val_),
        i_val_(opt.i_val_),
        b_val_(opt.b_val_) {}

  Option(String const &name, float const &value, OptionType const &type)
      : name_(name), type_(type) {
    f_val_ = value;

    if (type != OptionType::FLOAT) {
      throw OptionException("supplied a float option but with the wrong type");
    }
  }

  Option(String const &name, String const &value, OptionType const &type)
      : name_(name), type_(type) {
    s_val_ = value;

    if (type_ != OptionType::STRING) {
      throw OptionException("supplied a string option but with the wrong type");
    }
  }

  Option(String const &name, char const *value, OptionType const &type)
      : name_(name), type_(type) {
    s_val_ = String(value);

    if (type_ != OptionType::STRING) {
      throw OptionException("supplied a string option but with the wrong type");
    }
  }

  Option(char const *name, char const *value, OptionType const &type)
      : name_(String(name)), type_(type) {
    s_val_ = String(value);

    if (type_ != OptionType::STRING) {
      throw OptionException("supplied a string option but with the wrong type");
    }
  }

  Option(String const &name, bool const &value, OptionType const &type)
      : name_(name), type_(type) {
    if (type_ != OptionType::BOOL) {
      throw OptionException("calling wrong constructor for Option");
    }

    b_val_ = value;
  }

  Option(String const &name, int const &value, OptionType const &type)
      : name_(name), type_(type) {
    if (type_ == OptionType::INT) {
      i_val_ = value;
    } else if (type_ == OptionType::FLOAT) {
      f_val_ = float(value);
    } else if (type_ == OptionType::BOOL) {
      b_val_ = bool(value);
    } else if (type_ == OptionType::STRING) {
      throw OptionException("calling wrong constructor for Option");
    }
  }

  virtual ~Option() {}

 public:
  float get_float() {
    if (type_ == OptionType::STRING) {
      throw OptionException("attemped to get float but option is a string");
    }

    if (type_ == OptionType::BOOL) {
      throw OptionException("attemped to get float but option is a bool");
    }

    if (type_ == OptionType::INT) {
      return (float)i_val_;
    }

    return f_val_;
  }

  int get_int() {
    if (type_ == OptionType::STRING) {
      throw OptionException("attemped to get int but option is a string");
    }

    if (type_ == OptionType::BOOL) {
      throw OptionException("attemped to get int but option is a bool");
    }

    if (type_ == OptionType::FLOAT) {
      return f_val_;
    }

    return i_val_;
  }

  String const &get_string() {
    if (type_ == OptionType::INT) {
      throw OptionException("attemped to get string but option is a int");
    }
    if (type_ == OptionType::FLOAT) {
      throw OptionException("attemped to get string but option is a float");
    }
    if (type_ == OptionType::BOOL) {
      throw OptionException("attemped to get string but option is a bool");
    }

    return s_val_;
  }

  bool const &get_bool() {
    if (type_ == OptionType::INT) {
      throw OptionException("attemped to get bool but option is a int");
    }
    if (type_ == OptionType::FLOAT) {
      throw OptionException("attemped to get bool but option is a float");
    }
    if (type_ == OptionType::STRING) {
      throw OptionException("attemped to get bool but option is a string");
    }
    return b_val_;
  }

  void value(String const &);

  void value(char const *c_str) { return value(String(c_str)); }

  void value(float const &);

  void value(int const &);

  void value(bool const &);

 public:  // getters
  inline String const &name() const { return name_; }

  inline OptionType const &type() const { return type_; }

 public:
  inline String type_name() {
    if (type_ == OptionType::BOOL) {
      return "BOOL";
    } else if (type_ == OptionType::INT) {
      return "INT";
    } else if (type_ == OptionType::STRING) {
      return "STRING";
    } else if (type_ == OptionType::FLOAT) {
      return "FLOAT";
    } else {
      throw OptionException("unknown option type");
    }
  }

 protected:
  String name_;
  String s_val_;
  int i_val_;
  float f_val_;
  bool b_val_;

  OptionType type_;
};

typedef std::shared_ptr<Option> OptionOP;

class Options {
 public:
  Options() : options_(std::vector<OptionOP>()), locked_(false) {}

  Options(String const &name)
      : options_(std::vector<OptionOP>()), locked_(false) {}

  Options(Options const &opts)
      : locked_(opts.locked_), options_(std::vector<OptionOP>()) {
    for (auto const &opt : opts.options_) {
      auto new_opt = std::make_shared<Option>(*opt);
      options_.push_back(new_opt);
    }
  }

  ~Options() {}

 public:
  typedef std::vector<OptionOP>::iterator iterator;
  typedef std::vector<OptionOP>::const_iterator const_iterator;

  iterator begin() { return options_.begin(); }

  iterator end() { return options_.end(); }

  const_iterator begin() const { return options_.begin(); }

  const_iterator end() const { return options_.end(); }

 public:
  size_t size() { return options_.size(); }

  inline void lock_option_adding() { locked_ = true; }

  template <typename T>
  void add_option(String const &name, T const &val, OptionType const &type) {
    if (locked_) {
      throw OptionException("options are locked you cannot add any new ones");
    }

    auto opt = std::make_shared<Option>(name, val, type);
    options_.push_back(opt);
  }

  inline int get_int(String const &name) {
    auto opt = _find_option(name);
    return opt->get_int();
  }

  inline float get_float(String const &name) {
    auto opt = _find_option(name);
    return opt->get_float();
  }

  inline String const &get_string(String const &name) {
    auto opt = _find_option(name);
    return opt->get_string();
  }

  inline bool get_bool(String const &name) {
    auto opt = _find_option(name);
    return opt->get_bool();
  }

  inline bool has_option(String const &name) {
    for (auto const &opt : options_) {
      if (opt->name() == name) {
        return true;
      }
    }

    return false;
  }

  template <typename T>
  void set_value(String const &name, T const &val) {
    auto opt = _find_option(name);
    opt->value(val);
  }

 private:
  OptionOP const &_find_option(String const &name) {
    for (auto const &opt : options_) {
      if (opt->name() == name) {
        return opt;
      }
    }

    throw OptionException("cannot find option with name " + name);
  }

 private:
  bool locked_;
  String name_;
  std::vector<OptionOP> options_;
};

// interface
class OptionClass {
 public:
  virtual ~OptionClass() {}
};

}  // namespace base

#endif /* defined(__REDESIGNC__Option__) */
