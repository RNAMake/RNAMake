
//
//  string.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//
#include <algorithm>
#include <functional>

// RNAMake Headers
#include "base/string.h"

namespace base {

Strings split_str_by_delimiter(String s, String delimiter) {
  String token;
  std::vector<String> tokens;
  size_t pos = 0;
  while ((pos = s.find(delimiter)) != String::npos) {
    token = s.substr(0, pos);
    tokens.push_back(token);
    s.erase(0, pos + delimiter.length());
  }

  if (s.length() > 0) {
    tokens.push_back(s);
  }

  return tokens;
}

Strings tokenize_line(String const& raw_line) {
  // note that this is specific to the CIF format grammar. DO NOT use this
  // elsewhere
  auto tokens = Strings{};
  auto token = String{};
  const auto whitespace = std::string(" \f\v\t\r\n");
  const auto line_len = raw_line.size();

  for (auto it = 0; it < line_len;) {
    char ch = raw_line[it];
    const auto is_ws = whitespace.find(ch) != std::string::npos;

    if ((ch == ';' || ch == '\'') && token.empty()) {
      const auto matching = raw_line.find(ch, it + 1);
      tokens.push_back(raw_line.substr(it + 1, matching - it - 1));
      it = matching + 1;
      token = "";

    } else if (is_ws || (it == (line_len - 1))) {
      if (!token.empty()) {
        if (!is_ws && it == line_len - 1) {
          token += ch;
        }
        tokens.push_back(token);
        token = "";
      }
      ++it;
    } else if (!is_ws) {
      token += ch;
      ++it;
    }
  }

  return tokens;
}

String join_by_delimiter(Strings const& strs, String const& delimiter) {
  String return_s;
  int i = 0;
  for (auto const& s : strs) {
    return_s += s;
    if (i != strs.size() - 1) {
      return_s += delimiter;
    }
  }

  return return_s;
}

String filename(String const& path) {
  Strings path_spl = split_str_by_delimiter(path, "/");
  return path_spl.back();
}

String base_dir(String const& path) {
  auto path_spl = split_str_by_delimiter(path, "/");
  auto base_path = String();
  for (int i = 0; i < path_spl.size() - 1; i++) {
    base_path += path_spl[i] + "/";
  }
  if (base_path == "") {
    base_path = "./";
  }

  return base_path;
}

bool is_number(String const& s) {
  for (auto const& c : s) {
    if (!std::isdigit(c)) {
      return false;
    }
  }
  return true;

  // return !s.empty() && std::find_if(s.begin(),
  //                                   s.end(), [](char c) { return
  //                                   !std::isdigit(c); }) == s.end();
}

// adapted from
// https://www.geeksforgeeks.org/check-given-string-valid-number-integer-floating-point/
DataType determine_string_data_type(String const& str) {
  int i = 0, j = str.length() - 1;

  // Handling whitespaces
  while (i < str.length() && str[i] == ' ') {
    i++;
  }
  while (j >= 0 && str[j] == ' ') {
    j--;
  }

  if (i > j) {
    return DataType::STRING;
  }

  // if string is of length 1 and the only
  // character is not a digit
  if (i == j && !(str[i] >= '0' && str[i] <= '9')) {
    return DataType::STRING;
  }

  // If the 1st char is not '+', '-', '.' or digit
  if (str[i] != '.' && str[i] != '+' && str[i] != '-' &&
      !(str[i] >= '0' && str[i] <= '9')) {
    return DataType::STRING;
  }

  // To check if a '.' or 'e' is found in given
  // string. We use this flag to make sure that
  // either of them appear only once.
  bool flagDotOrE = false;
  bool has_dot = false;
  int char_count = 0;

  for (; i <= j; i++) {
    // If any of the char does not belong to
    // {digit, +, -, ., e}
    if (str[i] != 'e' && str[i] != '.' && str[i] != '+' && str[i] != '-' &&
        !(str[i] >= '0' && str[i] <= '9')) {
      return DataType::STRING;
    }

    if (str[i] == '.') {
      has_dot = true;
      // checks if the char 'e' has already
      // occurred before '.' If yes, return 0.
      if (flagDotOrE == true) {
        return DataType::STRING;
      }

      // If '.' is the last character.
      if (i + 1 > str.length()) {
        return DataType::STRING;
      }

      // if '.' is not followed by a digit.
      if (!(str[i + 1] >= '0' && str[i + 1] <= '9')) {
        return DataType::STRING;
      }
    }

    else if (str[i] == 'e') {
      // set flagDotOrE = 1 when e is encountered.
      flagDotOrE = true;

      // if there is no digit before 'e'.
      if (!(str[i - 1] >= '0' && str[i - 1] <= '9')) {
        return DataType::STRING;
      }

      // If 'e' is the last Character
      if (i + 1 > str.length()) {
        return DataType::STRING;
      }

      // if e is not followed either by
      // '+', '-' or a digit
      if (str[i + 1] != '+' && str[i + 1] != '-' &&
          (str[i + 1] >= '0' && str[i] <= '9')) {
        return DataType::STRING;
      }
    }

    if (char_count > 0 && (str[i] == '+' || str[i] == '-')) {
      return DataType::STRING;
    }

    char_count++;
  }

  if (flagDotOrE || has_dot) {
    return DataType::FLOAT;
  } else {
    return DataType::INT;
  }
}

/**
 * @brief Left Trim
 * Trims whitespace from the left end of the provided String
 * @param[out] s The String to trim
 * @return The modified String&
 */
String& ltrim(String& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
          }));

  return s;
}

/**
 * @brief Right Trim
 * Trims whitespace from the right end of the provided String
 * @param[out] s The String to trim
 * @return The modified String&
 */
String& rtrim(String& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); })
              .base(),
          s.end());

  return s;
}

/**
 * @brief Trim
 * Trims whitespace from both ends of the provided String
 * @param[out] s The String to trim
 * @return The modified String&
 */
String& trim(String& s) { return ltrim(rtrim(s)); }

bool is_char_in_string(char c, String const& s) {
  int pos = s.find(c);
  if (pos == std::string::npos) {
    return false;
  } else {
    return true;
  }
}

String quoted_string(String const& s) { return String("'") + s + String("'"); }

String string_map_to_string(StringStringMap const& ssm) {
  auto s = String("");
  for (auto const& kv : ssm) {
    s += kv.first + "=" + kv.second + " ";
  }
  return s;
}

String& replace_all(String& context, String const& from, String const& to) {
  std::size_t lookHere = 0;
  std::size_t foundHere;
  while ((foundHere = context.find(from, lookHere)) != String::npos) {
    context.replace(foundHere, from.size(), to);
    lookHere = foundHere + to.size();
  }
  return context;
}

}  // namespace base
