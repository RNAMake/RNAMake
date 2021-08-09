
//
// Created by Joseph Yesselman on 2/15/17.
//

#include <primitives/pose.h>

namespace primitives {


  String
  get_dot_bracket_from_end_id(
          String const & end_id) {
      auto s = String("");
      auto spl = base::split_str_by_delimiter(end_id, "_");
      for(int i = 0; i < spl.size()-1; i+=2) {
          for(auto const & e : spl[i+1]) {
              if     (e == 'L') { s += '('; }
              else if(e == 'R') { s += ')'; }
              else if(e == 'U') { s += '.'; }
              else {
                  std::cout << end_id << std::endl;
                  std::cout << e << std::endl;
                  throw std::runtime_error(String(e, 1) + " is not a supported secondary structure symbol");
              }
          }
          if(i != spl.size()-2) { s += "&"; }
      }

      return s;
  }


}
