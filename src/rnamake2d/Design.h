#ifndef __RNAMAKE_DESIGN_H__
#define __RNAMAKE_DESIGN_H__

#include <utility>
#include <vector>

#include <base/types.h>

namespace rnamake2d {
   
    struct Design {

       static int design_ct;
       
       double score{0.};
       const String target;
       String sequence;
       const int id;

       explicit
       Design(String  target) : id(design_ct++),target(std::move(target)) {

       }

       void
       accept() {

       }

       void
       reject() {


       }


   };


}


using Designs = std::vector<rnamake2d::Design>;

#endif  // __RNAMAKE_DESIGN_H__
