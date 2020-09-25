#ifndef __RNAMAKE_DESIGN_H__
#define __RNAMAKE_DESIGN_H__

#include <utility>
#include <vector>


namespace rnamake2d {
   struct Design {

       double score{0.};
       const String target;
       String sequence;
       const int id;

       explicit
       Design(String  target) : id(count++),target(std::move(target)) {

       }

       void
       accept() {

       }

       void
       reject() {


       }

   private:
       static int count;

   };

}


using Designs = std::vector<rnamake2d::Design>;

#endif  // __RNAMAKE_DESIGN_H__