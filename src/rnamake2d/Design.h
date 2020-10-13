#ifndef __RNAMAKE_DESIGN_H__
#define __RNAMAKE_DESIGN_H__

#include <iostream>
#include <utility>
#include <vector>

#include <base/types.h>

namespace rnamake2d {
   
    struct Design {

       static int design_ct;
       
       double score{0.};
       const String target;
       String sequence; // the current sequence
       String candiate; // the candidate
       const int id;

       explicit
       Design(String targ, String seq = "") : id(design_ct++),
                                                target(targ),
                                                sequence(seq){
           while(target.size() - sequence.size()) {
                sequence.push_back('N');
            }
       }

       void
       accept() {
            sequence = candiate;
       }

       void
       reject() {
            candiate = sequence;
       }


   };


}


using Designs = std::vector<rnamake2d::Design>;

#endif  // __RNAMAKE_DESIGN_H__
