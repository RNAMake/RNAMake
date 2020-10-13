#ifndef __RNAMAKE_RULE_H__
#define __RNAMAKE_RULE_H__

#include <memory>

#include <base/types.h>
#include <rnamake2d/Design.h>

namespace rnamake2d {

   class Rule {
   private:
       String name_{"base_rule"};
       Reals params_{};

   public:
        virtual
        double
        score(Design const& ) const  = 0;

   public:
       virtual
       double
       score_mutation(Design const& ) const  = 0;

   public:
       String
       name()  const {
           return name_;
       }
   };



}

using RuleOP = std::shared_ptr<rnamake2d::Rule>;
#endif // __RNAMAKE_RULE_H__