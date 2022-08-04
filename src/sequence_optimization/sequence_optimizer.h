//
//  sequence_optimization.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//



#ifndef __RNAMake__sequence_optimizer__
#define __RNAMake__sequence_optimizer__

//#include "base/option.h"
#include "eternabot/sequence_designer.h"
#include "motif_data_structure/motif_graph.h"
#include "motif_data_structure/motif_tree.h"
#include <stdio.h>

/*
namespace sequence_optimization {

struct SequenceOptimizerResult {
  inline SequenceOptimizerResult(motif_data_structure::MotifTreeOP const &nmt,
                                 float nscore)
      : motif_tree(nmt), score(nscore) {}

  motif_data_structure::MotifTreeOP motif_tree;
  float score;
};

struct OptimizedSequence {
  inline OptimizedSequence(String const &nsequence,
                           float const &nclose_distance,
                           float const &neternabot_score)
      : sequence(nsequence), close_distance(nclose_distance),
        eternabot_score(neternabot_score) {}

  String sequence;
  float close_distance, eternabot_score;
};

typedef std::shared_ptr<SequenceOptimizerResult> SequenceOptimizerResultOP;
typedef std::shared_ptr<OptimizedSequence> OptimizedSequenceOP;
typedef std::vector<OptimizedSequenceOP> OptimizedSequenceOPs;

class SequenceOptimizer {
public:
  SequenceOptimizer() : _options(Options()) { setup_options(); }

  ~SequenceOptimizer() {}

public:
  OptimizedSequenceOPs
  get_optimized_sequences(motif_data_structure::MotifGraphOP &mg,
                          util::Uuid const &uuid_1, util::Uuid const &uuid_2,
                          int end_i, int end_j);

  String get_final_sequence(String const &, String const &);

public: // option wrappers
  inline Options &options() { return _options; }

  inline float get_int_option(String const &name) {
    return _options.get_int(name);
  }

  inline float get_float_option(String const &name) {
    return _options.get_float(name);
  }

  inline String get_string_option(String const &name) {
    return _options.get_string(name);
  }

  inline bool get_bool_option(String const &name) {
    return _options.get_bool(name);
  }

  inline bool has_option(String const &name) {
    return _options.has_option(name);
  }

  template <typename T>
  void set_option_value(String const &name, T const &val) {
    _options.set_value(name, val);
    update_var_options();
  }

protected:
  void setup_options();

  void update_var_options();

private:
  Options _options;
  motif_data_structure::MotifTreeOP _mt;
  eternabot::SequenceDesigner _designer;
  eternabot::SequenceDesignerResultOPs _designer_results;
  // options
  bool _sub_sequence;
  int _start;
  int _end;
  int _solutions;
};

} // namespace sequence_optimization
*/

#endif /* defined(__RNAMake__sequence_optimizer__) */


