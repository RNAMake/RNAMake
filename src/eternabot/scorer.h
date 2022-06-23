//
//  scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__scorer__
#define __RNAMake__scorer__

#include <stdio.h>

#include <eternabot/feature_generator.h>
#include <eternabot/strategy.h>
#include <eternabot/strategy/a_basic_test.h>
#include <eternabot/strategy/berex_test.h>
#include <eternabot/strategy/clear_plot.h>
#include <eternabot/strategy/direction_of_gc.h>
#include <eternabot/strategy/num_of_yellow.h>

namespace eternabot {

class Scorer {
public:
  Scorer()
      : generator_(FeatureGenerator()), strategies_(StrategyOPs()),
        weights_(Reals()) {

    strategies_.push_back(std::make_shared<ABasicTest>());
    strategies_.push_back(std::make_shared<CleanPlotStackCapsandSafeGC>());
    strategies_.push_back(std::make_shared<DirectionofGCPairsinMultiLoops>());
    strategies_.push_back(std::make_shared<BerexTest>());
    strategies_.push_back(
        std::make_shared<NumofYellowNucleotidesperLengthofString>());

    weights_.push_back(0.09281782);
    weights_.push_back(0.1250677);
    weights_.push_back(0.2156337);
    weights_.push_back(0.3661276);
    weights_.push_back(0.2230357);

    scores_ = Reals(strategies_.size());
    mean_ = 84.8005952381;
    stdev_ = 16.4725276237;
  }

  ~Scorer() {}

public:
  void setup(secondary_structure::PoseOP const &);

  float score_secondary_structure(secondary_structure::PoseOP const &);

public:
  Reals const &scores() { return scores_; }

private:
  FeatureGenerator generator_;
  FeaturesOP features_;
  StrategyOPs strategies_;
  Reals weights_, scores_;
  float mean_, stdev_, total_score_ = 0.0;
};

} // namespace eternabot

#endif /* defined(__RNAMake__scorer__) */
