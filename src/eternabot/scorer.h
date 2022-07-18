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
      : _generator(FeatureGenerator()), _strategies(StrategyOPs()),
        _weights(Reals()) {

    _strategies.push_back(std::make_shared<ABasicTest>());
    _strategies.push_back(std::make_shared<CleanPlotStackCapsandSafeGC>());
    _strategies.push_back(std::make_shared<DirectionofGCPairsinMultiLoops>());
    _strategies.push_back(std::make_shared<BerexTest>());
    _strategies.push_back(
        std::make_shared<NumofYellowNucleotidesperLengthofString>());

    _weights.push_back(0.09281782);
    _weights.push_back(0.1250677);
    _weights.push_back(0.2156337);
    _weights.push_back(0.3661276);
    _weights.push_back(0.2230357);

    _scores = Reals(_strategies.size());
    _mean = 84.8005952381;
    _stdev = 16.4725276237;
  }

  ~Scorer() {}

public:
  void setup(secondary_structure::PoseOP const &);

  float score_secondary_structure(secondary_structure::PoseOP const &);

public:
  Reals const &scores() { return _scores; }

private:
  FeatureGenerator _generator;
  FeaturesOP _features;
  StrategyOPs _strategies;
  Reals _weights, _scores;
  float _mean, _stdev, _total_score = 0.0;
};

} // namespace eternabot

#endif /* defined(__RNAMake__scorer__) */
