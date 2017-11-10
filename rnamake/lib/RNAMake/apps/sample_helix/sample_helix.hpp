//
//  sample_helix.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 6/7/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef sample_helix_hpp
#define sample_helix_hpp

#include <stdio.h>
#include "base/application.hpp"
#include "thermo_fluctuation/thermo_fluc_sampler.h"

class SampleHelixApp : public Application {
public:
    SampleHelixApp();

public: // application setups functions
    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

public:

    void
    run();

private:
    String
    _generate_structure(
            String const &);

    MotifOPs
    get_motifs_from_seq_and_ss(
            String const &,
            String const &);

private:
    ThermoFlucSampler sampler_;

};


#endif /* sample_helix_hpp */
