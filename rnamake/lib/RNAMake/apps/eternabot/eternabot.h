//
//  eternabot.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__eternabot__
#define __RNAMake__eternabot__

#include <stdio.h>
#include <base/application.hpp>

class EternabotApp : public base::Application {

public:

    EternabotApp() : base::Application() {}

    ~EternabotApp() {}

public:

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:
    struct Parameters {
        String seq, ss;
        int steps, n;
    };

private:
    Parameters parameters_;

};



#endif /* defined(__RNAMake__eternabot__) */
