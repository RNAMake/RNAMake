#ifndef __RNAMAKE_NEMO_SAMPLER_H__
#define __RNAMAKE_NEMO_SAMPLER_H__

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <climits>
#include <unistd.h>
#include <limits>
#include <filesystem>

#include <rnamake2d/design.h>
#include <base/types.h>

extern "C" {
    #include <RNAstruct.h>
    #include <fold_vars.h>
    #include <fold.h>
    #include <params.h>
    #include <part_func.h>
    #include <utils.h>
    #include <convert_epars.h>
    #include <read_epars.h>
    #include <MEA.h>
}

namespace rnamake2d {

    class NemoSampler {
    private:
        // general variables
        const double worst_score = -100000.0;
        bool   standard_nmcs = false;
        int    iter = 2500;
        int    npairs = 0;
        bool   found = false;
        bool   force_boost = false;
        int    rng_seed = 1996;
        int    len = 0;
        int    bpd = std::numeric_limits<int>::max();
        int    current_ = -1;
    private:
        // strings
        char*  target = nullptr;
        char*  start = nullptr;
        char*  seed = nullptr;
        char*  position = nullptr;
        char*  closest_seq = nullptr;
    private:
        // vectors
        std::vector <short> pt;
        std::vector <short> mt ;
        std::vector <int> lt;
        std::vector <int> jct ;
        std::vector <int> smap;
        std::vector <int> shuffle ;

    public:
        NemoSampler()   {
            init_rand();
            srand48( rng_seed ^ 0x5A5A5A5A );
        }

    public:
        ~NemoSampler() {
            reset_();
        }

    public:
        void
        initialize_design(Design& );

    public:
        void
        mutate(Design& )  ;

    public:
        int
        current() const {
            return current_;
        }

    public:
        void
        current (int curr) {
            current_ = curr;
        }

    private:
        bool
        play_move(char *position, const char *bases,  int z ) ;    private:

    private:
        int
        nemo_main( Design& ) ;

    private:
        void
        reset_() {
            // reset variables
            bpd = std::numeric_limits<int>::max();
            current_ = -1;
            // free dynamic Cstrings
            if( target ) free( target );
            if( start ) free( start );
            if( seed ) free( seed );
            if( closest_seq )  free( closest_seq );
            if( position ) free( position );
            // clear out the vectors
            if( !pt.empty() ) pt.clear();
            if( !mt.empty() ) mt.clear();
            if( !lt.empty() ) lt.clear();
            if( !jct.empty() ) jct.clear();
            if( !smap.empty() ) smap.clear();
            if( !shuffle.empty() ) shuffle.clear();
        }

    private:
        double
        sample( char * )   ;

    private:
        void
        config_eterna() {
            // convert_parameter_file("vrna185.par", "vrna185x2.par", VRNA_CONVERT_OUTPUT_ALL);
            if(!std::filesystem::exists("vrna185x.par")) {
                std::cout<<"Error: File \"vrna185x.par\" could not be found. Unable to update parameters. Continuing...\n";
                return;
            }
            read_parameter_file("vrna185x.par");
            dangles = 1;
            update_fold_params();
        }

    private:
        double
        nested(char *position,  int ) ;

    };
}

#endif // __RNAMAKE_NEMO_SAMPLER_H__