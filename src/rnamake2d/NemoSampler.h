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

#include <rnamake2d/Design.h>
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
        const double worst_score = -100000.0;
        int    verbosity = 0;
        bool   standard_nmcs = false;
        int    iter = 2500;
        char*  target = nullptr;
        char*  start = nullptr;
        char*  seed = nullptr;
        char*  position = nullptr;
        char*  closest_seq = nullptr;
        std::vector<short> pt;
        //short* pt = nullptr;       // Pair Table
        std::vector<short> mt ;       // Mismatch Table
        //short* mt = nullptr;       // Mismatch Table
        std::vector<int> lt;       // Loop Table
        //int*   lt = nullptr;       // Loop Table
        std::vector<int> jct ;      // Junction Table
        //int*   jct = nullptr;      // Junction Table
        std::vector<int> smap;     // Strength Map
        //int*   smap = nullptr;     // Strength Map
        int    npairs = 0;
        std::vector<int> shuffle ;
        //int*   shuffle = nullptr;
        bool   found = false;
        bool   force_boost = false;
        int    rng_seed = 1996;
        int    len = 0;
        int    bpd = std::numeric_limits<int>::max();
        int    current_ = -1;

    public:
        NemoSampler()   {
            init_rand();
            srand48( rng_seed ^ 0x5A5A5A5A );
        }

    public:
        ~NemoSampler() {
            reset_(); return;
            shuffle.clear();
            //if( shuffle ) free( shuffle );
            smap.clear();
            //if( smap ) free( smap );
            jct.clear();
            //if( jct ) free( jct );
            lt.clear();
            //if( lt ) free( lt );
            mt.clear();
            //if( mt ) free( mt );
            pt.clear();
            //if( pt ) free( pt );
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
            // make sure that when I give over the information the copies are "owned" by the sampler
            if( target ) free( target );
            if( start ) free( start );
            if( seed ) free( seed );
            pt.clear();
            //if( pt ) free( pt );       // Pair Table
            mt.clear();
            //if( mt ) free( mt );       // Mismatch Table
            lt.clear();       // Loop Table
            //if( lt ) free( lt );       // Loop Table
            jct.clear();      // Junction Table
            //if( jct ) free( jct );      // Junction Table
            smap.clear();
            //if( smap ) free( smap );     // Strength Map
            shuffle.clear();
            //if( shuffle  ) free( shuffle );
            if( closest_seq )  free( closest_seq );
            if( position ) free( position );

            bpd = std::numeric_limits<int>::max();
            current_ = -1;
        }

    private:
        double
        sample( char * )   ;

    private:
        void
        config_eterna(void) {
            // convert_parameter_file("vrna185.par", "vrna185x2.par", VRNA_CONVERT_OUTPUT_ALL);
            if(!std::filesystem::exists("vrna185x.par")) {
                std::cout<<"Error: File \"vrna185x.par\" could not be found. Unable to update parameters. Continuing...\n";
                return;
            }
            read_parameter_file("vrna185x.par");
            dangles = 1;
            update_fold_params(); // useful ?
        }

    private:
        double
        nested(char *position,  int ) ;

    };
}

#endif // __RNAMAKE_NEMO_SAMPLER_H__