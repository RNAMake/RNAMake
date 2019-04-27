//
// Created by Joseph Yesselman on 2019-04-23.
//

#ifndef RNAMAKE_NEW_THERMO_FLUC_GRAPH_LOGGING_H
#define RNAMAKE_NEW_THERMO_FLUC_GRAPH_LOGGING_H

#include <motif_data_structure/motif_state_graph.hpp>
#include <data_structure/graph_base.h>

namespace thermo_fluctuation {
namespace graph {
namespace logging {

class Logger {
public:
    Logger(
            String const & log_file) {
        _out = std::make_shared<std::ofstream>();
        _out->open(log_file);
    }

    ~Logger() {}

    virtual
    Logger *
    clone() const = 0;

public:
    virtual
    void
    setup(
            motif_data_structure::MotifStateGraphOP,
            data_structure::NodeIndexandEdge const &,
            data_structure::NodeIndexandEdge const &) = 0;

    virtual
    void
    log(
            motif_data_structure::MotifStateGraphOP,
            float) = 0;

protected:
    std::shared_ptr<std::ofstream> _out;
};

/**
 * Logs the rotation and translation of the two target basepairs to will overlap to generate a hit
 * Log order: target bp_1 origin, target bp_1 rotation, target bp_2 origin, target bp_2 rotation, score
 *
 */
class TargetBPInfoLogger : public Logger {
public:
    TargetBPInfoLogger(
            String const & log_file):
            Logger(log_file) {}

    ~TargetBPInfoLogger() {}

   Logger *
    clone() const { return new TargetBPInfoLogger(*this); };

public:

    void
    setup(
            motif_data_structure::MotifStateGraphOP msg,
            data_structure::NodeIndexandEdge const & start,
            data_structure::NodeIndexandEdge const & end) {
        *_out << "d_bp1,r_bp1,d_bp2,r_bp2,score" << std::endl;
        _start = start;
        _end = end;
    }

    void
    log(
            motif_data_structure::MotifStateGraphOP msg,
            float score) {


    }

private:
    data_structure::NodeIndexandEdge _start, _end;
};


}
}
}


#endif //RNAMAKE_NEW_THERMO_FLUC_GRAPH_LOGGING_H
