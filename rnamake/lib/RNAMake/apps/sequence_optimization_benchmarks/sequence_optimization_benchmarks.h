//
// Created by Joseph Yesselman on 3/7/19.
//

#ifndef TEST_SEQUENCE_OPTIMIZATION_BENCHMARKS_H
#define TEST_SEQUENCE_OPTIMIZATION_BENCHMARKS_H

#include "base/application.hpp"
#include "util/basic_io.hpp"
#include "util/steric_lookup.hpp"
#include "motif_data_structure/motif_state_graph.hpp"
#include "motif_search/motif_state_monte_carlo.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"

struct NodeIndexandEdge {
    int ni; //node index
    int ei; //end index
};

class Timer {
public:
    Timer() {}

    ~Timer() {}

public:
    void
    start() {
        start_ =  std::chrono::steady_clock::now();
    }

    double
    end() {
        auto end = std::chrono::steady_clock::now();
        auto time = std::chrono::duration<double, std::milli>(end - start_).count();
        return time;
    }

private:
    std::chrono::steady_clock::time_point start_;
};

struct SequenceOptProblem {
public:
    inline
    SequenceOptProblem(
            motif_data_structure::MotifStateGraphOP n_msg,
            NodeIndexandEdge const & n_start,
            NodeIndexandEdge const & n_end,
            String const & n_start_name,
            String const & n_end_name,
            util::StericLookupNewOP n_lookup,
            bool n_target_an_aligned_end):
            msg(n_msg),
            start(n_start),
            end(n_end),
            start_name(n_start_name),
            end_name(n_end_name),
            lookup(n_lookup),
            target_an_aligned_end(n_target_an_aligned_end) {}

public:
    motif_data_structure::MotifStateGraphOP msg;
    NodeIndexandEdge start, end;
    String start_name, end_name;
    util::StericLookupNewOP lookup;
    bool target_an_aligned_end;
};

typedef std::shared_ptr<SequenceOptProblem> SequenceOptProblemOP;

class SequenceOptProblemFactory {
public:
    SequenceOptProblemFactory() {}

    virtual
    ~SequenceOptProblemFactory() {}

public:
    virtual
    SequenceOptProblemOP
    get_problem() = 0;

protected:
    void
    _setup_sterics(
            motif_data_structure::MotifStateGraphOP msg,
            util::StericLookupNewOP lookup) {

        auto beads = math::Points();
        for (auto & n : *msg) {
            for(auto const & b : n->data()->cur_state->beads()) {
                beads.push_back(b);
            }
        }
        lookup->add_points(beads);
    }
};

typedef std::shared_ptr<SequenceOptProblemFactory> SequenceOptProblemFactoryOP;

/*class TTRProblemFactory : public SequenceOptProblemFactory {
public:
    TTRProblemFactory():
            SequenceOptProblemFactory(),
            has_setup_(false) {}

    ~TTRProblemFactory() {}

public:
    SequenceOptProblemOP
    get_problem() {
        auto ttr = resources::Manager::instance().motif("GAAA_tetraloop", "", "A229-A245");
        auto bp_step_1 = resources::Manager::instance().bp_step("GC_LL_GC_RR");
        auto bp_step_2 = resources::Manager::instance().bp_step("CG_LL_CG_RR");
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->add_state(bp_step_1->get_state());
        msg->add_state(bp_step_2->get_state());
        msg->add_state(ttr->get_state());

        if (!has_setup_) {
            lookup_ = std::make_shared<util::StericLookup>();
            this->_setup_sterics(msg, lookup_);

            start_name_ = "A222-A251";
            end_name_ = "A149-A154";

            auto start_end_pos = msg->get_node(2)->data()->get_end_index(start_name_);
            auto end_end_pos = msg->get_node(2)->data()->get_end_index(end_name_);

            start_ = NodeIndexandEdge{2, start_end_pos};
            end_ = NodeIndexandEdge{2, end_end_pos};

            target_an_aligned_end_ = false;
            if(end_end_pos == msg->get_node(2)->data()->block_end_add()) {
                target_an_aligned_end_ = true;
            }
        }

        has_setup_ = true;

        auto p = std::make_shared<SequenceOptProblem>(msg, start_, end_, start_name_, end_name_,
                                                      lookup_, target_an_aligned_end_);
        return p;
    }

private:
    bool has_setup_;
    util::StericLookupOP lookup_;
    NodeIndexandEdge start_, end_;
    String start_name_, end_name_;
    bool target_an_aligned_end_;
};

class Add3WAYProblemFactory : public SequenceOptProblemFactory {
public:
    Add3WAYProblemFactory():
            SequenceOptProblemFactory(),
            has_setup_(false),
            mf_(motif::MotifFactory()) {}

    ~Add3WAYProblemFactory() {}

public:
    SequenceOptProblemOP
    get_problem() {
        if (!has_setup_) {
            String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/apt_new_interface/resources/";
            resources::Manager::instance().add_motif(base_path+"pRNA_3WJ.pdb", "prna");

            auto path = base::base_dir() + "/rnamake/lib/RNAMake/apps/sequence_optimization_benchmarks/resources/start.pdb";
            auto scaffold_rm = resources::Manager::instance().get_structure(path, "scaffold", 3);
            auto scaffold_m = std::make_shared<motif::Motif>(*scaffold_rm);
            mf_._setup_secondary_structure(scaffold_m);
            resources::Manager::instance().register_motif(scaffold_m);

            prna_ = resources::Manager::instance().motif("prna", "", "A7-C10");
            scaffold_ = resources::Manager::instance().motif("scaffold");
        }

        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->set_option_value("sterics", false);
        msg->add_state(scaffold_->get_state());
        msg->add_state(resources::Manager::instance().motif_state("HELIX.IDEAL.1"), 0, 2);
        msg->add_state(prna_->get_state());

        if (!has_setup_) {
            lookup_ = std::make_shared<util::StericLookup>();
            this->_setup_sterics(msg, lookup_);

            start_name_ = "B14-C7";
            end_name_ = "A41-A87";

            auto start_end_pos = msg->get_node(2)->data()->get_end_index(start_name_);
            auto end_end_pos = msg->get_node(0)->data()->get_end_index(end_name_);

            start_ = NodeIndexandEdge{2, start_end_pos};
            end_ = NodeIndexandEdge{0, end_end_pos};

            target_an_aligned_end_ = false;
            if(end_end_pos == msg->get_node(0)->data()->block_end_add()) {
                target_an_aligned_end_ = true;
            }
        }

        has_setup_ = true;

        auto p = std::make_shared<SequenceOptProblem>(msg, start_, end_, start_name_, end_name_,
                                                      lookup_, target_an_aligned_end_);
        return p;

    }

private:
    bool has_setup_;
    util::StericLookupOP lookup_;
    NodeIndexandEdge start_, end_;
    String start_name_, end_name_;
    bool target_an_aligned_end_;
    motif::MotifFactory mf_;
    motif::MotifOP scaffold_, prna_;

};
*/

class RibosomeTetherProblemFactory : public SequenceOptProblemFactory {
public:
    RibosomeTetherProblemFactory():
            SequenceOptProblemFactory(),
            has_setup_(false) {}

    ~RibosomeTetherProblemFactory() {}

public:
    SequenceOptProblemOP
    get_problem() {
        if(! has_setup_) {
            auto path = base::base_dir() + "/rnamake/lib/RNAMake/apps/sequence_optimization_benchmarks/resources/";
            resources::Manager::instance().add_motif(path + "short.out.1.pdb", "scaffold", util::MotifType::TWOWAY);
            scaffold_ = resources::Manager::instance().motif("scaffold", "", "A1019-A3915");
        }

        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->set_option_value("sterics", false);
        msg->add_state(scaffold_->get_state());

        if(! has_setup_) {
            lookup_ = std::make_shared<util::StericLookupNew>();
            //this->_setup_sterics(msg, lookup_);

            auto points = math::Points();
            for(auto const & b : scaffold_->beads()) {
                points.push_back(b.center());
            }
            lookup_->add_points(points);

            start_name_ = "A1003-A3926";
            end_name_ = "A1019-A3915";

            auto start_end_pos = msg->get_node(0)->data()->get_end_index(start_name_);
            auto end_end_pos = msg->get_node(0)->data()->get_end_index(end_name_);

            start_ = NodeIndexandEdge{0, start_end_pos};
            end_ = NodeIndexandEdge{0, end_end_pos};

            target_an_aligned_end_ = false;
            if(end_end_pos == msg->get_node(0)->data()->block_end_add()) {
                target_an_aligned_end_ = true;
            }
        }

        has_setup_ = true;

        auto p = std::make_shared<SequenceOptProblem>(msg, start_, end_, start_name_, end_name_,
                                                      lookup_, target_an_aligned_end_);
        return p;
    }

private:
    bool has_setup_;
    util::StericLookupNewOP lookup_;
    NodeIndexandEdge start_, end_;
    String start_name_, end_name_;
    bool target_an_aligned_end_;
    motif::MotifOP scaffold_;
};


struct SequenceOptimizationParameters {
public:
    inline
    SequenceOptimizationParameters() {
        problem = "TTR";
        helices = "ideal_helices";
        rounds = 1;
        motifs = 3;
        min_helix_size = 6;
        max_helix_size = 18;
    }

public:
    String problem, helices;
    int rounds, motifs;
    int min_helix_size, max_helix_size;
};

class SequenceOptimizationBenchmarks : public base::Application {
public:
    SequenceOptimizationBenchmarks();

    ~SequenceOptimizationBenchmarks() {}

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

    std::vector<motif::MotifStateOPs>
    _get_libraries();

    motif_search::MotifStateMonteCarloOP
    _get_search(
            SequenceOptProblemOP,
            std::vector<motif::MotifStateOPs> const &);

    SequenceOptimizer3DOP
    _get_optimizer(
            SequenceOptProblemOP,
            motif_data_structure::MotifGraphOP);

    String const &
    _get_motif_names(
            motif_data_structure::MotifGraphOP);

private:
    SequenceOptimizationParameters parameters_;
    String motif_names_;

};


#endif //TEST_SEQUENCE_OPTIMIZATION_BENCHMARKS_H
