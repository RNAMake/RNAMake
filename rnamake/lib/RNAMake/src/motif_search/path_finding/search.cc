//
// Created by Joseph Yesselman on 3/17/19.
//

#include <base/log.h>
#include <motif_search/path_finding/search.h>
#include <base/global_constants.h>

namespace motif_search {
namespace path_finding {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// option functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Search::setup_options() {
    options_.add_option("sterics", true, base::OptionType::BOOL);
    options_.add_option("max_node_level", 100, base::OptionType::INT);
    options_.add_option("min_node_level", -1, base::OptionType::INT);
    options_.add_option("min_size", 0, base::OptionType::INT);
    options_.add_option("max_size", 1000000, base::OptionType::INT);
    options_.add_option("max_solutions", 1, base::OptionType::INT);
    options_.add_option("accept_score", 10, base::OptionType::FLOAT);
    options_.add_option("min_ss_score", 10000, base::OptionType::FLOAT);
    options_.add_option("max_steps", 1000000000, base::OptionType::FLOAT);
    options_.add_option("return_best", false, base::OptionType::BOOL);
    options_.add_option("helix_end", false, base::OptionType::BOOL);
    options_.lock_option_adding();
}

void
Search::update_var_options() {
    parameters_.sterics = options_.get_bool("sterics");
    parameters_.min_size = options_.get_int("min_size");
    parameters_.max_size = options_.get_int("max_size");
    parameters_.max_solutions = options_.get_int("max_solutions");
    parameters_.max_node_level = options_.get_int("max_node_level");
    parameters_.min_node_level = options_.get_int("min_node_level");
    parameters_.accept_score = options_.get_float("accept_score");
    parameters_.min_ss_score = options_.get_float("min_ss_score");
    parameters_.max_steps = options_.get_float("max_steps");
    parameters_.helix_end = options_.get_bool("helix_end");
    parameters_.return_best = options_.get_bool("return_best");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// search functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
Search::setup(
        motif_search::ProblemOP p)  {

    auto ms = std::make_shared<motif::MotifState>(
            "start", Strings{"start", "start"}, Strings{"", ""},
            structure::BasepairStateOPs {p->start, p->start},
            math::Points(), 0, 0, 0);

    auto start_n = std::make_shared<Node>(ms, nullptr, 100000, 1, 0, -1);
    queue_ = NodeQueue();
    queue_.push(start_n);

    scorer_->set_target(p->end, p->target_an_aligned_end);

    if(p->lookup->size() != 0) {
        lookup_ = p->lookup;
        using_lookup_ = true;
    }
    else {
        using_lookup_ = false;
    }
}

motif_search::SolutionOP
Search::next() {

    auto current = NodeOP(nullptr);
    auto child = NodeOP(nullptr);
    auto best_node = NodeOP(nullptr);
    auto selector_data= SelectorNodeDataOP(nullptr);
    int steps = 0;
    float score = 0, best = 10000, new_score = 0;
    while (! queue_.empty()) {
        current = queue_.top();
        queue_.pop();

        steps += 1;

        if(steps > parameters_.max_steps) { LOG_VERBOSE << "reached max steps"; break; }

        score = scorer_->accept_score(*current);
        if(score < best) {
            best = score;
            best_node = current;
            LOG_VERBOSE << "best_score=" << best << " motifs_in_solution=" << current->level() << " steps=" << steps;
        }

        // accept solution
        if(score < parameters_.accept_score && _accept_node(*current)) {
            _get_solution_motif_names(current, motif_names_);
            if(! solution_filter_->accept(motif_names_)) { continue; }

            LOG_VERBOSE << "solution found!";
            auto msg = _graph_from_node(current);
            return std::make_shared<Solution>(msg, score);
        }

        if(current->level()+1 >= parameters_.max_node_level) { continue; }

        int j = 0;
        selector_->start(current->node_type());
        while(! selector_->finished()) {
            selector_data = selector_->next();
            j = -1;
            for(auto const & end : current->state()->end_states()) {
                j++;
                if(j == 0) { continue; }

                for (auto & ms : selector_data->motif_states) {
                    if(current->size() + ms->size() > parameters_.max_size) { continue; }

                    aligner_.get_aligned_motif_state(end, ms);
                    new_score = scorer_->score(*ms, *current);

                    if(new_score > current->score()) { continue; }
                    if(parameters_.sterics && _steric_clash(*ms, *current)) { continue; }

                    child = std::make_shared<Node>(std::make_shared<motif::MotifState>(*ms),
                                                   current, new_score, current->level() + 1,
                                                   j, selector_data->type);
                    queue_.push(child);
                }
            }
        }
    }

    if(parameters_.return_best) {
        LOG_VERBOSE << "no solution met constraints, returning best";
        auto msg = _graph_from_node(best_node);
        return std::make_shared<Solution>(msg, best);
    }

    LOG_VERBOSE << "search ran out of options";
    return SolutionOP(nullptr);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

motif_data_structure::MotifStateGraphOP
Search::_graph_from_node(
        NodeOP node) {
    auto current = node;
    auto nodes = NodeOPs();
    while(current != nullptr) {
        nodes.push_back(current);
        current = current->parent();
    }
    nodes.pop_back();
    std::reverse(nodes.begin(), nodes.end());
    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);
    int i = 0, j = 0;
    for(auto const & n : nodes) {
        auto ms = std::make_shared<motif::MotifState>(*n->state());
        ms->new_uuids();
        if(i == 0) {
            msg->add_state(ms);
        }
        else {
            j = msg->add_state(ms, -1, n->parent_end_index());
            if(j == -1) {
                LOG_ERROR << "something went horribly wrong, cannot build solution";
            }
        }
    }
    return msg;
}

bool
Search::_steric_clash(
        motif::MotifState const & ms,
        Node const & n) {
    auto clash = false;
    if(using_lookup_) {
        clash = lookup_->clash(ms.beads());
        if(clash) { return true; }
    }

    auto current = n.parent();
    if(current == nullptr) { return false;}
    float dist = 0;
    // first node will be start node, so no reason to compute the start node as its a dummy
    // might need to add a check to see if state has beads or this will seg fault ... -- JDY
    while(current->parent() != nullptr) {
        // quick check to see how far these two motifs are from each other. If too far sterics calc is probably
        // not worth, need to benchamark -- JDY
        dist = ms.beads()[0].distance(current->state()->beads()[0]);
        if(dist > 25) {
            current = current->parent();
            continue;
        }
        for(auto const & b1 : ms.beads()) {
            for(auto const & b2 : current->state()->beads()) {
                dist = b1.distance(b2);
                if(dist < STERIC_CLASH_RADIUS) {
                    return true;
                }
            }
        }
        current = current->parent();

    }
    return false;
}

void
Search::_get_solution_motif_names(
        NodeOP n,
        Strings & motif_names) {
    motif_names_.resize(0);
    motif_names.push_back(n->state()->name());
    while(n->parent() != nullptr) {
        n = n->parent();
        motif_names_.push_back(n->state()->name());
    }
}


}
}

















