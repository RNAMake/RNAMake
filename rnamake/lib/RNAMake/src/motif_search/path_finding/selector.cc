//
// Created by Joseph Yesselman on 3/18/19.
//

#include <motif_search/path_finding/selector.h>

namespace motif_search {
namespace path_finding {


void
Selector::connect(
        String const & name_i,
        String const & name_j) {

    int i = -1, j = -1;
    for (auto const & n : graph_) {
        if (n->data()->name == name_i) { i = n->index(); }
        if (n->data()->name == name_j) { j = n->index(); }
    }

    if (i == -1 || j == -1) {
        throw SelectorException("could not connect nodes: " + name_i + " " + name_j);
    }
    graph_.connect(i, j);
}


SelectorOP
default_selector() {
    auto s = std::make_shared<RoundRobinSelector>();
    s->add("ideal_helices");
    s->add("twoway");
    return s;
}


}
}