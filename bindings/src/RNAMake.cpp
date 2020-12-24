// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// binding_includes
#include <base.hpp>
#include <math.hpp>
#include <data_structure.hpp>
#include <util.hpp>
#include <secondary_structure.hpp>
#include <eternabot.hpp>
#include <structure.hpp>
#include <motif.hpp>
#include <motif_tools.hpp>
#include <resources.hpp>
#include <motif_data_structures.hpp>
#include <thermo_fluctuation.hpp>
#include <motif_search.hpp>

namespace py = pybind11;

PYBIND11_MODULE(RNAMake, m) {

    auto base = m.def_submodule("base");
    base::add_bindings(base);

    auto math = m.def_submodule("math");
    math::add_bindings(math);

    auto data_structure = m.def_submodule("data_structure");
    data_structure::add_bindings(data_structure);

    auto util = m.def_submodule("util");
    util::add_bindings(util);

    auto secondary_structure = m.def_submodule("secondary_structure");
    secondary_structure::add_bindings(secondary_structure);

    auto eternabot = m.def_submodule("eternabot");
    eternabot::add_bindings(eternabot);

    auto structure = m.def_submodule("structure");
    structure::add_bindings(structure);

    auto motif = m.def_submodule("motif");
    motif::add_bindings(motif);

    auto motif_tools = m.def_submodule("motif_tools");
    motif_tools::add_bindings(motif_tools);

    auto resources = m.def_submodule("resources");
    resources::add_bindings(resources);

    auto motif_data_structure = m.def_submodule("motif_data_structure");
    motif_data_structure::add_bindings(motif_data_structure);

    auto thermo_fluctuation = m.def_submodule("thermo_fluctuation");
    thermo_fluctuation::add_bindings(thermo_fluctuation);

    auto motif_search = m.def_submodule("motif_search");
    motif_search::add_bindings(motif_search);
}