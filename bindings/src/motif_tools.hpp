#ifndef PYBIND11_MOTIF_TOOLS_HPP
#define PYBIND11_MOTIF_TOOLS_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// motif_tools includes
#include <motif_tools/segmenter.h>

namespace motif_tools {
    namespace py = pybind11;

    void
    add_bindings(py::module_ & m) {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// motif_tools
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // classes

        py::class_<Pair, std::shared_ptr<Pair>>(m, "Pair")
                // ctors
                .def(py::init<structure::ResidueOP const &,structure::ResidueOP const &,int>())
                        // methods
                .def("contains",[] (Pair const & ptr, structure::ResidueOP const & res) -> int  {
                    return ptr.contains(res); } )
                        // public attributes
                .def_readwrite("res1", &Pair::res1)
                .def_readwrite("res2", &Pair::res2)
                .def_readwrite("dist", &Pair::dist)
                ;

        py::class_<PairSearch, std::shared_ptr<PairSearch>>(m, "PairSearch")
                // ctors
                .def(py::init<>())
                        // methods
                .def("search",[] (PairSearch  & ptr, structure::ResidueOPs const & res, PairOPs pairs, PairOPs end_pairs) -> PairSearchNodes const & {
                    return ptr.search(res, pairs, end_pairs); } )
                ;

        py::class_<PairSearchNode, std::shared_ptr<PairSearchNode>>(m, "PairSearchNode")
                // ctors
                .def(py::init<motif_tools::PairOPs>())
                        // methods
                .def("contains",[] (PairSearchNode const & ptr, structure::ResidueOP const & res) ->  int {
                    return ptr.contains(res); } )
                .def("contains_pair",[] (PairSearchNode  & ptr, PairOP pair) ->  int {
                    return ptr.contains_pair(pair); } )
                        // public attributes
                .def_readwrite("pairs", &PairSearchNode::pairs)
                .def_readwrite("score", &PairSearchNode::score)
                ;
/*
        py::class_<PairSearchNodeCompare, std::shared_ptr<PairSearchNodeCompare>>(m, "PairSearchNodeCompare")
		// operators
		.def(py::self () PairSearchNode const &)
		;
*/
        py::class_<Segmenter, std::shared_ptr<Segmenter>>(m, "Segmenter")
                // ctors
                .def(py::init<>())
                        // methods
                .def("apply",[] (Segmenter  & ptr, structure::RNAStructureOP const & m,structure::BasepairOPs const & bps) -> SegmentsOP {
                    return ptr.apply(m, bps); } )
                ;

        py::class_<Segments, std::shared_ptr<Segments>>(m, "Segments")
                // ctors
                .def(py::init<motif::MotifOP const &,motif::MotifOP const &>())
                        // public attributes
                .def_readwrite("removed", &Segments::removed)
                .def_readwrite("remaining", &Segments::remaining)
                ;


    }

}


#endif // PYBIND11_MOTIF_TOOLS_HPP