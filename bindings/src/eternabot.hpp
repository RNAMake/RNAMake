#ifndef PYBIND11_ETERNABOT_HPP
#define PYBIND11_ETERNABOT_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>

#include <memory>
// eternabot includes
#include <eternabot/feature_generator.h>
#include <eternabot/scorer.h>
#include <eternabot/sequence_designer.h>
#include <eternabot/strategy.h>
#include <eternabot/strategy/a_basic_test.h>
#include <eternabot/strategy/berex_test.h>
#include <eternabot/strategy/clear_plot.h>
#include <eternabot/strategy/direction_of_gc.h>
#include <eternabot/strategy/num_of_yellow.h>

namespace eternabot {
    namespace py = pybind11;
    void
    add_bindings(py::module_& m) {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// eternabot
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // classes

        py::class_<ABasicTest, std::shared_ptr<ABasicTest>>(m, "ABasicTest")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (ABasicTest  & ptr, FeaturesOP const & features) ->  float {
                    return ptr.score(features); } )
                        // inherited methods

                .def("mean",[] (Strategy const & ptr) ->  float {
                    return ptr.mean(); } )
                .def("stdev",[] (Strategy const & ptr) ->  float {
                    return ptr.stdev(); } )
                ;

        py::class_<BerexTest, std::shared_ptr<BerexTest>>(m, "BerexTest")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (BerexTest  & ptr, FeaturesOP const & features) -> float {
                    return ptr.score(features); } )
                        // inherited methods
                .def("mean",[] (Strategy const & ptr) ->  float {
                    return ptr.mean(); } )
                .def("stdev",[] (Strategy const & ptr) ->  float {
                    return ptr.stdev(); } )
                ;

        py::class_<CleanPlotStackCapsandSafeGC, std::shared_ptr<CleanPlotStackCapsandSafeGC>>(m, "CleanPlotStackCapsandSafeGC")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (CleanPlotStackCapsandSafeGC  & ptr, FeaturesOP const & features) -> float {
                    return ptr.score(features); } )
                        // inherited methods
                .def("mean",[] (Strategy const & ptr) ->  float {
                    return ptr.mean(); } )
                .def("stdev",[] (Strategy const & ptr) ->  float {
                    return ptr.stdev(); } )
                ;

        py::class_<DirectionofGCPairsinMultiLoops, std::shared_ptr<DirectionofGCPairsinMultiLoops>>(m, "DirectionofGCPairsinMultiLoops")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (DirectionofGCPairsinMultiLoops  & ptr, FeaturesOP const & features) -> float {
                    return ptr.score(features); } )
                        // inherited methods
                .def("mean",[] (Strategy const & ptr) ->  float {
                    return ptr.mean(); } )
                .def("stdev",[] (Strategy const & ptr) ->  float {
                    return ptr.stdev(); } )
                ;

        py::class_<FeatureGenerator, std::shared_ptr<FeatureGenerator>>(m, "FeatureGenerator")
                // ctors
                .def(py::init<>())
                        // methods
                .def("get_features",[] (FeatureGenerator  & ptr, secondary_structure::PoseOP const & p) -> FeaturesOP {
                    return ptr.get_features(p); } )
                .def("update_features",[] (FeatureGenerator  & ptr, FeaturesOP & features, secondary_structure::PoseOP const & p) {
                    ptr.update_features(features, p); } )
                ;

        py::class_<Features, std::shared_ptr<Features>>(m, "Features")
                // ctors
                .def(py::init<>())
                        // public attributes
                .def_readwrite("length", &Features::length)
                .def_readwrite("a_count", &Features::a_count)
                .def_readwrite("c_count", &Features::c_count)
                .def_readwrite("g_count", &Features::g_count)
                .def_readwrite("u_count", &Features::u_count)
                .def_readwrite("gu", &Features::gu)
                .def_readwrite("gc", &Features::gc)
                .def_readwrite("ua", &Features::ua)
                .def_readwrite("meltpoint", &Features::meltpoint)
                .def_readwrite("fe", &Features::fe)
                .def_readwrite("dotplot", &Features::dotplot)
                .def_readwrite("pairmap", &Features::pairmap)
                .def_readwrite("helices", &Features::helices)
                ;

        py::class_<NumofYellowNucleotidesperLengthofString, std::shared_ptr<NumofYellowNucleotidesperLengthofString>>(m, "NumofYellowNucleotidesperLengthofString")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (NumofYellowNucleotidesperLengthofString  & ptr, FeaturesOP const & features) -> float {
                    return ptr.score(features); } )
                        // inherited methods
                .def("mean",[] (Strategy const & ptr) ->  float {
                    return ptr.mean(); } )
                .def("stdev",[] (Strategy const & ptr) ->  float {
                    return ptr.stdev(); } )
                ;
/*
        py::class_<Parameters, std::shared_ptr<Parameters>>(m, "Parameters")
		;
*/
        py::class_<eternabot::Scorer, std::shared_ptr<eternabot::Scorer>>(m, "EternabotScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("setup",[] (eternabot::Scorer  & ptr, secondary_structure::PoseOP const & p) {
                    ptr.setup(p); } )
                .def("score_secondary_structure",[] (eternabot::Scorer  & ptr, secondary_structure::PoseOP const & p) -> float {
                    return ptr.score_secondary_structure(p); } )
                .def("scores",[] (eternabot::Scorer  & ptr) -> Floats const & {
                    return ptr.scores(); } )
                ;

        py::class_<SequenceDesigner, std::shared_ptr<SequenceDesigner>>(m, "SequenceDesigner")
                // ctors
                .def(py::init<>())
                        // methods
                .def("setup",[] (SequenceDesigner  & ptr) {
                    ptr.setup(); } )
                .def("design",[] (SequenceDesigner  & ptr, secondary_structure::PoseOP const & p) -> SequenceDesignerResultOPs const &  {
                    return ptr.design(p); } )
                .def("get_int_option",[] (SequenceDesigner  & ptr, String const & name) ->  float {
                    return ptr.get_int_option(name); } )
                .def("get_float_option",[] (SequenceDesigner  & ptr, String const & name) ->  float {
                    return ptr.get_float_option(name); } )
                .def("get_string_option",[] (SequenceDesigner  & ptr, String const & name) ->  String {
                    return ptr.get_string_option(name); } )
                .def("get_bool_option",[] (SequenceDesigner  & ptr, String const & name) ->  bool {
                    return ptr.get_bool_option(name); } )
                .def("has_option",[] (SequenceDesigner  & ptr, String const & name) ->  bool {
                    return ptr.has_option(name); } )
                .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
                ;

        py::class_<SequenceDesignerResult, std::shared_ptr<SequenceDesignerResult>>(m, "SequenceDesignerResult")
                // ctors
                .def(py::init<String const &,float>())
        // public attributes
        .def_readwrite("sequence", &SequenceDesignerResult::sequence)
                .def_readwrite("score", &SequenceDesignerResult::score)
                ;
/*
        py::class_<Strategy, std::shared_ptr<Strategy>>(m, "Strategy")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (Strategy  & ptr, FeaturesOP const &) -> float {
		 return ptr.score(&); } )
		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;
        py::class_<sequence_designer_result_less_than_key, std::shared_ptr<sequence_designer_result_less_than_key>>(m, "sequence_designer_result_less_than_key")
		// operators
		.def(py::self () eternabot::SequenceDesignerResultOP)
		;

*/
    }
}


#endif // PYBIND11_ETERNABOT_HPP