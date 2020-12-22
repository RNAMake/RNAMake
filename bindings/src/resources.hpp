#ifndef PYBIND11_RESOURCES_HPP
#define PYBIND11_RESOURCES_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// resources includes
#include <resources/added_motif_library.h>
#include <resources/motif_ensemble_sqlite_connection.h>
#include <resources/motif_sqlite_connection.h>
#include <resources/motif_sqlite_library.h>
#include <resources/motif_state_ensemble_sqlite_library.h>
#include <resources/motif_state_sqlite_library.h>
#include <resources/resource_manager.h>
#include <resources/sqlite_library.h>

namespace resources {
    namespace py = pybind11;

    void
    add_bindings(py::module_ & m) {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// resources
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // free functions

        m.def("build_sqlite_library", [] (String const & path, std::vector<Strings>const & data, Strings const & keys, String const & primary_key) {
                  build_sqlite_library(path, data, keys, primary_key); },
              py::arg("path"),
              py::arg("data"),
              py::arg("keys"),
              py::arg("primary_key")
        );

        m.def("get_motif_from_resource_manager", [] (String const & name, String const & end_id, String const & end_name) ->  motif::MotifOP {
                  return get_motif_from_resource_manager(name, end_id, end_name); },
              py::arg("name") = dummy_name,
              py::arg("end_id") = dummy_end_id,
              py::arg("end_name") = dummy_name
        );

        m.def("sqlite3_escape", [] (String & unescaped_string) {
                  sqlite3_escape(unescaped_string); },
              py::arg("unescaped_string")
        );
        // classes

        py::class_<AddedMotifLibrary, std::shared_ptr<AddedMotifLibrary>>(m, "AddedMotifLibrary")
                // ctors
                .def(py::init<>())
                        // methods
                .def("add_motif",[] (AddedMotifLibrary  & ptr, motif::MotifOP const & m) {
                    ptr.add_motif(m); } )
                .def("get",[] (AddedMotifLibrary  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifOP {
                    return ptr.get(name, end_id, end_name); } )
                .def("get_multi",[] (AddedMotifLibrary  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifOPs {
                    return ptr.get_multi(name, end_id, end_name); } )
                .def("contains",[] (AddedMotifLibrary  & ptr, String const & name, String const & end_id, String const & end_name) -> int {
                    return ptr.contains(name, end_id, end_name); } )
                ;

        py::class_<Manager, std::unique_ptr<Manager,py::nodelete>>(m, "Manager")
                .def(py::init([](){
                    return std::unique_ptr<Manager,py::nodelete>(&Manager::instance());
                }))
                        //.def(py::init([](){return Manager::instance();}))
                        //.def("__init__",[]() {return Manager::instance();})
                        //.def(py::init([](){return Manager::instance();}))
                        // methods
                        //.def(py::init(),py::return_value_policy::automatic_reference)
                .def("instance",[] (Manager  & ptr) -> Manager & {
                    return ptr.instance(); } )
                .def("bp_step",[] (Manager  & ptr, String const & end_id) -> motif::MotifOP {
                    return ptr.bp_step(end_id); } )
                .def("bp_step_state",[] (Manager  & ptr, String const & end_id) -> motif::MotifStateOP {
                    return ptr.bp_step_state(end_id); } )
                .def("motif",[] (Manager  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifOP {
                    return ptr.motif(name, end_id, end_name); } )
                .def("motif_state",[] (Manager  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifStateOP {
                    return ptr.motif_state(name, end_id, end_name); } )
                .def("motif_state_ensemble",[] (Manager  & ptr, String const & name) -> motif::MotifStateEnsembleOP {
                    return ptr.motif_state_ensemble(name); } )
                .def("get_structure",[] (Manager  & ptr, String const & path, String name, int force_num_chains) -> structure::RNAStructureOP {
                    return ptr.get_structure(path, name, force_num_chains); } )
                .def("add_motif",[] (Manager  & ptr, String const & path, String name, util::MotifType mtype) {
                    ptr.add_motif(path, name, mtype); } )
                .def("add_motif",[] (Manager  & ptr, motif::MotifOP const & m, String name) {
                    ptr.add_motif(m, name); } )
                .def("register_motif",[] (Manager  & ptr, motif::MotifOP const & m) {
                    ptr.register_motif(m); } )
                .def("register_extra_motif_ensembles",[] (Manager  & ptr, String const & f_name) {
                    ptr.register_extra_motif_ensembles(f_name); } )
                .def("has_supplied_motif_ensemble",[] (Manager  & ptr, String const & m_name, String const & end_name) -> int {
                    return ptr.has_supplied_motif_ensemble(m_name, end_name); } )
                .def("get_supplied_motif_ensemble",[] (Manager  & ptr, String const & m_name, String const & end_name) -> motif::MotifEnsembleOP const &  {
                    return ptr.get_supplied_motif_ensemble(m_name, end_name); } )
                .def("get_helix_name",[] (Manager  & ptr, String const & str) -> String {
                    return ptr.get_helix_name(str); } )
                .def("get_motif_from_state",[] (Manager  & ptr, motif::MotifStateOP ms) -> motif::MotifOP {
                    return ptr.get_motif_from_state(ms); } )
                ;

        py::class_<MotifEnsembleSqliteConnection, std::shared_ptr<MotifEnsembleSqliteConnection>>(m, "MotifEnsembleSqliteConnection")
                // ctors
                .def(py::init<>())
                .def(py::init<String const &>())
        // methods
        .def("next",[] (MotifEnsembleSqliteConnection  & ptr) -> MotifEnsembleSqliteDataOP const & {
            return ptr.next(); } )
                ;

        py::class_<MotifEnsembleSqliteData, std::shared_ptr<MotifEnsembleSqliteData>>(m, "MotifEnsembleSqliteData")
                // ctors
                .def(py::init<>())
                        // public attributes
                .def_readwrite("data", &MotifEnsembleSqliteData::data)
                .def_readwrite("name", &MotifEnsembleSqliteData::name)
                .def_readwrite("id", &MotifEnsembleSqliteData::id)
                ;

        py::class_<MotifSqliteConnection, std::shared_ptr<MotifSqliteConnection>>(m, "MotifSqliteConnection")
                // ctors
                .def(py::init<>())
                .def(py::init<String const &>())
        // methods
        .def("next",[] (MotifSqliteConnection  & ptr) -> MotifSqliteDataOP const & {
            return ptr.next(); } )
                .def("contains",[] (MotifSqliteConnection  & ptr) -> MotifSqliteDataOP const & {
                    return ptr.contains(); } )
                .def("clear",[] (MotifSqliteConnection  & ptr) {
                    ptr.clear(); } )
                ;

        py::class_<MotifSqliteData, std::shared_ptr<MotifSqliteData>>(m, "MotifSqliteData")
                // ctors
                .def(py::init<>())
                        // public attributes
                .def_readwrite("data", &MotifSqliteData::data)
                .def_readwrite("name", &MotifSqliteData::name)
                .def_readwrite("end_name", &MotifSqliteData::end_name)
                .def_readwrite("end_id", &MotifSqliteData::end_id)
                .def_readwrite("id", &MotifSqliteData::id)
                ;

        py::class_<MotifSqliteLibrary, std::shared_ptr<MotifSqliteLibrary>>(m, "MotifSqliteLibrary")
                // ctors
                .def(py::init<String const &>())
        .def(py::init<int,std::filesystem::path const &>())
        // methods
        .def("begin",[] (MotifSqliteLibrary  & ptr) -> MotifSqliteLibrary::iterator {
            return ptr.begin(); } )
                .def("end",[] (MotifSqliteLibrary  & ptr) -> MotifSqliteLibrary::iterator {
                    return ptr.end(); } )
                .def("get_libnames",[] (MotifSqliteLibrary  & ptr) -> StringStringMap {
                    return ptr.get_libnames(); } )
                .def("get",[] (MotifSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifOP {
                    return ptr.get(name, end_id, end_name, id); } )
                .def("get_multi",[] (MotifSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifOPs {
                    return ptr.get_multi(name, end_id, end_name, id); } )
                .def("contains",[] (MotifSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> int {
                    return ptr.contains(name, end_id, end_name, id); } )
                .def("get_random",[] (MotifSqliteLibrary  & ptr) -> motif::MotifOP {
                    return ptr.get_random(); } )
                .def("load_all",[] (MotifSqliteLibrary  & ptr, int limit) {
                    ptr.load_all(limit); } )
                ;

        py::class_<MotifStateEnsembleSqliteLibrary, std::shared_ptr<MotifStateEnsembleSqliteLibrary>>(m, "MotifStateEnsembleSqliteLibrary")
                // ctors
                .def(py::init<String const &>())
        // methods
        .def("begin",[] (MotifStateEnsembleSqliteLibrary  & ptr) -> MotifStateEnsembleSqliteLibrary::iterator {
            return ptr.begin(); } )
                .def("end",[] (MotifStateEnsembleSqliteLibrary  & ptr) -> MotifStateEnsembleSqliteLibrary::iterator {
                    return ptr.end(); } )
                .def("get_libnames",[] (MotifStateEnsembleSqliteLibrary  & ptr) ->  StringStringMap {
                    return ptr.get_libnames(); } )
                .def("get",[] (MotifStateEnsembleSqliteLibrary  & ptr, String const & name, String const & id) -> motif::MotifStateEnsembleOP {
                    return ptr.get(name, id); } )
                .def("get_multi",[] (MotifStateEnsembleSqliteLibrary  & ptr, String const & name, String const & id) -> motif::MotifStateEnsembleOPs {
                    return ptr.get_multi(name, id); } )
                .def("contains",[] (MotifStateEnsembleSqliteLibrary  & ptr, String const & name, String const & id) -> int {
                    return ptr.contains(name, id); } )
                .def("get_random",[] (MotifStateEnsembleSqliteLibrary  & ptr) -> motif::MotifStateEnsembleOP {
                    return ptr.get_random(); } )
                .def("load_all",[] (MotifStateEnsembleSqliteLibrary  & ptr, int limit) {
                    ptr.load_all(limit); } )
                ;

        py::class_<MotifStateSqliteLibrary, std::shared_ptr<MotifStateSqliteLibrary>>(m, "MotifStateSqliteLibrary")
                // ctors
                .def(py::init<String const &>())
        // methods
        .def("begin",[] (MotifStateSqliteLibrary  & ptr) -> MotifStateSqliteLibrary::iterator {
            return ptr.begin(); } )
                .def("end",[] (MotifStateSqliteLibrary  & ptr) -> MotifStateSqliteLibrary::iterator {
                    return ptr.end(); } )
                .def("get_libnames",[] (MotifStateSqliteLibrary  & ptr) -> StringStringMap {
                    return ptr.get_libnames(); } )
                .def("get",[] (MotifStateSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifStateOP {
                    return ptr.get(name, end_id, end_name, id); } )
                .def("get_multi",[] (MotifStateSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifStateOPs {
                    return ptr.get_multi(name, end_id, end_name, id); } )
                .def("contains",[] (MotifStateSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> int {
                    return ptr.contains(name, end_id, end_name, id); } )
                .def("get_random",[] (MotifStateSqliteLibrary  & ptr) -> motif::MotifStateOP {
                    return ptr.get_random(); } )
                .def("load_all",[] (MotifStateSqliteLibrary  & ptr, int limit) {
                    ptr.load_all(limit); } )
                .def("get_name",[] (MotifStateSqliteLibrary  & ptr) -> String const & {
                    return ptr.get_name(); } )
                ;

        py::class_<SqliteLibrary, std::shared_ptr<SqliteLibrary>>(m, "SqliteLibrary")
                // ctors
                .def(py::init<>())
                ;
/*
        py::class_<iterator, std::shared_ptr<iterator>>(m, "iterator")
		// ctors
		.def(py::init<std::map<String, motif::MotifStateOP>::iterator const &>())
		// operators
		.def(py::self == py::self)
		.def(py::self != py::self)
		;

        py::class_<iterator, std::shared_ptr<iterator>>(m, "iterator")
		// ctors
		.def(py::init<std::map<String, motif::MotifStateEnsembleOP>::iterator const &>())
		// operators
		.def(py::self == py::self)
		.def(py::self != py::self)
		;

        py::class_<iterator, std::shared_ptr<iterator>>(m, "iterator")
		// ctors
		.def(py::init<std::map<String, motif::MotifOP>::iterator const &>())
		// operators
		.def(py::self == py::self)
		.def(py::self != py::self)
		;
*/
    }
}

#endif // PYBIND11_RESOURCES_HPP