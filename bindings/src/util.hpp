#ifndef PYBIND11_UTIL_HPP
#define PYBIND11_UTIL_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// util includes
#include <util/basic_io.hpp>
#include <util/cartesian_product.h>
#include <util/monte_carlo.h>
#include <util/motif_type.h>
#include <util/random_number_generator.h>
#include <util/sqlite3_connection.h>
#include <util/steric_lookup.hpp>
#include <util/uuid.h>
#include <util/x3dna.h>

namespace util {
    namespace py = pybind11;
    void
    add_bindings(py::module_ & m) {

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// util
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Exceptions
        py::register_exception<Sqlite3ConnectionException>(m, "Sqlite3ConnectionException");
        py::register_exception<X3dnaException>(m, "X3dnaException");


        // free functions

        m.def("compare_bps", [] (X3dna::X3Basepairs & lhs, X3dna::X3Basepairs & rhs) -> String {
                  return compare_bps(lhs, rhs); },
              py::arg("lhs"),
              py::arg("rhs")
        );

        m.def("get_str_from_x3dna_type", [] (X3dnaBPType type) -> String {
                  return get_str_from_x3dna_type(type); },
              py::arg("type")
        );

        m.def("get_x3dna_by_type", [] (String const & name) -> X3dnaBPType {
                  return get_x3dna_by_type(name); },
              py::arg("name")
        );

        m.def("json_cleanup", [] () {
            json_cleanup(); }
        );

        m.def("points_to_pdb", [] (String const & filename, math::Points const & points) {
                  points_to_pdb(filename, points); },
              py::arg("filename"),
              py::arg("points")
        );

        m.def("points_to_pdb_str", [] (math::Points const & points) -> String {
                  return points_to_pdb_str(points); },
              py::arg("points")
        );

        m.def("str_to_type", [] (String const s) -> MotifType const {
                  return str_to_type(s); },
              py::arg("s")
        );

        m.def("type_to_str", [] (MotifType const mtype) -> String const {
                  return type_to_str(mtype); },
              py::arg("mtype")
        );
        // classes
#define Value double
#define Values std::vector<Value>
        py::class_<CartesianProduct<Value>, std::shared_ptr<CartesianProduct<Value>>>(m, "CartesianProductFloat")
                // ctors
                .def(py::init<>())
                .def(py::init<std::vector<Values> const &>(), py::arg("values"))
                        // methods
                .def("setup",[] (CartesianProduct<Value>  & ptr, std::vector<Values> const & values) {
                    ptr.setup(values); }, py::arg("values") )
                .def("end",[] (CartesianProduct<Value>  & ptr) ->  int const {
                    return ptr.end(); } )
                .def("next",[] (CartesianProduct<Value>  & ptr) -> Values const & {
                    return ptr.next(); } )
                ;
#undef Values
#undef Value
        py::class_<MonteCarlo, std::shared_ptr<MonteCarlo>>(m, "MonteCarlo")
                // ctors
                .def(py::init<float>(), py::arg("temperature") = 1.0f)
                        // methods
                .def("accept",[] (MonteCarlo  & ptr, float current, float next) ->  int {
                    return ptr.accept(current, next); }, py::arg("current"), py::arg("next") )
                .def("set_temperature",[] (MonteCarlo  & ptr, float new_temp) {
                    ptr.set_temperature(new_temp); }, py::arg("new_temp") )
                .def("get_temperature",[] (MonteCarlo  & ptr) ->  float {
                    return ptr.get_temperature(); } )
                .def("scale_temperature",[] (MonteCarlo  & ptr, float scale) {
                    ptr.scale_temperature(scale); }, py::arg("scale") )
                ;

        py::class_<RandomNumberGenerator, std::shared_ptr<RandomNumberGenerator>>(m, "RandomNumberGenerator")
                // ctors
                .def(py::init<>())
                        // methods
                .def("rand",[] (RandomNumberGenerator  & ptr) ->  float {
                    return ptr.rand(); } )
                .def("randrange",[] (RandomNumberGenerator  & ptr, int i) ->  int {
                    return ptr.randrange(i); }, py::arg("i") )
                ;

        py::class_<Sqlite3Connection, std::shared_ptr<Sqlite3Connection>>(m, "Sqlite3Connection")
                // ctors
                .def(py::init<>())
                .def(py::init<String const>(), py::arg("name"))
                        // methods
                .def("query",[] (Sqlite3Connection  & ptr, String const & query_statement) {
                    ptr.query(query_statement); }, py::arg("query_statement") )
                .def("count",[] (Sqlite3Connection  & ptr) -> int {
                    return ptr.count(); } )
                .def("fetch_one",[] (Sqlite3Connection  & ptr, String const & query_statement) -> Strings {
                    return ptr.fetch_one(query_statement); }, py::arg("query_statement") )
                ;

        py::class_<StericLookup, std::shared_ptr<StericLookup>>(m, "StericLookup")
                // ctors
                .def(py::init<>())
                .def(py::init<float,float,int>(), py::arg("grid_size"), py::arg("cutoff"), py::arg("radius"))
                        // methods
                .def("add_point",[] (StericLookup  & ptr, math::Point const & p) {
                    ptr.add_point(p); }, py::arg("p") )
                .def("add_points",[] (StericLookup  & ptr, math::Points const & points) {
                    ptr.add_points(points); }, py::arg("points") )
                .def("clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
                    return ptr.clash(p); }, py::arg("p") )
                .def("clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
                    return ptr.clash(p); } , py::arg("p"))
                .def("better_clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
                    return ptr.better_clash(p); } , py::arg("p"))
                .def("total_clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
                    return ptr.total_clash(p); } , py::arg("p"))
                .def("total_clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
                    return ptr.total_clash(p); } , py::arg("p"))
                ;

        py::class_<StericLookupNew, std::shared_ptr<StericLookupNew>>(m, "StericLookupNew")
                // ctors
                .def(py::init<>())
                        // methods
                .def("add_point",[] (StericLookupNew  & ptr, math::Point const & p) {
                    ptr.add_point(p); } , py::arg("p"))
                .def("add_points",[] (StericLookupNew  & ptr, math::Points const & points) {
                    ptr.add_points(points); } , py::arg("points"))
                .def("clash",[] (StericLookupNew  & ptr, math::Point const & p) -> bool {
                    return ptr.clash(p); } , py::arg("p"))
                .def("clash",[] (StericLookupNew  & ptr, math::Point const & p) -> bool {
                    return ptr.clash(p); } , py::arg("p"))
                .def("to_pdb",[] (StericLookupNew  & ptr, String const & pdb_name) {
                    ptr.to_pdb(pdb_name); } , py::arg("pdb_name"))
                .def("size",[] (StericLookupNew  & ptr) -> int {
                    return ptr.size(); } )
                ;

        py::class_<Uuid, std::shared_ptr<Uuid>>(m, "Uuid")
                // ctors
                .def(py::init<>())
                        // methods
                .def("s_uuid",[] (Uuid const & ptr) ->  String const & {
                    return ptr.s_uuid(); } )
                        // operators
                .def(py::self == py::self)
                .def(py::self == py::self)
                .def(py::self != py::self)
                .def(py::self != py::self)
                .def(py::self < py::self)
                ;

        py::class_<UuidCompare, std::shared_ptr<UuidCompare>>(m, "UuidCompare")
	        .def(py::init<>())
            .def("compare", [] (UuidCompare & comp,  Uuid const & uuid1 , Uuid const & uuid2 ) {
                return comp(uuid1, uuid2); }, py::arg("uuid1"), py::arg("uuid2"))
        // operators
//		.def(py::self () Uuid const &)
		;
        py::class_<X3dna::X3BPInfo, std::shared_ptr<X3dna::X3BPInfo>>(m, "X3BPInfo")
                // public attributes
                .def_readwrite("res1_type_name", &X3dna::X3BPInfo::res1_type_name)
                .def_readwrite("res2_type_name", &X3dna::X3BPInfo::res2_type_name)
                .def_readwrite("res1_name", &X3dna::X3BPInfo::res1_name)
                .def_readwrite("res2_name", &X3dna::X3BPInfo::res2_name)
                .def_readwrite("res1_chain_id", &X3dna::X3BPInfo::res1_chain_id)
                .def_readwrite("res2_chain_id", &X3dna::X3BPInfo::res2_chain_id)
                .def_readwrite("res1_num", &X3dna::X3BPInfo::res1_num)
                .def_readwrite("res2_num", &X3dna::X3BPInfo::res2_num)
                ;

        py::class_<X3dna::X3Basepair, std::shared_ptr<X3dna::X3Basepair>>(m, "X3Basepair")
                .def(py::init<>())
                // methods
                .def("valid",[] (X3dna::X3Basepair const & ptr) -> bool {
                    return ptr.valid(); } )
                .def("to_string",[] (X3dna::X3Basepair const & ptr) -> String {
                    return ptr.to_string(); } )
                .def("key",[] (X3dna::X3Basepair const & ptr) -> unsigned int {
                    return ptr.key(); } )
                        // public attributes
                .def_readwrite("res1", &X3dna::X3Basepair::res1)
                .def_readwrite("res2", &X3dna::X3Basepair::res2)
                .def_readwrite("d", &X3dna::X3Basepair::d)
                .def_readwrite("r", &X3dna::X3Basepair::r)
                .def_readwrite("bp_type", &X3dna::X3Basepair::bp_type)
                ;

        py::class_<X3dna::X3Motif, std::shared_ptr<X3dna::X3Motif>>(m, "X3Motif")
                // public attributes
                .def(py::init<>())
                .def_readwrite("residues", &X3dna::X3Motif::residues)
                .def_readwrite("mtype", &X3dna::X3Motif::mtype)
                ;

        py::class_<X3dna::X3Residue, std::shared_ptr<X3dna::X3Residue>>(m, "X3Residue")
                // methods
                .def(py::init<>())
                .def(py::init<int, char, char>(), py::arg("nnum") = -1, py::arg("nchain_id") = ' ', py::arg("ni_code") = ' ')
                .def("valid",[] (X3dna::X3Residue const & ptr) -> bool {
                    return ptr.valid(); } )
                        // operators
                .def(py::self == py::self)
                        // public attributes
                .def_readwrite("num", &X3dna::X3Residue::num)
                .def_readwrite("chain_id", &X3dna::X3Residue::chain_id)
                .def_readwrite("i_code", &X3dna::X3Residue::i_code)
                ;

        py::class_<X3dna, std::shared_ptr<X3dna>>(m, "X3dna")
                // ctors
                .def(py::init<>())
                        // methods
                .def("set_envs",[] (X3dna  & ptr) {
                    ptr.set_envs(); } )
                .def("get_basepairs",[] (X3dna const & ptr, String const & pdb_path) -> X3dna::X3Basepairs {
                    return ptr.get_basepairs(pdb_path); }, py::arg("pdb_path") )
                .def("get_basepairs_json",[] (X3dna const & ptr, String const & pdb_path) -> X3dna::X3Basepairs  {
                    return ptr.get_basepairs_json(pdb_path); }, py::arg("pdb_path") )
                .def("get_motifs",[] (X3dna const & ptr, String const & strings) -> X3dna::X3Motifs {
                    return ptr.get_motifs(strings); }, py::arg("strings") )
                .def("get_motifs",[] (X3dna const & ptr, String const & str) -> X3dna::X3Motifs {
                    return ptr.get_motifs(str); }, py::arg("str") )
                .def("set_rebuild_files",[] (X3dna const & ptr, bool rebuild_files) {
                    ptr.set_rebuild_files(rebuild_files); }, py::arg("rebuild_files")  )
                ;

        // enums
        py::enum_<MotifType>(m,"MotifType")
                .value("TWOWAY",MotifType::TWOWAY)
                .value("NWAY",MotifType::NWAY)
                .value("HAIRPIN",MotifType::HAIRPIN)
                .value("TCONTACT_HP_HP",MotifType::TCONTACT_HP_HP)
                .value("TCONTACT_H_HP",MotifType::TCONTACT_H_HP)
                .value("TCONTACT_H_H",MotifType::TCONTACT_H_H)
                .value("T_T",MotifType::T_T)
                .value("T_T_T",MotifType::T_T_T)
                .value("TWOWAY_SEGMENTS",MotifType::TWOWAY_SEGMENTS)
                .value("HELIX",MotifType::HELIX)
                .value("SSTRAND",MotifType::SSTRAND)
                .value("TCONTACT",MotifType::TCONTACT)
                .value("UNKNOWN",MotifType::UNKNOWN)
                .value("ALL",MotifType::ALL)
                ;
        py::enum_<X3dnaBPType>(m,"X3dnaBPType")
                .value("cmU",X3dnaBPType::cmU)
                .value("cMUM",X3dnaBPType::cMUM)
                .value("tWPW",X3dnaBPType::tWPW)
                .value("cDPM",X3dnaBPType::cDPM)
                .value("DWPW",X3dnaBPType::DWPW)
                .value("tWUM",X3dnaBPType::tWUM)
                .value("tmUM",X3dnaBPType::tmUM)
                .value("cWPM",X3dnaBPType::cWPM)
                .value("DWUW",X3dnaBPType::DWUW)
                .value("cMPD",X3dnaBPType::cMPD)
                .value("cDUm",X3dnaBPType::cDUm)
                .value("cMPW",X3dnaBPType::cMPW)
                .value("tMPm",X3dnaBPType::tMPm)
                .value("tMUW",X3dnaBPType::tMUW)
                .value("cmUm",X3dnaBPType::cmUm)
                .value("cMUW",X3dnaBPType::cMUW)
                .value("cWUW",X3dnaBPType::cWUW)
                .value("cDUM",X3dnaBPType::cDUM)
                .value("cmPM",X3dnaBPType::cmPM)
                .value("cmUM",X3dnaBPType::cmUM)
                .value("DDDD",X3dnaBPType::DDDD)
                .value("cmUW",X3dnaBPType::cmUW)
                .value("tMUm",X3dnaBPType::tMUm)
                .value("cDUW",X3dnaBPType::cDUW)
                .value("cMPm",X3dnaBPType::cMPm)
                .value("cMUm",X3dnaBPType::cMUm)
                .value("cDDD",X3dnaBPType::cDDD)
                .value("tWPm",X3dnaBPType::tWPm)
                .value("cDPm",X3dnaBPType::cDPm)
                .value("tmPm",X3dnaBPType::tmPm)
                .value("tWPD",X3dnaBPType::tWPD)
                .value("tmPW",X3dnaBPType::tmPW)
                .value("tDDD",X3dnaBPType::tDDD)
                .value("cWUD",X3dnaBPType::cWUD)
                .value("cWUM",X3dnaBPType::cWUM)
                .value("tDUW",X3dnaBPType::tDUW)
                .value("tMPM",X3dnaBPType::tMPM)
                .value("tDUM",X3dnaBPType::tDUM)
                .value("cMUD",X3dnaBPType::cMUD)
                .value("cWUm",X3dnaBPType::cWUm)
                .value("tDPm",X3dnaBPType::tDPm)
                .value("tMUD",X3dnaBPType::tMUD)
                .value("cmPW",X3dnaBPType::cmPW)
                .value("cMPM",X3dnaBPType::cMPM)
                .value("cmPD",X3dnaBPType::cmPD)
                .value("cmUD",X3dnaBPType::cmUD)
                .value("cDUD",X3dnaBPType::cDUD)
                .value("cWPW",X3dnaBPType::cWPW)
                .value("tDUD",X3dnaBPType::tDUD)
                .value("tDPW",X3dnaBPType::tDPW)
                .value("tmUm",X3dnaBPType::tmUm)
                .value("cWPD",X3dnaBPType::cWPD)
                .value("tmPD",X3dnaBPType::tmPD)
                .value("tDPD",X3dnaBPType::tDPD)
                .value("cDPD",X3dnaBPType::cDPD)
                .value("tDUm",X3dnaBPType::tDUm)
                .value("tDPM",X3dnaBPType::tDPM)
                .value("tWUD",X3dnaBPType::tWUD)
                .value("tmUW",X3dnaBPType::tmUW)
                .value("tMUM",X3dnaBPType::tMUM)
                .value("tMPD",X3dnaBPType::tMPD)
                .value("cDPW",X3dnaBPType::cDPW)
                .value("tmPM",X3dnaBPType::tmPM)
                .value("tWUm",X3dnaBPType::tWUm)
                .value("cWPm",X3dnaBPType::cWPm)
                .value("tmUD",X3dnaBPType::tmUD)
                .value("tWPM",X3dnaBPType::tWPM)
                .value("DWPm",X3dnaBPType::DWPm)
                .value("tMPW",X3dnaBPType::tMPW)
                .value("DDPm",X3dnaBPType::DDPm)
                .value("tWUW",X3dnaBPType::tWUW)
                .value("cmPm",X3dnaBPType::cmPm)
                .value("DWUm",X3dnaBPType::DWUm)
                .value("DMPm",X3dnaBPType::DMPm)
                .value("DWPM",X3dnaBPType::DWPM)
                .value("DMPM",X3dnaBPType::DMPM)
                .value("DmPW",X3dnaBPType::DmPW)
                .value("DWUM",X3dnaBPType::DWUM)
                .value("DmPm",X3dnaBPType::DmPm)
                .value("DDUM",X3dnaBPType::DDUM)
                .value("DMUm",X3dnaBPType::DMUm)
                .value("DDUm",X3dnaBPType::DDUm)
                .value("DMPW",X3dnaBPType::DMPW)
                .value("DMPD",X3dnaBPType::DMPD)
                .value("DMUM",X3dnaBPType::DMUM)
                .value("DmUm",X3dnaBPType::DmUm)
                .value("DMUW",X3dnaBPType::DMUW)
                .value("DWUD",X3dnaBPType::DWUD)
                ;
    }

}

#endif // PYBIND11_UTIL_HPP
