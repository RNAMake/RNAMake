// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

#include <sstream>
// base includes
#include <base/application.hpp>
#include <base/backtrace.h>
#include <base/cl_option.h>
#include <base/command_line_parser.hpp>
#include <base/env_manager.h>
#include <base/exception.h>
#include <base/file_io.h>
#include <base/log.h>
#include <base/option.h>
#include <base/settings.h>
#include <base/string.h>
#include <base/sys_interface.h>
#include <base/vector_container.h>

// math includes
#include <math/euler.h>
#include <math/hashing.h>
#include <math/numerical.h>
#include <math/quaternion.h>
#include <math/stats.h>
#include <math/transform.h>
#include <math/xyz_matrix.h>
#include <math/xyz_vector.h>

using namespace math;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_math,m) {

	// free functions

    m.def("are_floats_equal", [] (double const a, double const b, double tol) -> int {
		 return are_floats_equal(a, b, tol); },
	py::arg("a"),
	py::arg("b"),
	py::arg("tol") = 0.001 
    );

    m.def("are_xyzMatrix_equal", [] (Matrix const & m, Matrix const & mc) -> int {
		 return are_xyzMatrix_equal(m, mc); },
	py::arg("m"),
	py::arg("mc") 
    );

    m.def("are_xyzVector_equal", [] (Vector const & vec, Vector const & correct_vec, float tol) -> int {
		 return are_xyzVector_equal(vec, correct_vec, tol); },
	py::arg("vec"),
	py::arg("correct_vec"),
	py::arg("tol") = 0.001 
    );

    m.def("are_xyzVectors_equal", [] (Vectors const & v, Vectors const & vc) -> int {
		 return are_xyzVectors_equal(v, vc); },
	py::arg("v"),
	py::arg("vc") 
    );

    m.def("avg_unsigned_diff", [] (std::vector<double> const & x, std::vector<double> const & y) -> double {
		 return avg_unsigned_diff(x, y); },
	py::arg("x"),
	py::arg("y") 
    );

    m.def("axis_angle_from_matrix", [] (Matrix & m, AxisAngle & aa) {
		axis_angle_from_matrix(m, aa); },
	py::arg("m"),
	py::arg("aa") 
    );

    m.def("calc_euler", [] (Matrix & M, Vector & euler) {
		calc_euler(M, euler); },
	py::arg("M"),
	py::arg("euler") 
    );

    m.def("degrees", [] (float radians) -> float {
		 return degrees(radians); },
	py::arg("radians") 
    );

    m.def("dot_vector", [] (math::Matrix m, Vector const & v, Vector & vr) {
		dot_vector(m, v, vr); },
	py::arg("m"),
	py::arg("v"),
	py::arg("vr") 
    );

    m.def("dot_vectors", [] (math::Matrix m, Vectors const & v, Vectors & vr) {
		dot_vectors(m, v, vr); },
	py::arg("m"),
	py::arg("v"),
	py::arg("vr") 
    );

    m.def("get_quaternion_from_matrix", [] (Matrix const & m) -> Quaternion {
		 return get_quaternion_from_matrix(m); },
	py::arg("m") 
    );

    m.def("get_random_quaternion", [] () -> Quaternion {
		 return get_random_quaternion(); } 
    );

    m.def("matrix_from_str", [] (std::string const & s) ->  const Matrix {
		 return matrix_from_str(s); },
	py::arg("s") 
    );

    m.def("matrix_to_str", [] (math::Matrix m) ->  String {
		 return matrix_to_str(m); },
	py::arg("m") 
    );

    m.def("mean", [] (std::vector<double> const & a) -> double {
		 return mean(a); },
	py::arg("a") 
    );

    m.def("norm", [] (std::vector<double> const & v) ->  double {
		 return norm(v); },
	py::arg("v") 
    );

    m.def("pearson_coeff", [] (std::vector<double> const & x, std::vector<double> const & y) -> double {
		 return pearson_coeff(x, y); },
	py::arg("x"),
	py::arg("y") 
    );

    m.def("power_iteration", [] (std::vector<std::vector<double> > const & A, std::vector<double> & eigen_values, int num_simulations) {
		power_iteration(A, eigen_values, num_simulations); },
	py::arg("A"),
	py::arg("eigen_values"),
	py::arg("num_simulations") 
    );

    m.def("sqsum", [] (std::vector<double> const & a) -> double {
		 return sqsum(a); },
	py::arg("a") 
    );

    m.def("stdev", [] (std::vector<double> const & nums) -> double {
		 return stdev(nums); },
	py::arg("nums") 
    );

    m.def("sum", [] (std::vector<double> const & a) -> double {
		 return sum(a); },
	py::arg("a") 
    );

    m.def("vector_from_str", [] (std::string const & s) ->  const Vector {
		 return vector_from_str(s); },
	py::arg("s") 
    );

    m.def("vector_to_str", [] (math::Vector v) ->  const String {
		 return vector_to_str(v); },
	py::arg("v") 
    );

    m.def("vectors_from_str", [] (std::string const & s) ->  const Vectors {
		 return vectors_from_str(s); },
	py::arg("s") 
    );

    m.def("vectors_to_str", [] (math::Vectors vs) ->  String {
		 return vectors_to_str(vs); },
	py::arg("vs") 
    );
	// classes

        py::class_<AverageQuaternionCalculator, std::shared_ptr<AverageQuaternionCalculator>>(m, "AverageQuaternionCalculator")
		// ctors
		.def(py::init<>())
		// methods
		.def("add_quaternion",[] (AverageQuaternionCalculator  & ptr, Quaternion const & q) {
		ptr.add_quaternion(q); } )
		.def("get_average",[] (AverageQuaternionCalculator  & ptr) -> Quaternion {
		 return ptr.get_average(); } )
		;

        py::class_<AxisAngle, std::shared_ptr<AxisAngle>>(m, "AxisAngle")
		// public attributes
		.def_readwrite("angle", &AxisAngle::angle)
		.def_readwrite("axis", &AxisAngle::axis)
		;

        py::class_<Quaternion, std::shared_ptr<Quaternion>>(m, "Quaternion")
		// ctors
		.def(py::init<>())
		.def(py::init<double>())
		.def(py::init<double,double,double,double>())
		.def(py::init<Quaternion const &>())
		// methods
		.def("dot",[] (Quaternion const & ptr, Quaternion const & q) ->  double {
		 return ptr.dot(q); } )
		.def("get_rotation_matrix",[] (Quaternion  & ptr) -> Matrix {
		 return ptr.get_rotation_matrix(); } )
		.def("get_a",[] (Quaternion const & ptr) ->  double {
		 return ptr.get_a(); } )
		.def("get_b",[] (Quaternion const & ptr) ->  double {
		 return ptr.get_b(); } )
		.def("get_c",[] (Quaternion const & ptr) ->  double {
		 return ptr.get_c(); } )
		.def("get_d",[] (Quaternion const & ptr) ->  double {
		 return ptr.get_d(); } )
	    .def("to_str", [] (Quaternion const& q) -> String { auto ss = std::stringstream(); ss << q; return ss.str(); })
	    .def("at", [](Quaternion const& q, int index) -> decltype(q[index]) {return q[index];})
		 // operators
		//.def(py::self << std::ostream)
		//.def(py::self [] int())
		.def(py::self += double())
		.def(py::self *= double())
		;

        py::class_<SixDCoordinateBinner, std::shared_ptr<SixDCoordinateBinner>>(m, "SixDCoordinateBinner")
		// ctors
		.def(py::init<math::BoundingBox,math::Real6>())
		// methods
		.def("bin6",[] (SixDCoordinateBinner const & ptr, math::Real6 values) -> Bin6D {
		 return ptr.bin6(values); } )
		.def("bin_center_point",[] (SixDCoordinateBinner const & ptr, math::Bin6D bin) -> Real6 {
		 return ptr.bin_center_point(bin); } )
		.def("bin_to_values",[] (SixDCoordinateBinner const & ptr, math::Bin6D bin) -> Real6 {
		 return ptr.bin_to_values(bin); } )
		.def("bin_index",[] (SixDCoordinateBinner const & ptr, math::Real6 values) -> uint64_t {
		 return ptr.bin_index(values); } )
		.def("bin_from_index",[] (SixDCoordinateBinner  & ptr, uint64_t bin_index) -> Bin6D {
		 return ptr.bin_from_index(bin_index); } )
		.def("get_bounding_box",[] (SixDCoordinateBinner  & ptr) -> BoundingBox const & {
		 return ptr.get_bounding_box(); } )
		.def("get_bin_widths",[] (SixDCoordinateBinner  & ptr) -> Real6 const & {
		 return ptr.get_bin_widths(); } )
		;

        py::class_<SixDHistogram, std::shared_ptr<SixDHistogram>>(m, "SixDHistogram")
		// ctors
		.def(py::init<math::BoundingBox,math::Real6>())
		.def(py::init<Strings const &,SixDHistogramStrType const &>())
		.def(py::init<std::ifstream &>())
		// methods
		//.def("read",[] (SixDHistogram  & ptr ) { return ptr.read(>(&num), ); } )
		//.def("read",[] (SixDHistogram  & ptr, reinterpret_cast<char * >(&count), sizeof(count) ) -> in {
		 //return ptr.read(>(&count), ); } )
		.def("size",[] (SixDHistogram  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("add",[] (SixDHistogram  & ptr, math::Real6 values) {
		ptr.add(values); } )
		.def("contains",[] (SixDHistogram  & ptr, math::Real6 values) -> bool {
		 return ptr.contains(values); } )
		.def("get",[] (SixDHistogram  & ptr, math::Real6 values) -> uint64_t {
		 return ptr.get(values); } )
		.def("within_constraints",[] (SixDHistogram  & ptr, std::array<Real2, 6> const & constraints) -> uint64_t {
		 return ptr.within_constraints(constraints); } )
		.def("total_count",[] (SixDHistogram  & ptr) -> uint64_t {
		 return ptr.total_count(); } )
		.def("to_text_file",[] (SixDHistogram  & ptr, String const & fname) {
		ptr.to_text_file(fname); } )
		.def("to_binary_file",[] (SixDHistogram  & ptr, String const & fname, uint64_t cuttoff) {
		ptr.to_binary_file(fname, cuttoff); } )
		.def("output_binary",[] (SixDHistogram  & ptr, std::ofstream & out, uint64_t cuttoff) {
		ptr.output_binary(out, cuttoff); } )
		// public attributes
		//.def_readwrite("datas", &SixDHistogram::)
        ////.def_readwrite("", &SixDHistogram::)
		//.def_readwrite("lower", &SixDHistogram::lower)
		//.def_readwrite("", &SixDHistogram::)
		//.def_readwrite("upper", &SixDHistogram::upper)
		//.def_readwrite("bin_widths", &SixDHistogram::bin_widths)
		//.def_readwrite("", &SixDHistogram::)
		//.def_readwrite("binner_", &SixDHistogram::binner_)
		//.def_readwrite("num", &sum)
		//.def_readwrite("key", &SixDHistogram::key)
		//.def_readwrite("count", &SixDHistogram::count)
		//.def_readwrite("num", &SixDHistogram::num)
		//.def_readwrite("", &SixDHistogram::)
		;
//        py::class_<SixDHistogramStrType, std::shared_ptr<SixDHistogramStrType>>(m, "SixDHistogramStrType")
//		;
        py::enum_<SixDHistogramStrType>(m, "SixDHistogramStrType")
            .value("TEXT", SixDHistogramStrType::TEXT)
            .value("BINARY", SixDHistogramStrType::BINARY)
        ;

        py::class_<ThreeDCoordinateBinner, std::shared_ptr<ThreeDCoordinateBinner>>(m, "ThreeDCoordinateBinner")
		// ctors
		.def(py::init<math::BoundingBox,math::Real3>())
		// methods
		.def("bin3",[] (ThreeDCoordinateBinner const & ptr, Point const & values) -> Real3 {
		 return ptr.bin3(values); } )
		.def("bin_center_point",[] (ThreeDCoordinateBinner const & ptr, math::Real3 bin) -> Real3 {
		 return ptr.bin_center_point(bin); } )
		.def("bin_to_values",[] (ThreeDCoordinateBinner const & ptr, math::Bin6D bin) -> Real3 {
		 return ptr.bin_to_values(bin); } )
		.def("bin_index",[] (ThreeDCoordinateBinner const & ptr, Point const & values) -> uint32_t {
		 return ptr.bin_index(values); } )
		.def("bin_from_index",[] (ThreeDCoordinateBinner  & ptr, uint32_t bin_index) -> Real3 {
		 return ptr.bin_from_index(bin_index); } )
		.def("get_bounding_box",[] (ThreeDCoordinateBinner  & ptr) -> BoundingBox const & {
		 return ptr.get_bounding_box(); } )
		.def("get_bin_widths",[] (ThreeDCoordinateBinner  & ptr) -> Real3 const & {
		 return ptr.get_bin_widths(); } )
		;

        py::class_<ThreeDHistogram, std::shared_ptr<ThreeDHistogram>>(m, "ThreeDHistogram")
		// ctors
		.def(py::init<>())
		.def(py::init<math::BoundingBox,math::Real3>())
		// methods
		.def("setup",[] (ThreeDHistogram  & ptr, math::BoundingBox bounding_box, math::Real3 bin_widths) {
		ptr.setup(bounding_box, bin_widths); } )
		.def("size",[] (ThreeDHistogram  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("add",[] (ThreeDHistogram  & ptr, Point const & values) {
		ptr.add(values); } )
		.def("contains",[] (ThreeDHistogram  & ptr, Point const & values) -> bool {
		 return ptr.contains(values); } )
		.def("get",[] (ThreeDHistogram  & ptr, Point const & values) -> uint32_t {
		 return ptr.get(values); } )
		.def("write_histo_to_pdb",[] (ThreeDHistogram  & ptr, String const & pdb_name) {
		ptr.write_histo_to_pdb(pdb_name); } )
		;

        py::class_<Transform, std::shared_ptr<Transform>>(m, "Transform")
		// ctors
		.def(py::init<>())
		.def(py::init<Matrix const &,Point const &>())
		// methods
		//.def("dot",[] (Transform  & ptr, Matrix const & a, Matrix const & b, Transform & c) { ptr.dot(a, b, c); } )
		.def("xx",[] (Transform const & ptr) -> float {
		 return ptr.xx(); } )
		.def("xy",[] (Transform const & ptr) -> float {
		 return ptr.xy(); } )
		.def("xz",[] (Transform const & ptr) -> float {
		 return ptr.xz(); } )
		.def("yx",[] (Transform const & ptr) -> float {
		 return ptr.yx(); } )
		.def("yy",[] (Transform const & ptr) -> float {
		 return ptr.yy(); } )
		.def("yz",[] (Transform const & ptr) -> float {
		 return ptr.yz(); } )
		.def("zx",[] (Transform const & ptr) -> float {
		 return ptr.zx(); } )
		.def("zy",[] (Transform const & ptr) -> float {
		 return ptr.zy(); } )
		.def("zz",[] (Transform const & ptr) -> float {
		 return ptr.zz(); } )
		.def("px",[] (Transform const & ptr) -> float {
		 return ptr.px(); } )
		.def("py",[] (Transform const & ptr) -> float {
		 return ptr.py(); } )
		.def("pz",[] (Transform const & ptr) -> float {
		 return ptr.pz(); } )
		.def("xx",[] (Transform  & ptr) -> float {
		 return ptr.xx(); } )
		.def("xy",[] (Transform  & ptr) -> float {
		 return ptr.xy(); } )
		.def("xz",[] (Transform  & ptr) -> float {
		 return ptr.xz(); } )
		.def("yx",[] (Transform  & ptr) -> float {
		 return ptr.yx(); } )
		.def("yy",[] (Transform  & ptr) -> float {
		 return ptr.yy(); } )
		.def("yz",[] (Transform  & ptr) -> float {
		 return ptr.yz(); } )
		.def("zx",[] (Transform  & ptr) -> float {
		 return ptr.zx(); } )
		.def("zy",[] (Transform  & ptr) -> float {
		 return ptr.zy(); } )
		.def("zz",[] (Transform  & ptr) -> float {
		 return ptr.zz(); } )
		.def("px",[] (Transform  & ptr) -> float {
		 return ptr.px(); } )
		.def("py",[] (Transform  & ptr) -> float {
		 return ptr.py(); } )
		.def("pz",[] (Transform  & ptr) -> float {
		 return ptr.pz(); } )
		.def("xaxis",[] (Transform const & ptr) -> Vector {
		 return ptr.xaxis(); } )
		.def("yaxis",[] (Transform const & ptr) -> Vector {
		 return ptr.yaxis(); } )
		.def("zaxis",[] (Transform const & ptr) -> Vector {
		 return ptr.zaxis(); } )
		.def("translation",[] (Transform const & ptr) -> Vector {
		 return ptr.translation(); } )
		.def("rotation",[] (Transform const & ptr) ->  Matrix {
		 return ptr.rotation(); } )
		.def("rotation",[] (Transform  & ptr, Matrix const & m) {
		ptr.rotation(m); } )
		.def("translation",[] (Transform  & ptr, Vector const & v) {
		ptr.translation(v); } )
		;

        py::class_<_BoundingBox<Point>, std::shared_ptr<_BoundingBox<Point>>>(m, "_BoundingBox")
		// ctors
		.def(py::init<>())
		.def(py::init<Point const &>())
		.def(py::init<Point const &, Point const &>())
		.def(py::init<_BoundingBox<Point> const &>())
		// methods
//		.def("add",[] (_BoundingBox<Point>  & ptr, Point const & pp) {
//		ptr.add(pp); } )
		.def("reset",[] (_BoundingBox<Point>  & ptr, Point const & p) {
		ptr.reset(p); } )
		.def("expand",[] (_BoundingBox<Point>  & ptr, double const & scalar) {
		ptr.expand(scalar); } )
		.def("contract",[] (_BoundingBox<Point>  & ptr, double const & scalar) {
		ptr.contract(scalar); } )
		.def("translate",[] (_BoundingBox<Point>  & ptr, Point const & t) {
		ptr.translate(t); } )
		.def("intersects",[] (_BoundingBox<Point> const & ptr, _BoundingBox<Point> const & bb) ->  bool {
		 return ptr.intersects(bb); } )
		.def("contains",[] (_BoundingBox<Point> const & ptr, double const & x, double const & y, double const & z) ->  bool {
		 return ptr.contains(x, y, z); } )
		.def("contains",[] (_BoundingBox<Point> const & ptr, Point const & p) ->  bool {
		 return ptr.contains(p); } )
		.def("set_lower",[] (_BoundingBox<Point>  & ptr, Point const & p) {
		ptr.set_lower(p); } )
		.def("set_upper",[] (_BoundingBox<Point>  & ptr, Point const & p) {
		ptr.set_upper(p); } )
		.def("lower",[] (_BoundingBox<Point> const & ptr) ->  Point const & {
		 return ptr.lower(); } )
		.def("upper",[] (_BoundingBox<Point> const & ptr) ->  Point const & {
		 return ptr.upper(); } )
		// operators
		//.def(py::self = py::self)
		;

#define Value double
        py::class_<xyzMatrix<Value>, std::shared_ptr<xyzMatrix<Value>>>(m, "xyzMatrix")
		// ctors
		.def(py::init<>())
		.def(py::init<xyzMatrix<Value> const &>())
		.def(py::init<xyzMatrix<double> const &>())
		.def(py::init<const double &,const double &,const double &,const double &,const double &,const double &,const double &,const double &,const double &>())
		.def(py::init<Value const &>())
		.def(py::init<String const &>())
		// methods
		.def("transpose",[] (xyzMatrix<Value> & ptr) { return ptr.transpose(); } )
		//.def("dot",[] (xyzMatrix<Value>  & ptr, xyzMatrix<Value> const & a, xyzMatrix<Value> const & b, xyzMatrix<Value> & c) {
		//ptr.dot(a, b, c); } )
		.def("to_str",[] (xyzMatrix<Value> const & ptr) ->  String const {
		 return ptr.to_str(); } )
		.def("row",[] (xyzMatrix<Value>  & ptr, const int i, Vector const & v) {
		ptr.row(i, v); } )
		.def("row",[] (xyzMatrix<Value>  & ptr, const int i, std::vector<Value> const & v) {
		ptr.row(i, v); } )
		.def("identity",[] (xyzMatrix<Value>  & ptr) ->  xyzMatrix<Value> {
		 return ptr.identity(); } )
		.def("transpose",[] (xyzMatrix<Value>  & ptr) ->  xyzMatrix<Value> & {
		 return ptr.transpose(); } )
		.def("difference",[] (xyzMatrix<Value> const & ptr, xyzMatrix<Value> const & b) ->  float const {
		 return ptr.difference(b); } )
		.def("get_flip_orientation",[] (xyzMatrix<Value> const & ptr) ->  const xyzMatrix<Value> {
		 return ptr.get_flip_orientation(); } )
		.def("get_unitarize",[] (xyzMatrix<Value> const & ptr) ->  const xyzMatrix<Value> {
		 return ptr.get_unitarize(); } )
		.def("unitarize",[] (xyzMatrix<Value>  & ptr) {
		ptr.unitarize(); } )
		.def("xx",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.xx(); } )
		.def("xy",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.xy(); } )
		.def("xz",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.xz(); } )
		.def("yx",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.yx(); } )
		.def("yy",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.yy(); } )
		.def("yz",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.yz(); } )
		.def("zx",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.zx(); } )
		.def("zy",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.zy(); } )
		.def("zz",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
		 return ptr.zz(); } )
		.def("xx",[] (xyzMatrix<Value>  & ptr, Value const & xx_a) {
		ptr.xx(xx_a); } )
		.def("xy",[] (xyzMatrix<Value>  & ptr, Value const & xy_a) {
		ptr.xy(xy_a); } )
		.def("xz",[] (xyzMatrix<Value>  & ptr, Value const & xz_a) {
		ptr.xz(xz_a); } )
		.def("yx",[] (xyzMatrix<Value>  & ptr, Value const & yx_a) {
		ptr.yx(yx_a); } )
		.def("yy",[] (xyzMatrix<Value>  & ptr, Value const & yy_a) {
		ptr.yy(yy_a); } )
		.def("yz",[] (xyzMatrix<Value>  & ptr, Value const & yz_a) {
		ptr.yz(yz_a); } )
		.def("zx",[] (xyzMatrix<Value>  & ptr, Value const & zx_a) {
		ptr.zx(zx_a); } )
		.def("zy",[] (xyzMatrix<Value>  & ptr, Value const & zy_a) {
		ptr.zy(zy_a); } )
		.def("zz",[] (xyzMatrix<Value>  & ptr, Value const & zz_a) {
		ptr.zz(zz_a); } )
		.def("transposed",[] (xyzMatrix<Value> const & ptr) ->  xyzMatrix<Value> {
		 return ptr.transposed(); } )
		// operators
		.def(py::self += py::self )
		.def(py::self -= py::self )
		.def(py::self += Value() )
		.def(py::self -= Value() )
		.def(py::self + py::self)
		.def(py::self + py::self)
		.def(py::self + Value() )
		.def(py::self - py::self)
		.def(py::self - py::self)
		.def(py::self - Value())
		.def(py::self * py::self)
		.def(py::self * py::self)
		.def(py::self * Value())
		//.def(py::self / py::self)
		;

        py::class_<xyzVector<Value>, std::shared_ptr<xyzVector<Value>>>(m, "xyzVector<Value>")
		// ctors
		.def(py::init<>())
		.def(py::init<xyzVector<Value> const &>())
		.def(py::init<xyzVector<Value> const &>())
		.def(py::init<Value const &,Value const &,Value const &>())
		.def(py::init<std::vector<Value> const &>())
		.def(py::init<Value const &>())
		.def(py::init<String const &>())
		// methods
		.def("to_str",[] (xyzVector<Value> const & ptr) ->  String const {
		 return ptr.to_str(); } )
		.def("zero",[] (xyzVector<Value>  & ptr) ->  xyzVector<Value> & {
		 return ptr.zero(); } )
		.def("negate",[] (xyzVector<Value>  & ptr) ->  xyzVector<Value> & {
		 return ptr.negate(); } )
		.def("negated",[] (xyzVector<Value> const & ptr) ->  xyzVector<Value> {
		 return ptr.negated(); } )
		.def("negated",[] (xyzVector<Value> const & ptr, xyzVector<Value> & a) {
		ptr.negated(a); } )
//		.def("add",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b, xyzVector<Value> & r) {
//		ptr.add(a, b, r); } )
//		.def("add",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.add(v, t, r); } )
//		.def("add",[] (xyzVector<Value>  & ptr, Value const & t, xyzVector<Value> const & v, xyzVector<Value> & r) {
//		ptr.add(t, v, r); } )
//		.def("subtract",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b, xyzVector<Value> & r) {
//		ptr.subtract(a, b, r); } )
//		.def("subtract",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.subtract(v, t, r); } )
//		.def("subtract",[] (xyzVector<Value>  & ptr, Value const & t, xyzVector<Value> const & v, xyzVector<Value> & r) {
//		ptr.subtract(t, v, r); } )
//		.def("multiply",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.multiply(v, t, r); } )
//		.def("multiply",[] (xyzVector<Value>  & ptr, Value const & t, xyzVector<Value> const & v, xyzVector<Value> & r) {
//		ptr.multiply(t, v, r); } )
//		.def("divide",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.divide(v, t, r); } )
		.def("normalize",[] (xyzVector<Value>  & ptr) ->  xyzVector<Value> & {
		 return ptr.normalize(); } )
		.def("distance",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  Value {
		 return ptr.distance(v); } )
		.def("distance_squared",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  Value {
		 return ptr.distance_squared(v); } )
		.def("dot",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  Value {
		 return ptr.dot(v); } )
//		.def("dot_product",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b) -> friend  Value {
//		 return ptr.dot_product(a, b); } )
		.def("cross",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  xyzVector<Value> {
		 return ptr.cross(v); } )
//		.def("cross",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b) -> friend  xyzVector<Value> {
//		 return ptr.cross(a, b); } )
		.def("x",[] (xyzVector<Value> const & ptr) ->  Value const & {
		 return ptr.x(); } )
		.def("y",[] (xyzVector<Value> const & ptr) ->  Value const & {
		 return ptr.y(); } )
		.def("z",[] (xyzVector<Value> const & ptr) ->  Value const & {
		 return ptr.z(); } )
		.def("length",[] (xyzVector<Value> const & ptr) ->  Value {
		 return ptr.length(); } )
		.def("length_squared",[] (xyzVector<Value> const & ptr) ->  Value {
		 return ptr.length_squared(); } )
		.def("norm",[] (xyzVector<Value> const & ptr) ->  Value {
		 return ptr.norm(); } )
		.def("norm_squared",[] (xyzVector<Value> const & ptr) ->  Value {
		 return ptr.norm_squared(); } )
		.def("magnitude",[] (xyzVector<Value> const & ptr) ->  Value {
		 return ptr.magnitude(); } )
		.def("magnitude_squared",[] (xyzVector<Value> const & ptr) ->  Value {
		 return ptr.magnitude_squared(); } )
		.def("x",[] (xyzVector<Value>  & ptr, Value const & x_a) {
		ptr.x(x_a); } )
		.def("y",[] (xyzVector<Value>  & ptr, Value const & y_a) {
		ptr.y(y_a); } )
		.def("z",[] (xyzVector<Value>  & ptr, Value const & z_a) {
		ptr.z(z_a); } )
		// operators
//		.def(py::self = py::self)
//		.def(py::self = py::self  )
		.def(py::self += py::self  )
		.def(py::self -= py::self  )
		//.def(py::self = Value()  )
		.def(py::self += Value()  )
		.def(py::self -= Value()  )
		.def(py::self *= Value()  )
		.def(py::self /= Value()  )
		.def(py::self + py::self)
		.def(py::self + py::self)
		.def(py::self + Value()  )
		.def(py::self - py::self)
		.def(py::self - py::self)
		.def(py::self - Value()  )
//		.def(py::self * py::self)
		.def(py::self * Value()  )
//		.def(py::self / py::self)
        .def("at", [] (xyzVector<Value> const& vec, int index) -> Value { return vec[index]; })
//        .def(py::self [] int )
//		.def(py::self [] int )
//		.def(py::self () int )
//		.def(py::self () int )
		.def(py::self == py::self)
		.def(py::self != py::self)
		;
#undef Value


}