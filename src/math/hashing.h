//
// Created by Joseph Yesselman on 2/16/18.
//

#ifndef TEST_HASHING_H
#define TEST_HASHING_H

#include <stdint.h>
#include <stdio.h>

#include <array>
#include <map>

#include <math/vector_3.hpp>

namespace math {

template <typename T> class _BoundingBox {

public: // types
  typedef T PointPosition; // PointPosition is a synonym for T

public: // construct/destruct
  inline _BoundingBox() = default;

  inline _BoundingBox(PointPosition const &pp) : _lower(pp), _upper(pp) {}

  inline _BoundingBox(PointPosition const &lower, PointPosition const &upper)
      : _lower(lower), _upper(upper) {
    if (lower == upper) {
      String msg = "Boundaries are the same! Try different bounds.";
      base::log_and_throw<base::MathException>(msg);
    }
  }

  /// @brief copy constructor
  inline _BoundingBox(_BoundingBox const &bb)
      : _lower(bb._lower), _upper(bb._upper) {}

  /// @brief default destructor
  inline ~_BoundingBox() = default;

public: // assignment
  /// @brief copy assignment
  inline _BoundingBox &operator=(_BoundingBox const &bb) {
    if (this != &bb) {
      _lower = bb._lower;
      _upper = bb._upper;
    }
    return *this;
  }

public: // box management
  /// @brief add a point, expands bounds if necessary
  inline void add(PointPosition const &pp) {
    _lower.min(pp);
    _upper.max(pp);
  }

  /// @brief reset corners
  inline void reset(PointPosition const &p = PointPosition()) {
    _lower = p;
    _upper = p;
  }

  /// @brief expand box corners (additive)
  template <typename U> inline void expand(U const &scalar) {
    _lower -= scalar;
    _upper += scalar;
  }

  /// @brief contract box corners (subtractive)
  template <typename U> inline void contract(U const &scalar) {
    _lower += scalar;
    _upper -= scalar;
  }

  /// @brief translate bounding box
  inline void translate(PointPosition const &t) {
    _lower += t;
    _upper += t;
  }

public: // box query
  /// @brief intersects another bounding box?
  inline bool intersects(_BoundingBox const &bb) const {
    return !(_lower.x() > bb._upper.x() || bb._lower.x() > _upper.x() ||
             _lower.y() > bb._upper.y() || bb._lower.y() > _upper.y() ||
             _lower.z() > bb._upper.z() || bb._lower.z() > _upper.z());
  }

  /// @brief is point contained within this bounding box?
  template <typename U>
  inline bool contains(U const &x, U const &y, U const &z) const {
    return _lower.x() <= x && _lower.y() <= y && _lower.z() <= z &&
           x <= _upper.x() && y <= _upper.y() && z <= _upper.z();
  }

  /// @brief is point contained within this bounding box?
  inline bool contains(PointPosition const &p) const {
    return contains(p.x(), p.y(), p.z());
  }

public: // setters
  /// @brief set lower corner
  inline void set_lower(PointPosition const &p) { _lower = p; }

  /// @brief set upper corner
  inline void set_upper(PointPosition const &p) { _upper = p; }

public: // getters
  /// @brief get lower corner
  inline PointPosition const &lower() const { return _lower; }

  /// @brief get upper corner
  inline PointPosition const &upper() const { return _upper; }

private:                // data
  PointPosition _lower; // lower corner
  PointPosition _upper; // upper corner
};

typedef _BoundingBox<Vector3> BoundingBox;
typedef std::array<size_t, 3> Size3;
typedef std::array<size_t, 6> Size6;
typedef std::array<double, 2> Real2;
typedef std::array<double, 3> Real3;
typedef std::array<double, 6> Real6;
typedef std::array<double, 6> Bin6D;

class SixDCoordinateBinner {
public:
  SixDCoordinateBinner(BoundingBox const &bounding_box, Real6 const &bin_widths)
      : _bounding_box(bounding_box), _bin_widths(bin_widths),
        _bin_values(Real6()) {
    auto span = _bounding_box.upper() - _bounding_box.lower();
    auto new_upper = _bounding_box.upper();

    _dimsizes[0] = static_cast<size_t>(span.get_x() / _bin_widths[0]);
    _dimsizes[1] = static_cast<size_t>(span.get_y() / _bin_widths[1]);
    _dimsizes[2] = static_cast<size_t>(span.get_z() / _bin_widths[2]);

    if (_dimsizes[0] * _bin_widths[0] < span.get_x()) {
      _dimsizes[0] += 1;
      new_upper.set_x(_bounding_box.lower().get_x() +
                      _dimsizes[0] * _bin_widths[0]);
    } else if (_dimsizes[1] * _bin_widths[1] < span.get_y()) {
      _dimsizes[1] += 1;
      new_upper.set_y(_bounding_box.lower().get_y() +
                      _dimsizes[1] * _bin_widths[1]);
    } else if (_dimsizes[2] * _bin_widths[2] < span.get_z()) {
      _dimsizes[2] += 1;
      new_upper.set_z(_bounding_box.lower().get_z() +
                      _dimsizes[2] * _bin_widths[2]);
    }
    _bounding_box.set_upper(new_upper);

    for (int ii = 3; ii <= 5; ++ii) {
      _dimsizes[ii] = static_cast<size_t>(360.0 / _bin_widths[ii]);
    }
    if (_dimsizes[5] == 0) {
      _dimsizes[5] = 1;
    }

    _dimprods[5] = 1;
    for (int ii = 4; ii >= 0; --ii) {
      _dimprods[ii] = _dimprods[ii + 1] * _dimsizes[ii + 1];
    }

    for (int ii = 0; ii < 6; ++ii) {
      _halfbin_widths[ii] = _bin_widths[ii] / 2;
    }
  }

  Bin6D bin6(Real6 const &values) const {
    auto xyzcoord = Vector3(values[0], values[1], values[2]);
    assert(_bounding_box.contains(xyzcoord));

    auto from_corner = xyzcoord - _bounding_box.lower();
    Bin6D bins;

    bins[0] = static_cast<size_t>(from_corner.get_x() / _bin_widths[0]);
    if (bins[0] == _dimsizes[0]) {
      bins[0] -= 1;
    }

    bins[1] = static_cast<size_t>(from_corner.get_y() / _bin_widths[1]);
    if (bins[1] == _dimsizes[1]) {
      bins[1] -= 1;
    }

    bins[2] = static_cast<size_t>(from_corner.get_z() / _bin_widths[2]);
    if (bins[2] == _dimsizes[2]) {
      bins[2] -= 1;
    }

    auto euler = _wrap_euler_angles(values);

    bins[3] = static_cast<size_t>(euler[0] / _bin_widths[3]) % _dimsizes[3];
    bins[4] = static_cast<size_t>(euler[1] / _bin_widths[4]) % _dimsizes[4];
    bins[5] = static_cast<size_t>(euler[2] / _bin_widths[5]) % _dimsizes[5];
    return bins;
  }

  Real6 bin_center_point(Bin6D const &bin) const {
    Real6 center;

    center[0] = _bounding_box.lower().get_x() + bin[0] * _bin_widths[0] +
                _halfbin_widths[0];
    center[1] = _bounding_box.lower().get_y() + bin[1] * _bin_widths[1] +
                _halfbin_widths[1];
    center[2] = _bounding_box.lower().get_z() + bin[2] * _bin_widths[2] +
                _halfbin_widths[2];
    center[3] = bin[3] * _bin_widths[3] + _halfbin_widths[3];
    center[4] = bin[4] * _bin_widths[4] + _halfbin_widths[4];
    center[5] = bin[5] * _bin_widths[5] + _halfbin_widths[5];

    return center;
  }

  Real6 bin_to_values(Bin6D const &bin) const {
    Real6 values;

    values[0] = _bounding_box.lower().get_x() + bin[0] * _bin_widths[0];
    values[1] = _bounding_box.lower().get_y() + bin[1] * _bin_widths[1];
    values[2] = _bounding_box.lower().get_z() + bin[2] * _bin_widths[2];
    values[3] = bin[3] * _bin_widths[3];
    values[4] = bin[4] * _bin_widths[4];
    values[5] = bin[5] * _bin_widths[5];

    return values;
  }

  uint64_t bin_index(Real6 const &values) const {
    auto bin = bin6(values);

    uint64_t const A = bin[0] * _dimprods[0] + bin[1] * _dimprods[1] +
                       bin[2] * _dimprods[2] + bin[3] * _dimprods[3] +
                       bin[4] * _dimprods[4] + bin[5] * _dimprods[5];
    return A;
  }

  Bin6D bin_from_index(uint64_t bin_index) {
    Bin6D bin;
    for (int ii = 0; ii < 6; ++ii) {
      bin[ii] = bin_index / _dimprods[ii];
      bin_index = bin_index % _dimprods[ii];
    }

    return bin;
  }

public: // getters
  BoundingBox const &get_bounding_box() { return _bounding_box; }

  Real6 const &get_bin_widths() { return _bin_widths; }

private:
  Real3 _wrap_euler_angles(Real6 const &values) const {
    auto euler = Real3{values[3], values[4], values[5]};
    for (int i = 0; i < euler.size(); i++) {
      if (euler[i] > 360) {
        euler[i] -= 360;
      }
      if (euler[i] < 0) {
        euler[i] += 360;
      }
    }
    return euler;
  }

private:
  BoundingBox _bounding_box;
  Size6 _dimsizes;
  Size6 _dimprods;
  Real6 _bin_widths;
  Real6 _halfbin_widths;
  Real6 _bin_values;
};

class ThreeDCoordinateBinner {
public:

  /// @brief - constructor
  ThreeDCoordinateBinner(BoundingBox const &bounding_box,
                         Real3 const &bin_widths)
      : _bounding_box(bounding_box), _bin_widths(bin_widths) {
    auto span = _bounding_box.upper() - _bounding_box.lower();
    auto new_upper = _bounding_box.upper();

    _dimsizes[0] = static_cast<size_t>(span.get_x() / _bin_widths[0]);
    _dimsizes[1] = static_cast<size_t>(span.get_y() / _bin_widths[1]);
    _dimsizes[2] = static_cast<size_t>(span.get_z() / _bin_widths[2]);

    if (_dimsizes[0] * _bin_widths[0] < span.get_x()) {
      _dimsizes[0] += 1;
      new_upper.set_x(_bounding_box.lower().get_x() +
                      _dimsizes[0] * _bin_widths[0]);
    } else if (_dimsizes[1] * _bin_widths[1] < span.get_y()) {
      _dimsizes[1] += 1;
      new_upper.set_y(_bounding_box.lower().get_y() +
                      _dimsizes[1] * _bin_widths[1]);
    } else if (_dimsizes[2] * _bin_widths[2] < span.get_z()) {
      _dimsizes[2] += 1;
      new_upper.set_z(_bounding_box.lower().get_z() +
                      _dimsizes[2] * _bin_widths[2]);
    }
    _bounding_box.set_upper(new_upper);

    if (_dimsizes[2] == 0) {
      _dimsizes[2] = 1;
    }
    _dimprods[2] = 1;
    for (int ii = 1; ii >= 0; --ii) {
      _dimprods[ii] = _dimprods[ii + 1] * _dimsizes[ii + 1];
    }

    for (int ii = 0; ii < 3; ++ii) {
      _halfbin_widths[ii] = _bin_widths[ii] / 2;
    }
  }

  Real3 bin3(Vector3 const &values) const {
    // assert(_bounding_box.contains(values));

    auto from_corner = values - _bounding_box.lower();
    Real3 bins;

    bins[0] = static_cast<size_t>(from_corner.get_x() / _bin_widths[0]);

    if (bins[0] == _dimsizes[0]) {
      bins[0] -= 1;
    }

    bins[1] = static_cast<size_t>(from_corner.get_y() / _bin_widths[1]);

    if (bins[1] == _dimsizes[1]) {
      bins[1] -= 1;
    }

    bins[2] = static_cast<size_t>(from_corner.get_z() / _bin_widths[2]);

    if (bins[2] == _dimsizes[2]) {
      bins[2] -= 1;
    }

    return bins;
  }

  Real3 bin_center_point(Real3 const &bin) const {
    Real3 center;
    center[0] = _bounding_box.lower().get_x() + bin[0] * _bin_widths[0] +
                _halfbin_widths[0];
    center[1] = _bounding_box.lower().get_y() + bin[1] * _bin_widths[1] +
                _halfbin_widths[1];
    center[2] = _bounding_box.lower().get_z() + bin[2] * _bin_widths[2] +
                _halfbin_widths[2];
    return center;
  }

  Real3 bin_to_values(Bin6D const &bin) const {
    Real3 values;
    values[0] = _bounding_box.lower().get_x() + bin[0] * _bin_widths[0];
    values[1] = _bounding_box.lower().get_y() + bin[1] * _bin_widths[1];
    values[2] = _bounding_box.lower().get_z() + bin[2] * _bin_widths[2];
    return values;
  }

  uint32_t bin_index(Vector3 const &values) const {
    auto bin = bin3(values);

    uint32_t const A =
        bin[0] * _dimprods[0] + bin[1] * _dimprods[1] + bin[2] * _dimprods[2];
    return A;
  }

  Real3 bin_from_index(uint32_t bin_index) {
    Real3 bin;
    for (int ii = 0; ii < 3; ++ii) {
      bin[ii] = bin_index / _dimprods[ii];
      bin_index = bin_index % _dimprods[ii];
    }

    return bin;
  }

public: // getters
  BoundingBox const &get_bounding_box() { return _bounding_box; }

  Real3 const &get_bin_widths() { return _bin_widths; }

private:
  BoundingBox _bounding_box;
  Size3 _dimsizes;
  Size3 _dimprods;
  Real3 _bin_widths;
  Real3 _halfbin_widths;
};

enum class SixDHistogramStrType { TEXT, BINARY };

class SixDHistogram {
public:
  SixDHistogram(BoundingBox const &bounding_box, Real6 const &bin_widths)
      : _binner(SixDCoordinateBinner(bounding_box, bin_widths)) {}

  SixDHistogram(Strings const &s, SixDHistogramStrType const &type)
      : _binner(BoundingBox(), Real6{0.1, 0.1, 0.1, 0.1, 0.1, 0.1}) {
    if (type == SixDHistogramStrType::TEXT) {
      _setup_from_text(s);
    }
  }

  SixDHistogram(std::ifstream &in)
      : _binner(BoundingBox(), Real6{0.1, 0.1, 0.1, 0.1, 0.1, 0.1}) {
    auto datas = Real3();
    for (int i = 0; i < 3; i++) {
      in.read(reinterpret_cast<char *>(&datas[i]), sizeof(datas[i]));
    }
    auto lower = Vector3(datas[0], datas[1], datas[2]);
    for (int i = 0; i < 3; i++) {
      in.read(reinterpret_cast<char *>(&datas[i]), sizeof(datas[i]));
    }
    auto upper = Vector3(datas[0], datas[1], datas[2]);
    auto bin_widths = Real6();
    for (int i = 0; i < 6; i++) {
      in.read(reinterpret_cast<char *>(&bin_widths[i]), sizeof(bin_widths[i]));
    }
    _binner = SixDCoordinateBinner(BoundingBox(lower, upper), bin_widths);
    uint64_t num;
    uint64_t key, count;
    in.read(reinterpret_cast<char *>(&num), sizeof(num));
    for (int i = 0; i < num; i++) {
      in.read(reinterpret_cast<char *>(&key), sizeof(key));
      in.read(reinterpret_cast<char *>(&count), sizeof(count));
      _stored_values[key] = count;
    }
  }

private:
  void _setup_from_text(Strings const &s) {
    auto lower = vector_from_str(s[0]);
    auto upper = vector_from_str(s[1]);
    auto spl = base::string::split(s[2], " ");
    auto bin_widths = Real6();
    auto bb = BoundingBox(lower, upper);
    for (int i = 0; i < 6; i++) {
      bin_widths[i] = std::stod(spl[i]);
    }
    _binner = SixDCoordinateBinner(bb, bin_widths);
    auto values = Real6();
    for (int i = 4; i < s.size(); i++) {
      if (s[i].length() < 5) {
        break;
      }
      spl = base::string::split(s[i], ",");
      for (int i = 0; i < 6; i++) {
        values[i] = std::stod(spl[i]);
      }
      auto bin_index = _binner.bin_index(values);
      _stored_values[bin_index] = std::stoull(spl[6]);
    }
  }

public:
  inline size_t size() { return _stored_values.size(); }

public:
  void add(Real6 const &values) {
    auto bin_index = _binner.bin_index(values);
    if (_stored_values.find(bin_index) == _stored_values.end()) {
      _stored_values[bin_index] = 0;
    }
    _stored_values[bin_index] += 1;
  }

public:
  bool contains(Real6 const &values) {
    auto bin_index = _binner.bin_index(values);
    if (_stored_values.find(bin_index) == _stored_values.end()) {
      return false;
    } else {
      return true;
    }
  }

  uint64_t get(Real6 const &values) {
    auto bin_index = _binner.bin_index(values);
    return _stored_values[bin_index];
  }

  uint64_t within_constraints(std::array<Real2, 6> const &constraints) {
    auto values = Real6();
    auto bins = Real6();
    auto total = 0;
    auto fail = 0;
    double max_distance = 6;
    double dist = 0;
    for (auto const &kv : _stored_values) {
      bins = _binner.bin_from_index(kv.first);
      values = _binner.bin_to_values(bins);

      fail = 0;
      for (int i = 0; i < 6; i++) {
        if (i != 3 &&
            (constraints[i][0] > values[i] || values[i] > constraints[i][1])) {
          fail = 1;
          break;
        }
        if (i == 3 &&
            (constraints[i][0] < values[i] && values[i] < constraints[i][1])) {
          fail = 1;
          break;
        }
      }

      if (!fail) {
        //    dist = sqrt(values[0]*values[0] + values[1]*values[1] +
        //    values[2]*values[2]); if(dist > max_distance) { continue;}
        total += kv.second;
      }
    }
    return total;
  }

  uint64_t total_count() {
    uint64_t count = 0;
    for (auto const &kv : _stored_values) {
      count += kv.second;
    }
    return count;
  }

public:
  void to_text_file(String const &fname) {
    std::ofstream out;
    out.open(fname);
    auto &bb = _binner.get_bounding_box();
    out << bb.lower().get_str() << std::endl;
    out << bb.upper().get_str() << std::endl;
    for (auto const &bin_width : _binner.get_bin_widths()) {
      out << bin_width << " ";
    }
    out << std::endl;

    out << "x,y,z,a,b,g,count\n";
    for (auto const &kv : _stored_values) {
      auto bin_index = kv.first;
      auto bin = _binner.bin_from_index(bin_index);
      auto values = _binner.bin_to_values(bin);
      for (auto const &v : values) {
        out << v << ",";
      }
      out << kv.second << "\n";
    }
    out.close();
  }

  void to_binary_file(String const &fname, uint64_t cuttoff = 0) {
    std::ofstream out;
    out.open(fname, std::ios::binary);
    output_binary(out, cuttoff);
    out.close();
  }

  void output_binary(std::ofstream &out, uint64_t cuttoff) {
    auto lower = _binner.get_bounding_box().lower();
    auto upper = _binner.get_bounding_box().upper();
    double x = lower.get_x();
    double y = lower.get_y();
    double z = lower.get_z();
    double x1 = upper.get_x();
    double y1 = upper.get_y();
    double z1 = upper.get_z();
    out.write((const char *)(&x), sizeof(lower.get_x()));
    out.write((const char *)(&y), sizeof(lower.get_y()));
    out.write((const char *)(&z), sizeof(lower.get_z()));
    out.write((const char *)(&x1), sizeof(upper.get_x()));
    out.write((const char *)(&y1), sizeof(upper.get_y()));
    out.write((const char *)(&z1), sizeof(upper.get_z()));

    for (int i = 0; i < 6; i++) {
      out.write((const char *)&_binner.get_bin_widths()[i],
                sizeof(_binner.get_bin_widths()[i]));
    }

    uint64_t num = 0;
    for (auto const &kv : _stored_values) {
      // if(kv.second <= cuttoff) { continue; }
      num += 1;
    }
    out.write((const char *)&num, sizeof(num));
    for (auto const &kv : _stored_values) {
      // if(kv.second <= cuttoff) { continue; }
      out.write((const char *)&kv.first, sizeof(kv.first));
      out.write((const char *)&kv.second, sizeof(kv.second));
    }
  }

private:
  SixDCoordinateBinner _binner;
  std::map<uint64_t, uint64_t> _stored_values;
};

class ThreeDHistogram {
public:
  ThreeDHistogram() {}

  ThreeDHistogram(BoundingBox const &bounding_box, Real3 const &bin_widths)
      : _binner(std::make_shared<ThreeDCoordinateBinner>(bounding_box,
                                                         bin_widths)) {}

  void setup(BoundingBox const &bounding_box, Real3 const &bin_widths) {
    _binner =
        std::make_shared<ThreeDCoordinateBinner>(bounding_box, bin_widths);
  }

  inline size_t size() { return _stored_values.size(); }

  void add(Vector3 const &values) {
    auto bin_index = _binner->bin_index(values);
    if (_stored_values.find(bin_index) == _stored_values.end()) {
      _stored_values[bin_index] = 0;
    }
    _stored_values[bin_index] += 1;
  }

  bool contains(Vector3 const &values) {
    auto bin_index = _binner->bin_index(values);
    if (_stored_values.find(bin_index) == _stored_values.end()) {
      return false;
    } else {
      return true;
    }
  }

  uint32_t get(Vector3 const &values) {
    auto bin_index = _binner->bin_index(values);
    return _stored_values[bin_index];
  }

  void write_histo_to_pdb(String const &pdb_name) {
    int i = 1;
    std::ofstream out;
    out.open(pdb_name);
    auto bins = Real3();
    auto p = Real3();
    for (auto const &kv : _stored_values) {
      bins = _binner->bin_from_index(kv.first);
      p = _binner->bin_center_point(bins);
      char buffer[200];
      std::sprintf(
          buffer,
          "ATOM %6d  P   C   A   1 %11.3f%8.3f%8.3f  1.00 62.18           P\n",
          i, p[0], p[1], p[2]);
      out << String(buffer);
      i++;
    }
    out.close();
  }

private:
  std::shared_ptr<ThreeDCoordinateBinner> _binner;
  std::map<uint32_t, int> _stored_values;
};

} // namespace math

#endif // TEST_HASHING_H
