//
// Created by Joseph Yesselman on 2/16/18.
//

#ifndef TEST_HASHING_H
#define TEST_HASHING_H

#include <map>
#include <array>
#include "math/xyz_vector.h"

template< typename T>
class _BoundingBox {
public: // types
    typedef T PointPosition;

public: // construct/destruct
    inline
    _BoundingBox() = default;

    inline
    _BoundingBox(
            PointPosition const & pp):
            lower_( pp ),
            upper_( pp ) {}

    inline
    _BoundingBox(
            PointPosition const & lower,
            PointPosition const & upper):
            lower_( lower ),
            upper_( upper ) {}

    /// @brief copy constructor
    inline
    _BoundingBox(
            _BoundingBox const & bb):
            lower_( bb.lower_ ),
            upper_( bb.upper_ ) {}

    /// @brief default destructor
    inline
    ~_BoundingBox() = default;

public: // assignment

    /// @brief copy assignment
    inline
    _BoundingBox &
    operator =(
            _BoundingBox const & bb) {
        if ( this != &bb ) {
            lower_ = bb.lower_;
            upper_ = bb.upper_;
        }
        return *this;
    }

public: // box management

    /// @brief add a point, expands bounds if necessary
    inline
    void
    add(
            PointPosition const & pp) {
        lower_.min( pp );
        upper_.max( pp );
    }


    /// @brief reset corners
    inline
    void
    reset(
            PointPosition const & p = PointPosition()) {
        lower_ = p;
        upper_ = p;
    }

    /// @brief expand box corners (additive)
    template< typename U >
    inline
    void
    expand(
            U const & scalar
    )
    {
        lower_ -= scalar;
        upper_ += scalar;
    }

    // @brief contract box corners (subtractive)
    template< typename U >
    inline
    void
    contract(
            U const & scalar
    )
    {
        lower_ += scalar;
        upper_ -= scalar;
    }

    /// @brief translate bounding box
    inline
    void
    translate(
            PointPosition const & t) {
        lower_ += t;
        upper_ += t;
    }

public: // box query

    /// @brief intersects another bounding box?
    inline
    bool
    intersects(
            _BoundingBox const & bb) const {
        return !( lower_.x() > bb.upper_.x() || bb.lower_.x() > upper_.x() ||
                  lower_.y() > bb.upper_.y() || bb.lower_.y() > upper_.y() ||
                  lower_.z() > bb.upper_.z() || bb.lower_.z() > upper_.z() );
    }

    /// @brief is point contained within this bounding box?
    template< typename U >
    inline
    bool
    contains(
            U const & x,
            U const & y,
            U const & z) const {
        return lower_.x() <= x && lower_.y() <= y && lower_.z() <= z &&
               x <= upper_.x() && y <= upper_.y() && z <= upper_.z();
    }

    /// @brief is point contained within this bounding box?
    inline
    bool
    contains(
            PointPosition const & p) const {
        return contains( p.x(), p.y(), p.z() );
    }


public: // setters

    /// @brief set lower corner
    inline
    void
    set_lower(
            PointPosition const & p) {
        lower_ = p;
    }

    /// @brief set upper corner
    inline
    void
    set_upper(PointPosition const & p) {
        upper_ = p;
    }

public: // getters

    /// @brief get lower corner
    inline
    PointPosition const &
    lower() const {
        return lower_;
    }

    /// @brief get upper corner
    inline
    PointPosition const &
    upper() const {
        return upper_;
    }



private: // data

    PointPosition lower_; // lower corner
    PointPosition upper_; // upper corner

};

typedef _BoundingBox<Point> BoundingBox;
typedef std::array<size_t, 3> Size3;
typedef std::array<size_t, 6> Size6;
typedef std::array<double, 3> Real3;
typedef std::array<double, 6> Real6;
typedef std::array<double, 6> Bin6D;

class SixDCoordinateBinner {
public:
    SixDCoordinateBinner(
            BoundingBox const & bounding_box,
            Real6 const & bin_widths):
            bounding_box_(bounding_box),
            bin_widths_(bin_widths) {

        auto span = bounding_box_.upper() - bounding_box_.lower();
        auto new_upper = bounding_box_.upper();

        for (int ii = 0; ii < 3; ++ii) {
            dimsizes_[ii] = static_cast<size_t>( span[ii] / bin_widths_[ii] );
            if (dimsizes_[ii] * bin_widths_[ii] < span[ii]) {
                dimsizes_[ii] += 1;
                new_upper[ii] = bounding_box_.lower()[ii] + dimsizes_[ii] * bin_widths_[ii];
            }
        }
        bounding_box_.set_upper(new_upper);

        for (int ii = 3; ii <= 5; ++ii) {
            dimsizes_[ii] = static_cast<size_t> ( 360.0 / bin_widths_[ii] );
        }
        if (dimsizes_[5] == 0) { dimsizes_[5] = 1; }

        dimprods_[5] = 1;
        for (int ii = 4; ii >= 0; --ii) {
            dimprods_[ii] = dimprods_[ii + 1] * dimsizes_[ii + 1];
        }

        for (int ii = 0; ii < 6; ++ii) { halfbin_widths_[ii] = bin_widths_[ii] / 2; }

    }

    Bin6D
    bin6(
            Real6 const & values) const  {
        auto xyzcoord = Point(values[0], values[1], values[2]);
        assert(bounding_box_.contains(xyzcoord));

        auto from_corner = xyzcoord - bounding_box_.lower();
        Bin6D bins;

        bins[0] = static_cast< size_t > ( from_corner.x() / bin_widths_[0] );
        if ( bins[0] == dimsizes_[0] ) { bins[0] -= 1; }

        bins[1] = static_cast< size_t > ( from_corner.y() / bin_widths_[1] );
        if ( bins[1] == dimsizes_[1] ) { bins[1] -= 1; }

        bins[2] = static_cast< size_t > ( from_corner.z() / bin_widths_[2] );
        if ( bins[2] == dimsizes_[2] ) { bins[2] -= 1; }

        auto euler = _wrap_euler_angles( values );

        bins[3] = static_cast<size_t> ( euler[ 0 ] / bin_widths_[ 3 ] ) % dimsizes_[ 3 ];
        bins[4] = static_cast<size_t> ( euler[ 1 ] / bin_widths_[ 4 ] ) % dimsizes_[ 4 ];
        bins[5] = static_cast<size_t> ( euler[ 2 ] / bin_widths_[ 5 ] ) % dimsizes_[ 5 ];
        return bins;
    }

    Real6
    bin_center_point(
            Bin6D const & bin) const {
        Real6 center;
        for ( int ii = 0; ii < 3; ++ii ) {
            center[ii] = bounding_box_.lower()[ii] + bin[ii] * bin_widths_[ii] + halfbin_widths_[ii];
        }
        for ( int ii = 0; ii < 3; ++ii ) {
            center[ii + 3] = bin[ii + 3] * bin_widths_[ii + 3] + halfbin_widths_[ii + 3];
        }
        return center;
    }

    Real6
    bin_to_values(
            Bin6D const & bin) const {
        Real6 values;
        for ( int ii = 0; ii < 3; ++ii ) {
            values[ii] = bounding_box_.lower()[ii] + bin[ii] * bin_widths_[ii];
        }
        for ( int ii = 0; ii < 3; ++ii ) {
            values[ii + 3] = bin[ii + 3] * bin_widths_[ii + 3];
        }
        return values;
    }

    uint64_t
    bin_index(
            Real6 const & values) const {
        auto bin = bin6(values);

        uint64_t const A =
                bin[ 0 ] * dimprods_[ 0 ] +
                bin[ 1 ] * dimprods_[ 1 ] +
                bin[ 2 ] * dimprods_[ 2 ] +
                bin[ 3 ] * dimprods_[ 3 ] +
                bin[ 4 ] * dimprods_[ 4 ] +
                bin[ 5 ] * dimprods_[ 5 ];
        return A;
    }

    Bin6D
    bin_from_index(
            uint64_t bin_index) {
        Bin6D bin;
        for (int ii = 0; ii < 6; ++ii ) {
            bin[ii] = bin_index / dimprods_[ii];
            bin_index = bin_index % dimprods_[ii];
        }

        return bin;
    }

public: // getters

    BoundingBox const &
    get_bounding_box() { return bounding_box_; }

    Real6 const &
    get_bin_widths() { return bin_widths_; }


private:
    Real3
    _wrap_euler_angles(
            Real6 const & values) const {

        auto euler = Real3{values[3], values[4], values[5]};
        for(int i = 0; i < euler.size(); i++) {
            if(euler[i] > 360) { euler[i] -= 360; }
            if(euler[i] < 0)   { euler[i] += 360; }
        }


        /*if ( ( values[ 5 ] < euler_offsets_[ 2 ] || values[ 5 ] >= 180.0 - euler_offsets_[ 2 ] )) {
            /// The idea is that, when we're wrapping theta, half of the points have to stay fixed so they end up in the same
            /// bin: if all points wrapped, then none would land in the same bin.
            euler[ 0 ] = values[ 3 ] - euler_offsets_[ 0 ] - 180.0;
            euler[ 1 ] = values[ 4 ] - euler_offsets_[ 1 ] - 180.0;
            if ( euler[ 1 ] < 0 ) { euler[ 1 ] += 360; }
        } else {
            /// Leave phi/psi in their usual positions
            euler[ 0 ] = values[ 3 ] - euler_offsets_[ 0 ];
            if ( euler[ 0 ] < 0 ) { euler[ 0 ] += 360; }

            euler[ 1 ] = values[ 4 ] - euler_offsets_[ 1 ];
            if ( euler[ 1 ] < 0 ) { euler[ 1 ] += 360; }
        }*/
        return euler;


    }

private:
    BoundingBox bounding_box_;
    Size6 dimsizes_;
    Size6 dimprods_;
    Real6 bin_widths_;
    Real6 halfbin_widths_;


};

enum class SixDHistogramStrType {
    TEXT,
    BINARY
};

class SixDHistogram {
public:
    SixDHistogram(
            BoundingBox const & bounding_box,
            Real6 const & bin_widths):
            binner_ (SixDCoordinateBinner(bounding_box, bin_widths)) {}

    SixDHistogram(
            Strings const & s,
            SixDHistogramStrType const & type):
            binner_(BoundingBox(), Real6{0.1, 0.1, 0.1, 0.1, 0.1, 0.1}) {
        if(type == SixDHistogramStrType::TEXT) { _setup_from_text(s); }
    }

private:
    void
    _setup_from_text(
            Strings const & s) {
        auto lower = vector_from_str(s[0]);
        auto upper = vector_from_str(s[1]);
        auto spl = split_str_by_delimiter(s[2], " ");
        auto bin_widths = Real6();
        auto bb = BoundingBox(lower, upper);
        for(int i = 0; i < 6; i++) { bin_widths[i] = std::stod(spl[i]); }
        binner_ = SixDCoordinateBinner(bb, bin_widths);
        auto values = Real6();
        for(int i = 4; i < s.size(); i++) {
            if(s[i].length() < 5) { break; }
            spl = split_str_by_delimiter(s[i], ",");
            for(int i = 0; i < 6; i++) { values[i] = std::stod(spl[i]); }
            auto bin_index = binner_.bin_index(values);
            stored_values_[bin_index] = std::stoull(spl[6]);
        }
    }

public:
    inline
    size_t
    size() {
        return stored_values_.size();
    }

public:
    void
    add(
            Real6 const & values) {
        auto bin_index = binner_.bin_index(values);
        if(stored_values_.find(bin_index) == stored_values_.end()) {
            stored_values_[bin_index] = 0;
        }
        stored_values_[bin_index] += 1;
    }

public:

    bool
    contains(
            Real6 const & values) {
        auto bin_index = binner_.bin_index(values);
        if(stored_values_.find(bin_index) == stored_values_.end()) { return false; }
        else                                                       { return true; }

    }

public:
    void
    to_text_file(
            String const & fname) {
        std::ofstream out;
        out.open(fname);
        auto & bb = binner_.get_bounding_box();
        out << bb.lower().to_str() << std::endl;
        out << bb.upper().to_str() << std::endl;
        for(auto const & bin_width : binner_.get_bin_widths()) {
            out << bin_width << " ";
        }
        out << std::endl;

        out << "x,y,z,a,b,g,count\n";
        for(auto const & kv : stored_values_) {
            auto bin_index = kv.first;
            auto bin = binner_.bin_from_index(bin_index);
            auto values = binner_.bin_to_values(bin);
            for(auto const & v : values) {
                out << v << ",";
            }
            out << kv.second << "\n";
        }
        out.close();
    }

    void
    to_str() {
        auto s = Strings();
        s.push_back("");
        auto lower = binner_.get_bounding_box().lower();

    }

private:
    SixDCoordinateBinner binner_;
    std::map<u_int64_t, u_int64_t> stored_values_;

};


#endif //TEST_HASHING_H
