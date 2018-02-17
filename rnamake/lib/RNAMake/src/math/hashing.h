//
// Created by Joseph Yesselman on 2/16/18.
//

#ifndef TEST_HASHING_H
#define TEST_HASHING_H

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
            Real6 const & bin_widths,
            Size3 const & euler_offsets):
            bounding_box_(bounding_box),
            bin_widths_(bin_widths),
            dimsizes_(Size6()),
            dimprods_(Size6()) {

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

        for (int ii = 3; ii <= 4; ++ii) {
            dimsizes_[ii] = static_cast<size_t> ( 360.0 / bin_widths_[ii] );
        }
        dimsizes_[5] = static_cast<size_t> ( 180.0 / bin_widths_[5] );
        if (dimsizes_[5] == 0) { dimsizes_[5] = 1; }

        /// Add an extra bin so that values near 180 ( wi binwidth/2) and near 0 ( wi binwidth/2 ) can be joined.
        if (euler_offsets[2] != 0) {
            ++dimsizes_[5];
        }

        dimprods_[5] = 1;
        for (int ii = 4; ii >= 0; --ii) {
            dimprods_[ii] = dimprods_[ii + 1] * dimsizes_[ii + 1];
        }

        for (int ii = 0; ii < 3; ++ii) {
            if (euler_offsets[ii] != 0) {
                euler_offsets_[ii] = bin_widths_[ii + 3] / 2;
            } else {
                euler_offsets_[ii] = 0;
            }
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
            center[ii + 3] = euler_offsets_[ii] + bin[ii + 3] * bin_widths_[ii + 3] + halfbin_widths_[ii + 3];
        }
        return center;
    }


private:
    Real3
    _wrap_euler_angles(
            Real6 const & values) const {

        Real3 euler;

        if ( values[ 5 ] > euler_offsets_[ 2 ] ) {
            euler[ 2 ] = values[ 5 ] + euler_offsets_[ 2 ];

        } else {
            euler[ 2 ] = values[ 5 ];
        }

        if ( ( values[ 5 ] < euler_offsets_[ 2 ] || values[ 5 ] >= 180.0 - euler_offsets_[ 2 ] ) &&
             ( values[ 3 ] - euler_offsets_[ 0 ] > 180.0 ) ) {
            /// WRAP if phi > 180 to the region of negative theta rotation.
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
        }
        return euler;


    }

private:
    BoundingBox bounding_box_;
    Size6 dimsizes_;
    Size6 dimprods_;
    Real6 bin_widths_;
    Real6 halfbin_widths_;
    Size3 euler_offsets_;


};


#endif //TEST_HASHING_H
