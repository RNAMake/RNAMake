//
//  xyzMatrix.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__xyzMatrix__h
#define __RNAMake__xyzMatrix__h

#include <vector>
#include <math.h>
#include <sstream>

//RNAMake Headers
#include <base/types.h>
#include <math/xyz_vector.h>

namespace math {

    template<typename T>
    class xyzMatrix {

    public:
        typedef T Value;
        typedef xyzVector<T> Vector;
        typedef std::vector<Vector> Vectors;


    public:
        template<typename> friend
        class xyzMatrix;

        friend class Transform;


        friend
        inline
        void
        transpose(
                xyzMatrix<T> const & a,
                xyzMatrix<T> & b) {

            b.xx_ = a.xx_;
            b.yy_ = a.yy_;
            b.zz_ = a.zz_;

            b.yx_ = a.xy_;
            b.xy_ = a.yx_;
            b.zx_ = a.xz_;
            b.xz_ = a.zx_;
            b.yz_ = a.zy_;
            b.zy_ = a.yz_;
        }

        friend
        inline
        xyzVector< T >
        dot_vector(
                xyzMatrix< T > const & m,
                xyzVector< T > & v) {

            auto new_v = xyzVector< T >();
            new_v.x(m.xx_ * v.get_x() + m.yx_ * v.get_y() + m.zx_ * v.get_z());
            new_v.y(m.xy_ * v.get_x() + m.yy_ * v.get_y() + m.zy_ * v.get_z());
            new_v.z(m.xz_ * v.get_x() + m.yz_ * v.get_y() + m.zz_ * v.get_z());
            return new_v;
        }


        // TODO change xyzMatrix<double> to xyzMatrix<T>
        friend
        inline
        void
        dot_vectors(
                xyzMatrix<double> const & m,
                Vectors const & v,
                Vectors & vr) {

            int i;
            for (i = 0; i < v.size(); i++) {
                dot_vector(m, v[i], vr[i]);
            }

        }


    public: //creation
        /// @brief Default constructor
        /// @note  Values are uninitialized for efficiency
        inline
        xyzMatrix() {}

        /// @brief Copy constructor
        inline
        xyzMatrix(xyzMatrix const & m) :
                xx_(m.xx_), xy_(m.xy_), xz_(m.xz_),
                yx_(m.yx_), yy_(m.yy_), yz_(m.yz_),
                zx_(m.zx_), zy_(m.zy_), zz_(m.zz_) {}


        /// @brief Copy constructor
        template<typename U>
        inline
        xyzMatrix(xyzMatrix<U> const & m) :
                xx_(m.xx_), xy_(m.xy_), xz_(m.xz_),
                yx_(m.yx_), yy_(m.yy_), yz_(m.yz_),
                zx_(m.zx_), zy_(m.zy_), zz_(m.zz_) {}


        inline
        xyzMatrix(
                const T & xx, const T & xy, const T & xz,
                const T & yx, const T & yy, const T & yz,
                const T & zx, const T & zy, const T & zz) :
                xx_(xx), xy_(xy), xz_(xz),
                yx_(yx), yy_(yy), yz_(yz),
                zx_(zx), zy_(zy), zz_(zz) {}


        /// @brief Uniform value constructor
        inline
        explicit
        xyzMatrix(Value const & t) :
                xx_(t), xy_(t), xz_(t),
                yx_(t), yy_(t), yz_(t),
                zx_(t), zy_(t), zz_(t) {}


        inline
        xyzMatrix(
                String const & s) {

            auto v = base::split_str_by_delimiter(s, " ");
            assert(v.size() > 8);
            xx_ = std::stod(v[0]);
            xy_ = std::stod(v[1]);
            xz_ = std::stod(v[2]);
            yx_ = std::stod(v[3]);
            yy_ = std::stod(v[4]);
            yz_ = std::stod(v[5]);
            zx_ = std::stod(v[6]);
            zy_ = std::stod(v[7]);
            zz_ = std::stod(v[8]);
        }

        inline
        xyzMatrix(
                json::JSON & j) :
                xx_(j[0].ToFloat()),
                xy_(j[1].ToFloat()),
                xz_(j[2].ToFloat()),
                yx_(j[3].ToFloat()),
                yy_(j[4].ToFloat()),
                yz_(j[5].ToFloat()),
                zx_(j[6].ToFloat()),
                zy_(j[7].ToFloat()),
                zz_(j[8].ToFloat()) {}

        /// @brief Destructor
        inline
        ~xyzMatrix() {}

    public:

        inline
        String const
        get_str() const {
            return std::to_string(xx_) + " " + std::to_string(xy_) + " " + \
               std::to_string(xz_) + " " + std::to_string(yx_) + " " + \
               std::to_string(yy_) + " " + std::to_string(yz_) + " " + \
               std::to_string(zx_) + " " + std::to_string(zy_) + " " + \
               std::to_string(zz_);
        }

        inline
        String const
        get_str_readable() const {
            auto s = String();
            s += "[" + std::to_string(xx_) + ", " + std::to_string(xy_) + ", " + std::to_string(xz_) + ",\n";
            s += " " + std::to_string(yx_) + ", " + std::to_string(yy_) + ", " + std::to_string(yz_) + ",\n";
            s += " " + std::to_string(zx_) + ", " + std::to_string(zy_) + ", " + std::to_string(zz_) + "]";
            return s;
        }

        inline
        json::JSON
        get_json() const { return json::Array(xx_, xy_, xz_, yx_, yy_, yz_, zx_, zy_, zz_); }


    public:

        inline
        static
        xyzMatrix
        identity() {
            return xyzMatrix(
                    Value(1), Value(0), Value(0),
                    Value(0), Value(1), Value(0),
                    Value(0), Value(0), Value(1)
            );
        }

        /// @brief Copy assignment
        template<typename U>
        inline
        xyzMatrix &
        operator=(xyzMatrix<U> const & m) {
            xx_ = m.xx_;
            xy_ = m.xy_;
            xz_ = m.xz_;
            yx_ = m.yx_;
            yy_ = m.yy_;
            yz_ = m.yz_;
            zx_ = m.zx_;
            zy_ = m.zy_;
            zz_ = m.zz_;
            return *this;
        }

        /// @brief += xyzMatrix
        template<typename U>
        inline
        xyzMatrix &
        operator+=(xyzMatrix<U> const & m) {
            xx_ += m.xx_;
            xy_ += m.xy_;
            xz_ += m.xz_;
            yx_ += m.yx_;
            yy_ += m.yy_;
            yz_ += m.yz_;
            zx_ += m.zx_;
            zy_ += m.zy_;
            zz_ += m.zz_;
            return *this;
        }

        template<typename U>
        inline
        xyzMatrix &
        operator-=(xyzMatrix<U> const & m) {
            xx_ -= m.xx_;
            xy_ -= m.xy_;
            xz_ -= m.xz_;
            yx_ -= m.yx_;
            yy_ -= m.yy_;
            yz_ -= m.yz_;
            zx_ -= m.zx_;
            zy_ -= m.zy_;
            zz_ -= m.zz_;
            return *this;
        }

    public: // Assignment: scalar
        /// @brief = Value
        inline
        xyzMatrix &
        operator=(Value const & t) {
            xx_ = xy_ = xz_ = t;
            yx_ = yy_ = yz_ = t;
            zx_ = zy_ = zz_ = t;
            return *this;
        }


        /// @brief += Value
        inline
        xyzMatrix &
        operator+=(Value const & t) {
            xx_ += t;
            xy_ += t;
            xz_ += t;
            yx_ += t;
            yy_ += t;
            yz_ += t;
            zx_ += t;
            zy_ += t;
            zz_ += t;
            return *this;
        }


        /// @brief -= Value
        inline
        xyzMatrix &
        operator-=(Value const & t) {
            xx_ -= t;
            xy_ -= t;
            xz_ -= t;
            yx_ -= t;
            yy_ -= t;
            yz_ -= t;
            zx_ -= t;
            zy_ -= t;
            zz_ -= t;
            return *this;
        }

    public: // Methods: basic mathematical
        /// @brief xyzMatrix + xyzMatrix
        friend
        inline
        xyzMatrix
        operator+(xyzMatrix const & a, xyzMatrix const & b) {
            return xyzMatrix(
                    a.xx_ + b.xx_, a.xy_ + b.xy_, a.xz_ + b.xz_,
                    a.yx_ + b.yx_, a.yy_ + b.yy_, a.yz_ + b.yz_,
                    a.zx_ + b.zx_, a.zy_ + b.zy_, a.zz_ + b.zz_
            );
        }


        /// @brief xyzMatrix + Value
        friend
        inline
        xyzMatrix
        operator+(xyzMatrix const & m, Value const & t) {
            return xyzMatrix(
                    m.xx_ + t, m.xy_ + t, m.xz_ + t,
                    m.yx_ + t, m.yy_ + t, m.yz_ + t,
                    m.zx_ + t, m.zy_ + t, m.zz_ + t
            );
        }


        /// @brief Value + xyzMatrix
        friend
        inline
        xyzMatrix
        operator+(Value const & t, xyzMatrix const & m) {
            return xyzMatrix(
                    t + m.xx_, t + m.xy_, t + m.xz_,
                    t + m.yx_, t + m.yy_, t + m.yz_,
                    t + m.zx_, t + m.zy_, t + m.zz_
            );
        }


        /// @brief xyzMatrix - xyzMatrix
        friend
        inline
        xyzMatrix
        operator-(xyzMatrix const & a, xyzMatrix const & b) {
            return xyzMatrix(
                    a.xx_ - b.xx_, a.xy_ - b.xy_, a.xz_ - b.xz_,
                    a.yx_ - b.yx_, a.yy_ - b.yy_, a.yz_ - b.yz_,
                    a.zx_ - b.zx_, a.zy_ - b.zy_, a.zz_ - b.zz_
            );
        }


        /// @brief xyzMatrix - Value
        friend
        inline
        xyzMatrix
        operator-(xyzMatrix const & m, Value const & t) {
            return xyzMatrix(
                    m.xx_ - t, m.xy_ - t, m.xz_ - t,
                    m.yx_ - t, m.yy_ - t, m.yz_ - t,
                    m.zx_ - t, m.zy_ - t, m.zz_ - t
            );
        }


        /// @brief Value - xyzMatrix
        friend
        inline
        xyzMatrix
        operator-(Value const & t, xyzMatrix const & m) {
            return xyzMatrix(
                    t - m.xx_, t - m.xy_, t - m.xz_,
                    t - m.yx_, t - m.yy_, t - m.yz_,
                    t - m.zx_, t - m.zy_, t - m.zz_
            );
        }


        /// @brief xyzMatrix * xyzMatrix
        friend
        inline
        xyzMatrix
        operator*(xyzMatrix const & a, xyzMatrix const & b) {
            return xyzMatrix(
                    // First row
                    (a.xx_ * b.xx_) + (a.xy_ * b.yx_) + (a.xz_ * b.zx_),
                    (a.xx_ * b.xy_) + (a.xy_ * b.yy_) + (a.xz_ * b.zy_),
                    (a.xx_ * b.xz_) + (a.xy_ * b.yz_) + (a.xz_ * b.zz_),

                    // Second row
                    (a.yx_ * b.xx_) + (a.yy_ * b.yx_) + (a.yz_ * b.zx_),
                    (a.yx_ * b.xy_) + (a.yy_ * b.yy_) + (a.yz_ * b.zy_),
                    (a.yx_ * b.xz_) + (a.yy_ * b.yz_) + (a.yz_ * b.zz_),

                    // Third row
                    (a.zx_ * b.xx_) + (a.zy_ * b.yx_) + (a.zz_ * b.zx_),
                    (a.zx_ * b.xy_) + (a.zy_ * b.yy_) + (a.zz_ * b.zy_),
                    (a.zx_ * b.xz_) + (a.zy_ * b.yz_) + (a.zz_ * b.zz_)
            );
        }

        /// @brief xyzMatrix * Value
        friend
        inline
        xyzMatrix
        operator*(xyzMatrix const & m, Value const & t) {
            return xyzMatrix(
                    m.xx_ * t, m.xy_ * t, m.xz_ * t,
                    m.yx_ * t, m.yy_ * t, m.yz_ * t,
                    m.zx_ * t, m.zy_ * t, m.zz_ * t
            );
        }


        /// @brief Value * xyzMatrix
        friend
        inline
        xyzMatrix
        operator*(Value const & t, xyzMatrix const & m) {
            return xyzMatrix(
                    t * m.xx_, t * m.xy_, t * m.xz_,
                    t * m.yx_, t * m.yy_, t * m.yz_,
                    t * m.zx_, t * m.zy_, t * m.zz_
            );
        }


        /// @brief xyzMatrix / Value
        friend
        inline
        xyzMatrix
        operator/(xyzMatrix const & m, Value const & t) {
            assert(t != Value(0));
            Value const inv_t(Value(1) / t);
            return xyzMatrix(
                    m.xx_ * inv_t, m.xy_ * inv_t, m.xz_ * inv_t,
                    m.yx_ * inv_t, m.yy_ * inv_t, m.yz_ * inv_t,
                    m.zx_ * inv_t, m.zy_ * inv_t, m.zz_ * inv_t
            );
        }

    public:
        /// @brief Transpose
        inline
        xyzMatrix &
        transpose() {
            Value temp = xy_;
            xy_ = yx_;
            yx_ = temp;

            temp = xz_;
            xz_ = zx_;
            zx_ = temp;

            temp = yz_;
            yz_ = zy_;
            zy_ = temp;

            return *this;
        }

        inline
        double const
        difference(
                xyzMatrix<T> const & b) const {

            double dist = 0.0f;
            dist += std::abs(xx_ - b.xx_);
            dist += std::abs(xy_ - b.xy_);
            dist += std::abs(xz_ - b.xz_);
            dist += std::abs(yx_ - b.yx_);
            dist += std::abs(yy_ - b.yy_);
            dist += std::abs(yz_ - b.yz_);
            dist += std::abs(zx_ - b.zx_);
            dist += std::abs(zy_ - b.zy_);
            dist += std::abs(zz_ - b.zz_);

            return dist;

        }

        inline const
        xyzMatrix
        get_flip_orientation() const {
            return xyzMatrix(
                    xx_, xy_, xz_,
                    -yx_, -yy_, -yz_,
                    -zx_, -zy_, -zz_);

        }

        inline const
        xyzMatrix
        get_unitarize() const {


            auto m = xyzMatrix(xx_, xy_, xz_,
                               yx_, yy_, yz_,
                               zx_, zy_, zz_);

            //R[0] /= math.sqrt(R[0].dot(R[0]))
            double dot = sqrt(xx_ * xx_ + xy_ * xy_ + xz_ * xz_);
            m.xx_ /= dot;
            m.xy_ /= dot;
            m.xz_ /= dot;
            //R[1] -= R[1].dot(R[0]) * R[0]
            dot = yx_ * m.xx_ + yy_ * m.xy_ + yz_ * m.xz_;
            m.yx_ -= dot * m.xx_;
            m.yy_ -= dot * m.xy_;
            m.yz_ -= dot * m.xz_;
            //R[1] /= math.sqrt(R[1].dot(R[1]))
            dot = sqrt(m.yx_ * m.yx_ + m.yy_ * m.yy_ + m.yz_ * m.yz_);
            m.yx_ /= dot;
            m.yy_ /= dot;
            m.yz_ /= dot;
            //R[2] -= R[2].dot(R[0]) * R[0]
            dot = m.zx_ * m.xx_ + m.zy_ * m.xy_ + m.zz_ * m.xz_;
            m.zx_ -= dot * m.xx_;
            m.zy_ -= dot * m.xy_;
            m.zz_ -= dot * m.xz_;
            //R[2] -= R[2].dot(R[1]) * R[1]
            dot = m.zx_ * m.yx_ + m.zy_ * m.yy_ + m.zz_ * m.yz_;
            m.zx_ -= dot * m.yx_;
            m.zy_ -= dot * m.yy_;
            m.zz_ -= dot * m.yz_;
            //R[2] /= math.sqrt(R[2].dot(R[2]))
            dot = sqrt(m.zx_ * m.zx_ + m.zy_ * m.zy_ + m.zz_ * m.zz_);
            m.zx_ /= dot;
            m.zy_ /= dot;
            m.zz_ /= dot;


            return m;
        }

        inline
        void
        unitarize() {
            double dot = sqrt(xx_ * xx_ + xy_ * xy_ + xz_ * xz_);
            xx_ /= dot;
            xy_ /= dot;
            xz_ /= dot;
            //R[1] -= R[1].dot(R[0]) * R[0]
            dot = yx_ * xx_ + yy_ * xy_ + yz_ * xz_;
            yx_ -= dot * xx_;
            yy_ -= dot * xy_;
            yz_ -= dot * xz_;
            //R[1] /= math.sqrt(R[1].dot(R[1]))
            dot = sqrt(yx_ * yx_ + yy_ * yy_ + yz_ * yz_);
            yx_ /= dot;
            yy_ /= dot;
            yz_ /= dot;
            //R[2] -= R[2].dot(R[0]) * R[0]
            dot = zx_ * xx_ + zy_ * xy_ + zz_ * xz_;
            zx_ -= dot * xx_;
            zy_ -= dot * xy_;
            zz_ -= dot * xz_;
            //R[2] -= R[2].dot(R[1]) * R[1]
            dot = zx_ * yx_ + zy_ * yy_ + zz_ * yz_;
            zx_ -= dot * yx_;
            zy_ -= dot * yy_;
            zz_ -= dot * yz_;
            //R[2] /= math.sqrt(R[2].dot(R[2]))
            dot = sqrt(zx_ * zx_ + zy_ * zy_ + zz_ * zz_);
            zx_ /= dot;
            zy_ /= dot;
            zz_ /= dot;
        }


    public: // Properties: scalars


        /// @brief Value xx const
        inline
        Value const &
        get_xx() const { return xx_; }

        /// @brief Value xy const
        inline
        Value const &
        get_xy() const { return xy_; }

        /// @brief Value xz const
        inline
        Value const &
        get_xz() const { return xz_; }

        /// @brief Value yx const
        inline
        Value const &
        get_yx() const { return yx_; }

        /// @brief Value yy const
        inline
        Value const &
        get_yy() const { return yy_; }

        /// @brief Value yz const
        inline
        Value const &
        get_yz() const { return yz_; }

        /// @brief Value zx const
        inline
        Value const &
        get_zx() const { return zx_; }

        /// @brief Value zy const
        inline
        Value const &
        get_zy() const { return zy_; }

        /// @brief Value zz const
        inline
        Value const &
        get_zz() const { return zz_; }


    public: // Properties: value assignment


        /// @brief xx assignment
        inline
        void
        set_xx(Value const & xx_a) {
            xx_ = xx_a;
        }


        /// @brief xy assignment
        inline
        void
        set_xy(Value const & xy_a) {
            xy_ = xy_a;
        }


        /// @brief xz assignment
        inline
        void
        set_xz(Value const & xz_a) {
            xz_ = xz_a;
        }


        /// @brief yx assignment
        inline
        void
        set_yx(Value const & yx_a) {
            yx_ = yx_a;
        }


        /// @brief yy assignment
        inline
        void
        set_yy(Value const & yy_a) {
            yy_ = yy_a;
        }


        /// @brief yz assignment
        inline
        void
        set_yz(Value const & yz_a) {
            yz_ = yz_a;
        }


        /// @brief zx assignment
        inline
        void
        set_zx(Value const & zx_a) {
            zx_ = zx_a;
        }


        /// @brief zy assignment
        inline
        void
        set_zy(Value const & zy_a) {
            zy_ = zy_a;
        }


        /// @brief zz assignment
        inline
        void
        set_zz(Value const & zz_a) {
            zz_ = zz_a;
        }

        inline
        xyzMatrix<T>
        get_transposed() const {
            return xyzMatrix(
                    xx_, yx_, zx_,
                    xy_, yy_, zy_,
                    xz_, yz_, zz_
            );
        }


    private:
        Value xx_, xy_, xz_;
        Value yx_, yy_, yz_;
        Value zx_, zy_, zz_;


    };

    typedef xyzMatrix<double> Matrix;
    typedef std::vector<Matrix> Matrices;


    //TODO Remove old code and ask joe about transform_1
    template< typename T >
    inline
    xyzVector< T >
    operator *( xyzMatrix< T > const & m, xyzVector< T > const & v )
    {
        return xyzVector< T >(
                m.xx() * v.x() + m.get_xy() * v.y() + m.get_xz() * v.get_z(),
                m.get_yx() * v.x() + m.get_yy() * v.y() + m.get_yz() * v.get_z(),
                m.get_zx() * v.x() + m.get_zy() * v.y() + m.get_zz() * v.get_z());
    }

    template< typename T >
    inline
    xyzMatrix< T >
    transform_1( xyzMatrix < T > const & m ) {

        return xyzMatrix< T >(
                m.get_xx(),  m.get_xy(),  m.get_xz(),
                -m.get_yx(), -m.get_yy(), -m.get_yz(),
                -m.get_zx(), -m.get_zy(), -m.get_zz());


    }

    inline
    void
    dot_vector(
            Matrix const & m,
            Vector const & v,
            Vector & vr) {

        vr.set_x(m.get_xx() * v.get_x() + m.get_yx() * v.get_y() + m.get_zx() * v.get_z());
        vr.set_y(m.get_xy() * v.get_x() + m.get_yy() * v.get_y() + m.get_zy() * v.get_z());
        vr.set_z(m.get_xz() * v.get_x() + m.get_yz() * v.get_y() + m.get_zz() * v.get_z());

    }

    inline
    Vector
    dot_vector(
            Matrix const & m,
            Vector const & v) {

        auto vr = Vector();
        vr.set_x(m.get_xx() * v.get_x() + m.get_yx() * v.get_y() + m.get_zx() * v.get_z());
        vr.set_y(m.get_xy() * v.get_x() + m.get_yy() * v.get_y() + m.get_zy() * v.get_z());
        vr.set_z(m.get_xz() * v.get_x() + m.get_yz() * v.get_y() + m.get_zz() * v.get_z());
        return vr;
    }

    inline
    void
    dot(
            Matrix const & a,
            Matrix const & b,
            Matrix & c) {

        c.set_xx(a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() + a.get_xz() * b.get_zx());
        c.set_xy(a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() + a.get_xz() * b.get_zy());
        c.set_xz(a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() + a.get_xz() * b.get_zz());
        c.set_yx(a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() + a.get_yz() * b.get_zx());
        c.set_yy(a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() + a.get_yz() * b.get_zy());
        c.set_yz(a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() + a.get_yz() * b.get_zz());
        c.set_zx(a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() + a.get_zz() * b.get_zx());
        c.set_zy(a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() + a.get_zz() * b.get_zy());
        c.set_zz(a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() + a.get_zz() * b.get_zz());
    }

    inline
    Matrix
    dot(
            Matrix const & a,
            Matrix const & b) {
        return Matrix(
                a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() + a.get_xz() * b.get_zx(),
                a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() + a.get_xz() * b.get_zy(),
                a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() + a.get_xz() * b.get_zz(),
                a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() + a.get_yz() * b.get_zx(),
                a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() + a.get_yz() * b.get_zy(),
                a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() + a.get_yz() * b.get_zz(),
                a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() + a.get_zz() * b.get_zx(),
                a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() + a.get_zz() * b.get_zy(),
                a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() + a.get_zz() * b.get_zz());
    }
  template< typename T >
  std::ostream &
  operator <<( std::ostream & stream, xyzMatrix< T > const & v ) {
      stream << "(" << v.get_xx() << ", " << v.get_xy() << ", " << v.get_xz() << ")" << std::endl;
      stream << "(" << v.get_yx() << ", " << v.get_yy() << ", " << v.get_yz() << ")" << std::endl;
      stream << "(" << v.get_zx() << ", " << v.get_zy() << ", " << v.get_zz() << ")" << std::endl;
      return stream;
  }

  inline
  String
  matrix_to_str(Matrix const & m) {
      std::stringstream ss;
      ss << m.get_xx() << " " << m.get_xy() << " " << m.get_xz() << " ";
      ss << m.get_yx() << " " << m.get_yy() << " " << m.get_yz() << " ";
      ss << m.get_zx() << " " << m.get_zy() << " " << m.get_zz() << " ";
      return ss.str();
  }


}


#endif