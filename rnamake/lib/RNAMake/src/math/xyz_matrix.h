//
//  xyzMatrix.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef REDESIGNC_xyzMatrix_h
#define REDESIGNC_xyzMatrix_h

#include <vector>
#include <math.h>
#include <sstream>

//RNAMake Headers
#include "base/types.h"
#include "math/xyz_vector.h"


template< typename T >
class xyzMatrix {
	
public:
	typedef  T					 Value;
	typedef  xyzVector< T >		 Vector;
	typedef  std::vector<Vector> Vectors;


private:
	template< typename > friend class xyzMatrix;
	friend class Transform;

	
	friend
	inline
	void
	transpose(
		xyzMatrix< T > const & a,
		xyzMatrix< T >  & b ) {
		
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
	void
	dot_vector(
		xyzMatrix< T > const & m,
		xyzVector< T > const & v,
		xyzVector< T > & vr) {
		
		vr.x ( m.xx_ * v.x() + m.yx_ * v.y() + m.zx_ * v.z()) ;
		vr.y ( m.xy_ * v.x() + m.yy_ * v.y() + m.zy_ * v.z());
		vr.z ( m.xz_ * v.x() + m.yz_ * v.y() + m.zz_ * v.z());
		
	}
	
	friend
	inline
	void
	dot_vectors(
			   xyzMatrix< T > const & m,
			   Vectors const & v,
			   Vectors & vr) {
		
		int i;
		for(i = 0; i < v.size(); i++) {
			dot_vector(m,v[i],vr[i]);
		}
		
	}
	
	friend
	inline
	void
	dot(
		xyzMatrix< T > const & a,
		xyzMatrix< T > const & b,
		xyzMatrix< T > & c)
	{
		
		c.xx_ = a.xx_*b.xx_ + a.xy_*b.yx_ + a.xz_*b.zx_;
		c.xy_ = a.xx_*b.xy_ + a.xy_*b.yy_ + a.xz_*b.zy_;
		c.xz_ = a.xx_*b.xz_ + a.xy_*b.yz_ + a.xz_*b.zz_;
		c.yx_ = a.yx_*b.xx_ + a.yy_*b.yx_ + a.yz_*b.zx_;
		c.yy_ = a.yx_*b.xy_ + a.yy_*b.yy_ + a.yz_*b.zy_;
		c.yz_ = a.yx_*b.xz_ + a.yy_*b.yz_ + a.yz_*b.zz_;
		c.zx_ = a.zx_*b.xx_ + a.zy_*b.yx_ + a.zz_*b.zx_;
		c.zy_ = a.zx_*b.xy_ + a.zy_*b.yy_ + a.zz_*b.zy_;
		c.zz_ = a.zx_*b.xz_ + a.zy_*b.yz_ + a.zz_*b.zz_;
	}
	

public: //creation
	/// @brief Default constructor
	/// @note  Values are uninitialized for efficiency
	inline
	xyzMatrix()
	{}
	
	/// @brief Copy constructor
	inline
	xyzMatrix( xyzMatrix const & m ) :
	xx_( m.xx_ ), xy_( m.xy_ ), xz_( m.xz_ ),
	yx_( m.yx_ ), yy_( m.yy_ ), yz_( m.yz_ ),
	zx_( m.zx_ ), zy_( m.zy_ ), zz_( m.zz_ )
	{}
	
	
	/// @brief Copy constructor
	template< typename U >
	inline
	xyzMatrix( xyzMatrix< U > const & m ) :
	xx_( m.xx_ ), xy_( m.xy_ ), xz_( m.xz_ ),
	yx_( m.yx_ ), yy_( m.yy_ ), yz_( m.yz_ ),
	zx_( m.zx_ ), zy_( m.zy_ ), zz_( m.zz_ )
	{}
	
	
	inline
	xyzMatrix(
		const T & xx, const T & xy, const T & xz,
		const T & yx, const T & yy, const T & yz,
		const T & zx, const T & zy, const T & zz):
	xx_( xx ), xy_( xy ), xz_( xz ),
	yx_( yx ), yy_( yy ), yz_( yz ),
	zx_( zx ), zy_( zy ), zz_( zz )
	{}
	
	
	/// @brief Uniform value constructor
	inline
	explicit
	xyzMatrix( Value const & t ) :
	xx_( t ), xy_( t ), xz_( t ),
	yx_( t ), yy_( t ), yz_( t ),
	zx_( t ), zy_( t ), zz_( t )
	{}
	
	/// @brief Destructor
	inline
	~xyzMatrix()
	{}

public:
	
	inline
	void
	row(const int i,
		Vector const & v) {
		if(i == 0) { xx_ = v.x(); xy_ = v.y(); xz_ = v.z(); }
		if(i == 1) { yx_ = v.x(); yy_ = v.y(); yz_ = v.z(); }
		if(i == 2) { zx_ = v.x(); zy_ = v.y(); zz_ = v.z(); }

	}
	
	inline
	void
	row(const int i,
		std::vector< T > const & v) {
		if(i == 0) { xx_ = v[0]; xy_ = v[1]; xz_ = v[2]; }
		if(i == 1) { yx_ = v[0]; yy_ = v[1]; yz_ = v[2]; }
		if(i == 2) { zx_ = v[0]; zy_ = v[1]; zz_ = v[2]; }
		
	}
	
	inline
	static
	xyzMatrix
	identity() {
		return xyzMatrix(
			 Value( 1 ), Value( 0 ), Value( 0 ),
			 Value( 0 ), Value( 1 ), Value( 0 ),
			 Value( 0 ), Value( 0 ), Value( 1 )
			 );
	}
	
	/// @brief Copy assignment
	template< typename U >
	inline
	xyzMatrix &
	operator =( xyzMatrix< U > const & m )
	{
		xx_ = m.xx_; xy_ = m.xy_; xz_ = m.xz_;
		yx_ = m.yx_; yy_ = m.yy_; yz_ = m.yz_;
		zx_ = m.zx_; zy_ = m.zy_; zz_ = m.zz_;
		return *this;
	}
	
	/// @brief += xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	operator +=( xyzMatrix< U > const & m )
	{
		xx_ += m.xx_; xy_ += m.xy_; xz_ += m.xz_;
		yx_ += m.yx_; yy_ += m.yy_; yz_ += m.yz_;
		zx_ += m.zx_; zy_ += m.zy_; zz_ += m.zz_;
		return *this;
	}
	
	template< typename U >
	inline
	xyzMatrix &
	operator -=( xyzMatrix< U > const & m )
	{
		xx_ -= m.xx_; xy_ -= m.xy_; xz_ -= m.xz_;
		yx_ -= m.yx_; yy_ -= m.yy_; yz_ -= m.yz_;
		zx_ -= m.zx_; zy_ -= m.zy_; zz_ -= m.zz_;
		return *this;
	}
	
public: // Assignment: scalar
	/// @brief = Value
	inline
	xyzMatrix &
	operator =( Value const & t )
	{
		xx_ = xy_ = xz_ = t;
		yx_ = yy_ = yz_ = t;
		zx_ = zy_ = zz_ = t;
		return *this;
	}
	
	
	/// @brief += Value
	inline
	xyzMatrix &
	operator +=( Value const & t )
	{
		xx_ += t; xy_ += t; xz_ += t;
		yx_ += t; yy_ += t; yz_ += t;
		zx_ += t; zy_ += t; zz_ += t;
		return *this;
	}
	
	
	/// @brief -= Value
	inline
	xyzMatrix &
	operator -=( Value const & t )
	{
		xx_ -= t; xy_ -= t; xz_ -= t;
		yx_ -= t; yy_ -= t; yz_ -= t;
		zx_ -= t; zy_ -= t; zz_ -= t;
		return *this;
	}

public: // Methods: basic mathematical
	/// @brief xyzMatrix + xyzMatrix
	friend
	inline
	xyzMatrix
	operator +( xyzMatrix const & a, xyzMatrix const & b )
	{
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
	operator +( xyzMatrix const & m, Value const & t )
	{
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
	operator +( Value const & t, xyzMatrix const & m )
	{
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
	operator -( xyzMatrix const & a, xyzMatrix const & b )
	{
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
	operator -( xyzMatrix const & m, Value const & t )
	{
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
	operator -( Value const & t, xyzMatrix const & m )
	{
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
	operator *( xyzMatrix const & a, xyzMatrix const & b )
	{
		return xyzMatrix(
			 // First row
			 ( a.xx_ * b.xx_ ) + ( a.xy_ * b.yx_ ) + ( a.xz_ * b.zx_ ),
			 ( a.xx_ * b.xy_ ) + ( a.xy_ * b.yy_ ) + ( a.xz_ * b.zy_ ),
			 ( a.xx_ * b.xz_ ) + ( a.xy_ * b.yz_ ) + ( a.xz_ * b.zz_ ),
			 
			 // Second row
			 ( a.yx_ * b.xx_ ) + ( a.yy_ * b.yx_ ) + ( a.yz_ * b.zx_ ),
			 ( a.yx_ * b.xy_ ) + ( a.yy_ * b.yy_ ) + ( a.yz_ * b.zy_ ),
			 ( a.yx_ * b.xz_ ) + ( a.yy_ * b.yz_ ) + ( a.yz_ * b.zz_ ),
			 
			 // Third row
			 ( a.zx_ * b.xx_ ) + ( a.zy_ * b.yx_ ) + ( a.zz_ * b.zx_ ),
			 ( a.zx_ * b.xy_ ) + ( a.zy_ * b.yy_ ) + ( a.zz_ * b.zy_ ),
			 ( a.zx_ * b.xz_ ) + ( a.zy_ * b.yz_ ) + ( a.zz_ * b.zz_ )
			 );
	}
	
	/// @brief xyzMatrix * Value
	friend
	inline
	xyzMatrix
	operator *( xyzMatrix const & m, Value const & t )
	{
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
	operator *( Value const & t, xyzMatrix const & m )
	{
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
	operator /( xyzMatrix const & m, Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value( 1 ) / t );
		return xyzMatrix(
						 m.xx_ * inv_t, m.xy_ * inv_t, m.xz_ * inv_t,
						 m.yx_ * inv_t, m.yy_ * inv_t, m.yz_ * inv_t,
						 m.zx_ * inv_t, m.zy_ * inv_t, m.zz_ * inv_t
						 );
	}

	/// @brief Transpose
	inline
	xyzMatrix &
	transpose()
	{
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
	float const
	difference(
			 xyzMatrix < T > const & b) const {
		
		float dist = 0.0f;
		dist += fabs(xx_ - b.xx_);
		dist += fabs(xy_ - b.xy_);
		dist += fabs(xz_ - b.xz_);
		dist += fabs(yx_ - b.yx_);
		dist += fabs(yy_ - b.yy_);
		dist += fabs(yz_ - b.yz_);
		dist += fabs(zx_ - b.zx_);
		dist += fabs(zy_ - b.zy_);
		dist += fabs(zz_ - b.zz_);
		
		return dist;
		
	}

	inline const
	xyzMatrix
	get_flip_orientation() const{
		return xyzMatrix(
			 xx_,  xy_,  xz_,
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
        double dot = sqrt(xx_*xx_ + xy_*xy_ + xz_*xz_);
        m.xx_ /= dot; m.xy_ /= dot; m.xz_ /= dot;
        //R[1] -= R[1].dot(R[0]) * R[0]
        dot = yx_*m.xx_ + yy_*m.xy_ + yz_*m.xz_;
        m.yx_ -= dot*m.xx_; m.yy_ -= dot*m.xy_; m.yz_ -= dot*m.xz_;
        //R[1] /= math.sqrt(R[1].dot(R[1]))
        dot = sqrt(m.yx_*m.yx_ + m.yy_*m.yy_ + m.yz_*m.yz_);
        m.yx_ /= dot; m.yy_ /= dot; m.yz_ /= dot;
        //R[2] -= R[2].dot(R[0]) * R[0]
        dot = m.zx_*m.xx_ + m.zy_*m.xy_ + m.zz_*m.xz_;
        m.zx_ -= dot*m.xx_; m.zy_ -= dot*m.xy_; m.zz_ -= dot*m.xz_;
        //R[2] -= R[2].dot(R[1]) * R[1]
        dot = m.zx_*m.yx_ + m.zy_*m.yy_ + m.zz_*m.yz_;
        m.zx_ -= dot*m.yx_; m.zy_ -= dot*m.yy_; m.zz_ -= dot*m.yz_;
        //R[2] /= math.sqrt(R[2].dot(R[2]))
        dot = sqrt(m.zx_*m.zx_ + m.zy_*m.zy_ + m.zz_*m.zz_);
        m.zx_ /= dot; m.zy_ /= dot; m.zz_ /= dot;

        
        return m;
    }
	
public: // Properties: scalars
	
	
	/// @brief Value xx const
	inline
	Value const &
	xx() const
	{
		return xx_;
	}
	
	
	/// @brief Value xx
	inline
	Value &
	xx()
	{
		return xx_;
	}
	
	
	/// @brief Value xy const
	inline
	Value const &
	xy() const
	{
		return xy_;
	}
	
	
	/// @brief Value xy
	inline
	Value &
	xy()
	{
		return xy_;
	}
	
	
	/// @brief Value xz const
	inline
	Value const &
	xz() const
	{
		return xz_;
	}
	
	
	/// @brief Value xz
	inline
	Value &
	xz()
	{
		return xz_;
	}
	
	
	/// @brief Value yx const
	inline
	Value const &
	yx() const
	{
		return yx_;
	}
	
	
	/// @brief Value yx
	inline
	Value &
	yx()
	{
		return yx_;
	}
	
	
	/// @brief Value yy const
	inline
	Value const &
	yy() const
	{
		return yy_;
	}
	
	
	/// @brief Value yy
	inline
	Value &
	yy()
	{
		return yy_;
	}
	
	
	/// @brief Value yz const
	inline
	Value const &
	yz() const
	{
		return yz_;
	}
	
	
	/// @brief Value yz
	inline
	Value &
	yz()
	{
		return yz_;
	}
	
	
	/// @brief Value zx const
	inline
	Value const &
	zx() const
	{
		return zx_;
	}
	
	
	/// @brief Value zx
	inline
	Value &
	zx()
	{
		return zx_;
	}
	
	
	/// @brief Value zy const
	inline
	Value const &
	zy() const
	{
		return zy_;
	}
	
	
	/// @brief Value zy
	inline
	Value &
	zy()
	{
		return zy_;
	}
	
	
	/// @brief Value zz const
	inline
	Value const &
	zz() const
	{
		return zz_;
	}
	
	
	/// @brief Value zz
	inline
	Value &
	zz()
	{
		return zz_;
	}
	
public: // Properties: value assignment
	
	
	/// @brief xx assignment
	inline
	void
	xx( Value const & xx_a )
	{
		xx_ = xx_a;
	}
	
	
	/// @brief xy assignment
	inline
	void
	xy( Value const & xy_a )
	{
		xy_ = xy_a;
	}
	
	
	/// @brief xz assignment
	inline
	void
	xz( Value const & xz_a )
	{
		xz_ = xz_a;
	}
	
	
	/// @brief yx assignment
	inline
	void
	yx( Value const & yx_a )
	{
		yx_ = yx_a;
	}
	
	
	/// @brief yy assignment
	inline
	void
	yy( Value const & yy_a )
	{
		yy_ = yy_a;
	}
	
	
	/// @brief yz assignment
	inline
	void
	yz( Value const & yz_a )
	{
		yz_ = yz_a;
	}
	
	
	/// @brief zx assignment
	inline
	void
	zx( Value const & zx_a )
	{
		zx_ = zx_a;
	}
	
	
	/// @brief zy assignment
	inline
	void
	zy( Value const & zy_a )
	{
		zy_ = zy_a;
	}
	
	
	/// @brief zz assignment
	inline
	void
	zz( Value const & zz_a )
	{
		zz_ = zz_a;
	}
	
	inline
	xyzMatrix<T>
	transposed()
	{
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

inline
const
Matrix
matrix_from_str(
	std::string const & s) {
	
    std::vector<std::string> values = split_str_by_delimiter(s," ");
	std::vector<double> point;
	Matrix m(0);
	int j = 0;
	for (std::vector<std::string>::iterator i = values.begin();
		 i != values.end(); ++i) {
		
		point.push_back(atof(i->c_str()));
		if (point.size() == 3) {
			m.row(j,point);
			point = std::vector<double>();
			j += 1;

		}
	}
	return m;
}

template< typename T >
inline
xyzVector< T >
operator *( xyzMatrix< T > const & m, xyzVector< T > const & v )
{
	return xyzVector< T >(
		m.xx() * v.x() + m.xy() * v.y() + m.xz() * v.z(),
		m.yx() * v.x() + m.yy() * v.y() + m.yz() * v.z(),
		m.zx() * v.x() + m.zy() * v.y() + m.zz() * v.z());
}

template< typename T >
inline
xyzMatrix< T >
transform_1( xyzMatrix < T > const & m ) {
	
	return xyzMatrix< T >(
		m.xx(),  m.xy(),  m.xz(),
	   -m.yx(), -m.yy(), -m.yz(),
	   -m.zx(), -m.zy(), -m.zz());
	
	
}



template< typename T >
std::ostream &
operator <<( std::ostream & stream, xyzMatrix< T > const & v ) {
	stream << "(" << v.xx() << ", " << v.xy() << ", " << v.xz() << ")" << std::endl;
	stream << "(" << v.yx() << ", " << v.yy() << ", " << v.yz() << ")" << std::endl;
	stream << "(" << v.zx() << ", " << v.zy() << ", " << v.zz() << ")" << std::endl;
	return stream;
}

inline
String
matrix_to_str(Matrix const & m) {
    std::stringstream ss;
    ss << m.xx() << " " << m.xy() << " " << m.xz() << " ";
    ss << m.yx() << " " << m.yy() << " " << m.yz() << " ";
    ss << m.zx() << " " << m.zy() << " " << m.zz() << " ";
    return ss.str();
}




#endif
