//
//  Transform.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 10/1/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Transform__
#define __REDESIGNC__Transform__

#include <iostream>

//RNAMake Headers
#include "math/xyz_vector.h"
#include "math/xyz_matrix.h"


class Transform {
public:
    template< typename > friend class xyzMatrix;
    template< typename > friend class xyzVector;
    
public:
		
	friend
	inline
	void
	dot(
		Matrix const & a,
		Matrix const & b,
		Transform & c)
	{
		
		c.xx_ = a.xx()*b.xx() + a.xy()*b.yx() + a.xz()*b.zx();
		c.xy_ = a.xx()*b.xy() + a.xy()*b.yy() + a.xz()*b.zy();
		c.xz_ = a.xx()*b.xz() + a.xy()*b.yz() + a.xz()*b.zz();
		c.yx_ = a.yx()*b.xx() + a.yy()*b.yx() + a.yz()*b.zx();
		c.yy_ = a.yx()*b.xy() + a.yy()*b.yy() + a.yz()*b.zy();
		c.yz_ = a.yx()*b.xz() + a.yy()*b.yz() + a.yz()*b.zz();
		c.zx_ = a.zx()*b.xx() + a.zy()*b.yx() + a.zz()*b.zx();
		c.zy_ = a.zx()*b.xy() + a.zy()*b.yy() + a.zz()*b.zy();
		c.zz_ = a.zx()*b.xz() + a.zy()*b.yz() + a.zz()*b.zz();
	}
    
	
public: //creation

	/// @brief Default constructor
	inline
	Transform():
	xx_( 1.0 ),
	yx_( 0.0 ),
	zx_( 0.0 ),
	px_( 0.0 ),
	xy_( 0.0 ),
	yy_( 1.0 ),
	zy_( 0.0 ),
	py_( 0.0 ),
	xz_( 0.0 ),
	yz_( 0.0 ),
	zz_( 1.0 ),
	pz_( 0.0 )
	{}
    
    inline
    Transform(Matrix const & r, Point const & t ) {
        rotation(r);
        translation(t);
    }
	
public:
	//Accessors
	
	float xx() const { return xx_; }
	float xy() const { return xy_; }
	float xz() const { return xz_; }
	float yx() const { return yx_; }
	float yy() const { return yy_; }
	float yz() const { return yz_; }
	float zx() const { return zx_; }
	float zy() const { return zy_; }
	float zz() const { return zz_; }
	float px() const { return px_; }
	float py() const { return py_; }
	float pz() const { return pz_; }
	float xx()		 { return xx_; }
	float xy()       { return xy_; }
	float xz()       { return xz_; }
	float yx()       { return yx_; }
	float yy()       { return yy_; }
	float yz()       { return yz_; }
	float zx()       { return zx_; }
	float zy()       { return zy_; }
	float zz()       { return zz_; }
	float px()       { return px_; }
	float py()       { return py_; }
	float pz()       { return pz_; }
	
	Vector
	xaxis() const {
		return Vector( xx_, xy_, xz_ );
	}
	
	Vector
	yaxis() const {
		return Vector( yx_, yy_, yz_ );
	}
	
	Vector
	zaxis() const {
		return Vector( zx_, zy_, zz_ );
	}
	
	Vector
	translation() const {
		return Vector( px_, py_, pz_ );
	}
	
public:
    
    inline
    Matrix
    rotation() const {
        
        return Matrix(xx_, xy_, xz_,
                      yx_, yy_, yz_,
                      zx_, zy_, zz_);
    }
    
	inline
	void
	rotation(
		Matrix const & m) {
		xx_ = m.xx(); xy_ = m.xy(); xz_ = m.xz();
		yx_ = m.yx(); yy_ = m.yy(); yz_ = m.yz();
		zx_ = m.zx(); zy_ = m.zy(); zz_ = m.zz();
	}
	
	inline
	void
	translation(
		Vector const & v) {
		px_ = v.x(); py_ = v.y(); pz_ = v.z();
	}
	
	
	
private:
	float xx_, yx_, zx_, px_;
	float xy_, yy_, zy_, py_;
	float xz_, yz_, zz_, pz_;

	
};




#endif /* defined(__REDESIGNC__Transform__) */
