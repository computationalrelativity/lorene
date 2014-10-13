/*
 *  Definition of Lorene class Diff.
 *
 */

/*
 *   Copyright (c) 2005 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef __DIFF_H_ 
#define __DIFF_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:52:33  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2006/04/10 15:20:51  j_novak
 * Operators dsdx and sx can now be used in the nucleus.
 *
 * Revision 1.2  2005/01/11 15:16:58  j_novak
 * More Diff operators.
 *
 * Revision 1.1  2005/01/10 16:39:21  j_novak
 * New class for 1D mono-domain differential operators.
 *
 *
 *
 * $Header$
 *
 */

#include "type_parite.h"
#include "matrice.h"

namespace Lorene {
/**
 * Base (abstract) class for 1D spectral differential operators in one domain.
 * \ingroup (ellip)
 *
 * This class is intended as a base class for several classes of elementary 
 * differential operators in term of \f$\xi\f$, in a given domain. 
 * Their main purpose is to compute the matrix of the elementary operator,
 * for the given spectral base and number of coefficients. Some of the
 * operators are not defined in the shells (division by \f$\xi\f$).
 * 
 */
class Diff {

    // Data : 
    // -----
 public:
    /// Maximal number of matrices stored per base
    static const int max_points = 50 ; 

 protected:
    int base ; ///< Base in radial direction
    int npoints ; ///< Number of coefficients

    // Constructors - Destructor
    // -------------------------
 protected:
    Diff(int base_r, int nr) ;	///< Standard constructor
    Diff(const Diff& ) ;		///< Copy constructor

    virtual ~Diff() ;			///< Destructor
 

    // Mutators / assignment
    // ---------------------
 protected:
    /// Assignment to another Diff
    void operator=(const Diff&) ;	
    
    // Accessors
    // ---------
 public:
    /// Returns the base on which the operator is defined
    int get_base() const {return base ;} ;

    /// Returns the number of coefficients (size of the matrix)
    int get_npoints() const { return npoints ;} ;

    /// Conversion to a matrix
    operator Matrice() const { return get_matrice() ; }
    
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const = 0 ;

    // Outputs
    // -------
 public:
    /// Display
    friend ostream& operator<<(ostream& , const Diff& ) ;	
    
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const = 0 ;

};


/**
 * Class for the elementary differential operator 
 * \f$ \frac{\partial}{\partial \xi} \f$ (see the base class \c Diff ).
 * \ingroup (ellip)
 * 
 */
class Diff_dsdx : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    /** Standard constructor, the base is that of the functions 
     * the operator is acting on (starting base).
     */
    Diff_dsdx(int base_r, int nr) ;	 

    Diff_dsdx(const Diff_dsdx& ) ;	 ///< Copy constructor

    virtual ~Diff_dsdx() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_dsdx
    void operator=(const Diff_dsdx&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};


/**
 * Class for the elementary differential operator 
 * \f$ \frac{\partial^2}{\partial \xi^2} \f$ (see the base class \c Diff ).
 * \ingroup (ellip)
 * 
 */
class Diff_dsdx2 : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_dsdx2(int base_r, int nr) ;	 ///< Standard constructor
    Diff_dsdx2(const Diff_dsdx2& ) ;	 ///< Copy constructor

    virtual ~Diff_dsdx2() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_dsdx2
    void operator=(const Diff_dsdx2&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \c Identity  (see the base class \c Diff ).
 * \ingroup (ellip)
 * 
 */
class Diff_id : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_id(int base_r, int nr) ;	 ///< Standard constructor
    Diff_id(const Diff_id& ) ;	 ///< Copy constructor

    virtual ~Diff_id() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_id
    void operator=(const Diff_id&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator multiplication by
 * \f$ \xi \f$ (see the base class \c Diff ). In the compactified
 * external domain, it corresponds to a multiplication by \f$\xi -1\f$.
 * It is not defined in the nucleus.
 * \ingroup (ellip)
 * 
 */
class Diff_mx : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_mx(int base_r, int nr) ;	 ///< Standard constructor
    Diff_mx(const Diff_mx& ) ;	 ///< Copy constructor

    virtual ~Diff_mx() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_mx
    void operator=(const Diff_mx&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator multiplication by
 * \f$ \xi^2 \f$ (see the base class \c Diff ). In the compactified
 * external domain, it corresponds to a multiplication by \f$(\xi -1)^2\f$.
 * \ingroup (ellip)
 * 
 */
class Diff_mx2 : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_mx2(int base_r, int nr) ;	 ///< Standard constructor
    Diff_mx2(const Diff_mx2& ) ;	 ///< Copy constructor

    virtual ~Diff_mx2() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_mx2
    void operator=(const Diff_mx2&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator division by
 * \f$ \xi - 1\f$ (see the base class \c Diff ). It is only defined 
 * in the compactified external domain.
 * It is not defined in the shells.
 * \ingroup (ellip)
 * 
 */
class Diff_sx : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_sx(int base_r, int nr) ;	 ///< Standard constructor
    Diff_sx(const Diff_sx& ) ;	 ///< Copy constructor

    virtual ~Diff_sx() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_sx
    void operator=(const Diff_sx&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator division by
 * \f$ \xi^2 \f$ (see the base class \c Diff ). In the compactified
 * external domain, it corresponds to a division by \f$(\xi -1)^2\f$. Not 
 * defined in the shells.
 * \ingroup (ellip)
 * 
 */
class Diff_sx2 : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_sx2(int base_r, int nr) ;	 ///< Standard constructor
    Diff_sx2(const Diff_sx2& ) ;	 ///< Copy constructor

    virtual ~Diff_sx2() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_sx2
    void operator=(const Diff_sx2&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \xi \frac{\partial}{\partial \xi} \f$ (see the base class \c Diff ).
 * In the compactified external domain the operator reads
 * \f$ (\xi -1)  \frac{\partial}{\partial \xi} \f$.
 * \ingroup (ellip)
 * 
 */
class Diff_xdsdx : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_xdsdx(int base_r, int nr) ;	 ///< Standard constructor
    Diff_xdsdx(const Diff_xdsdx& ) ;	 ///< Copy constructor

    virtual ~Diff_xdsdx() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_xdsdx
    void operator=(const Diff_xdsdx&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \frac{1}{\xi} \frac{\partial}{\partial \xi} \f$ (see the base 
 * class \c Diff ). In the compactified external domain, it reads
 * \f$ \frac{1}{\xi-1} \frac{\partial}{\partial \xi} \f$.
 * Not defined in the shells.
 * \ingroup (ellip)
 * 
 */
class Diff_sxdsdx : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_sxdsdx(int base_r, int nr) ;	 ///< Standard constructor
    Diff_sxdsdx(const Diff_sxdsdx& ) ;	 ///< Copy constructor

    virtual ~Diff_sxdsdx() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_sxdsdx
    void operator=(const Diff_sxdsdx&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \xi^2 \frac{\partial^2}{\partial \xi^2} \f$ (see the base 
 * class \c Diff ). In the external compactified domain, it reads
 * \f$ (\xi - 1)^2 \frac{\partial^2}{\partial \xi^2} \f$.
 * \ingroup (ellip)
 * 
 */
class Diff_x2dsdx2 : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_x2dsdx2(int base_r, int nr) ;	 ///< Standard constructor
    Diff_x2dsdx2(const Diff_x2dsdx2& ) ;	 ///< Copy constructor

    virtual ~Diff_x2dsdx2() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_x2dsdx2
    void operator=(const Diff_x2dsdx2&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \xi \frac{\partial^2}{\partial \xi^2} \f$ (see the base 
 * class \c Diff ). In the external compactified domain, it reads
 * \f$ (\xi - 1) \frac{\partial^2}{\partial \xi^2} \f$.
 * Not defined in the nucleus.
 * \ingroup (ellip)
 * 
 */
class Diff_xdsdx2 : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_xdsdx2(int base_r, int nr) ;	 ///< Standard constructor
    Diff_xdsdx2(const Diff_xdsdx2& ) ;	 ///< Copy constructor

    virtual ~Diff_xdsdx2() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_xdsdx2
    void operator=(const Diff_xdsdx2&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \xi^2 \frac{\partial}{\partial \xi} \f$ (see the base class \c Diff ).
 * This operator is not defined in the nucleus. In the compactified external
 * domain it reads \f$ (\xi-1)^2 \frac{\partial}{\partial \xi} \f$.
 * \ingroup (ellip)
 * 
 */
class Diff_x2dsdx : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_x2dsdx(int base_r, int nr) ;	 ///< Standard constructor
    Diff_x2dsdx(const Diff_x2dsdx& ) ;	 ///< Copy constructor

    virtual ~Diff_x2dsdx() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_x2dsdx
    void operator=(const Diff_x2dsdx&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \xi^3 \frac{\partial}{\partial \xi} \f$ (see the base class \c Diff ).
 * In the compactified external domain the operator reads
 * \f$ (\xi -1)^3  \frac{\partial}{\partial \xi} \f$.
 * \ingroup (ellip)
 * 
 */
class Diff_x3dsdx : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_x3dsdx(int base_r, int nr) ;	 ///< Standard constructor
    Diff_x3dsdx(const Diff_x3dsdx& ) ;	 ///< Copy constructor

    virtual ~Diff_x3dsdx() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_x3dsdx
    void operator=(const Diff_x3dsdx&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \xi^3 \frac{\partial^2}{\partial \xi^2} \f$ (see the base 
 * class \c Diff ). In the external compactified domain, it reads
 * \f$ (\xi - 1)^3 \frac{\partial^2}{\partial \xi^2} \f$.
 * Not defined in the nucleus.
 * \ingroup (ellip)
 * 
 */
class Diff_x3dsdx2 : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_x3dsdx2(int base_r, int nr) ;	 ///< Standard constructor
    Diff_x3dsdx2(const Diff_x3dsdx2& ) ;	 ///< Copy constructor

    virtual ~Diff_x3dsdx2() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_x3dsdx2
    void operator=(const Diff_x3dsdx2&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};

/**
 * Class for the elementary differential operator 
 * \f$ \xi^4 \frac{\partial^2}{\partial \xi^2} \f$ (see the base 
 * class \c Diff ). In the external compactified domain, it reads
 * \f$ (\xi - 1)^4 \frac{\partial^2}{\partial \xi^2} \f$.
 * \ingroup (ellip)
 * 
 */
class Diff_x4dsdx2 : public Diff {

    // Constructors - Destructor
    // -------------------------
 public:
    Diff_x4dsdx2(int base_r, int nr) ;	 ///< Standard constructor
    Diff_x4dsdx2(const Diff_x4dsdx2& ) ;	 ///< Copy constructor

    virtual ~Diff_x4dsdx2() ;			///< Destructor

 private:
    void initialize() ;  ///< Initializes arrays
 
    // Mutators / assignment
    // ---------------------
 public:
    /// Assignment to another Diff_x4dsdx2
    void operator=(const Diff_x4dsdx2&) ;	
	
    // Computational routines
    //-----------------------
 public:
    /// Returns the matrix associated with the operator
    virtual const Matrice& get_matrice() const ;

    // Outputs
    // -------
 protected:
    ///Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream&) const ;
};


}
#endif
