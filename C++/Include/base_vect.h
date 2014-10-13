/*
 *  Definition of Lorene classes Base_vect
 *				 Base_vect_cart
 *				 Base_vect_spher
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


#ifndef __BASE_VECT_H_ 
#define __BASE_VECT_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/13 08:52:31  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:09:39  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2004/03/22 13:12:40  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.4  2002/10/16 14:36:28  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
 * Revision 1.2  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/02/28  15:42:15  eric
 * Ajout de la fonction Base_vect_cart::get_align().
 *
 * Revision 2.5  2000/02/09  14:45:40  eric
 * *** empty log message ***
 *
 * Revision 2.4  2000/02/09  13:22:50  eric
 * REFONTE COMPLETE DE LA CLASSE
 * L'identification n'est plus base sur un membre statique (numero
 *  d'instance) mais sur les caracteres physiques (rot_phi, etc...)
 * Ajout des constructeurs par copie et lecture de fichier.
 *
 * Revision 2.3  2000/01/12  16:27:13  eric
 * Les constructeurs a un seul argument sont declares explicit.
 *
 * Revision 2.2  2000/01/10  13:34:52  eric
 * Modif commentairez.
 *
 * Revision 2.1  2000/01/10  13:26:49  eric
 * Ajout de la fonction set_rot_phi.
 *
 * Revision 2.0  2000/01/10  12:43:15  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdio>
#include "headcpp.h"

// Lorene classes
namespace Lorene {
class Tenseur ; 

		    //-----------------------------------//
		    //	    class Base_vect (base class) //
		    //-----------------------------------//

/**
 * Vectorial bases (triads) with respect to which the tensorial components are
 * defined. \ingroup (tensor)
 *
 */
class Base_vect {
    
    // Data : 
    // -----
    protected:
	char name[100] ;		///< Name of the basis


    // Constructors - Destructor
    // -------------------------
	
    protected:
	Base_vect() ;			    ///< Standard constructor

	/// Standard constructor with name
	explicit Base_vect(const char* name_i) ; 

	Base_vect(const Base_vect& ) ;	///< Copy constructor 
	 
    protected:
	/** Constructor from a file.
	 *  This constructor is protected because any \c Base_vect 
	 *  construction from a file must be done via the function 
	 *  \c Base_vect::bvect_from_file . 
	 */
	explicit Base_vect(FILE* ) ; 

    public:
	virtual ~Base_vect() ;			///< Destructor

    // Mutator / Assignment
    // --------------------

    private:
	 /// Assignement operator (not implemented).
	void operator=(const Base_vect& ) ;	  


    public: 
	/// Sets the basis name
	void set_name(const char* name_i) ; 


    // Extraction of information
    // -------------------------
    public:
	const char* get_name() const ;	///< Returns the basis name
	
	/** Returns a number to identify the sub-classe of \c Base_vect the
	 *  object belongs to. 
	 */
	virtual int identify() const = 0 ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	///< Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Base_vect& ) ;	

    protected: 
	virtual ostream& operator>>(ostream &) const = 0 ;    ///< Operator >>

    // Miscellaneous
    // -------------
    
    public:
	/** Construction of a vectorial basis from a file 
	 * (see \c sauve(FILE* ) ).
	 */
	static Base_vect* bvect_from_file(FILE* ) ; 

	/// Comparison operator (egality)
	virtual bool operator==(const Base_vect& ) const = 0 ;  

	/// Comparison operator (difference)
	bool operator!=(const Base_vect& ) const ;  
	
	/// Change the basis in which the components of a tensor are expressed
	virtual void change_basis(Tenseur& ) const = 0 ; 

};

ostream& operator<<(ostream& , const Base_vect& ) ;	


		    //-----------------------------------//
		    //	    class Base_vect_cart	 //
		    //-----------------------------------//


/**
 * Cartesian vectorial bases (triads). \ingroup (tensor)
 *
 */
class Base_vect_cart : public Base_vect {
    
    // Data : 
    // -----
    private:
	/// Angle between the \e x --axis and the absolute frame \e X --axis
	double rot_phi ;	

	/**
	 * Indicator of alignment with respect to the absolute frame: \\
	 *   \c align = 1  : basis aligned with the absolute frame 
	 *			(\f${\tt rot\_phi = 0}\f$) \\
	 *   \c align = -1  : basis anti-aligned with the absolute frame 
	 *			(\f${\tt rot\_phi} = \pi\f$) \\
	 *   \c align = 0  : general case 
	 */
	int align ; 

    // Constructors - Destructor
    // -------------------------
	
    public:
	explicit Base_vect_cart(double rot_phi_i) ;   ///< Standard constructor

	/// Standard constructor with name
	Base_vect_cart(double rot_phi_i, const char* name_i) ; 

 	Base_vect_cart(const Base_vect_cart& ) ;    ///< Copy constructor
	 
    protected:
	/** Constructor from a file.
	 *  This constructor is protected because any \c Base_vect_cart  
	 *  construction from a file must be done via the function 
	 *  \c Base_vect::bvect_from_file . 
	 */
	explicit Base_vect_cart(FILE* ) ; 

	/// The construction function from a file
	friend Base_vect* Base_vect::bvect_from_file(FILE* ) ; 

    public:
	virtual ~Base_vect_cart() ;			///< Destructor

    // Mutators
    // --------

    public:
	/// Assignment to another \c Base_vect_cart 
	void operator=(const Base_vect_cart& ) ;

	/** Sets a new value to the angle \c rot_phi 
	 *  between the \e x --axis and the absolute frame \e X --axis
	 */
	void set_rot_phi(double rot_phi_i) ; 

    private: 
	/// Computes \c align  from the value of \c rot_phi
	void set_align() ; 
	
    // Miscellaneous
    // -------------

    public:    	
	/// Comparison operator (egality)
	virtual bool operator==(const Base_vect& ) const ; 

	/// Change the basis in which the components of a tensor are expressed
	virtual void change_basis(Tenseur& ) const ; 

	/** Returns a number to identify the sub-classe of \c Base_vect the
	 *  object belongs to. 
	 */
	virtual int identify() const ; 

	/** Returns the
	 *  indicator of alignment with respect to the absolute frame.
	 * 
	 *  @return 
	 *   1 : basis aligned with the absolute frame 
	 *			(\f${\tt rot\_phi = 0}\f$) \\
	 *   -1 : basis anti-aligned with the absolute frame 
	 *			(\f${\tt rot\_phi} = \pi\f$) \\
	 *   0 : general case 
	 */
	int get_align() const {return align;} ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


};

		    //-----------------------------------//
		    //	    class Base_vect_spher	 //
		    //-----------------------------------//


/**
 * Spherical orthonormal vectorial bases (triads).\ingroup (tensor)
 *
 */
class Base_vect_spher : public Base_vect {
    
    // Data : 
    // -----
    private:
	double ori_x ;		///< Absolute coordinate \e X of the origin
	double ori_y ;		///< Absolute coordinate \e Y  of the origin
	double ori_z ;		///< Absolute coordinate \e Z of the origin
    
	/// Angle between the \e x --axis and the absolute frame \e X --axis
	double rot_phi ;	

	

    // Constructors - Destructor
    // -------------------------
	
    public:
	Base_vect_spher(double xa0, double ya0, double za0,  
			double rot_phi_i) ;   ///< Standard constructor

	/// Standard constructor with name
	Base_vect_spher(double xa0, double ya0, double za0, 
			double rot_phi_i, const char* name_i) ; 

	Base_vect_spher(const Base_vect_spher& ) ;	///< Copy constructor
	 
    protected:
	/** Constructor from a file.
	 *  This constructor is protected because any \c Base_vect_spher  
	 *  construction from a file must be done via the function 
	 *  \c Base_vect::bvect_from_file . 
	 */
	explicit Base_vect_spher(FILE* ) ; 

	/// The construction function from a file
	friend Base_vect* Base_vect::bvect_from_file(FILE* ) ; 

    public:
	virtual ~Base_vect_spher() ;			///< Destructor

    // Mutators
    // --------
    public:
	/// Assignment to another \c Base_vect_spher 
	void operator=(const Base_vect_spher& ) ;

    public:
	/// Sets a new origin
	void set_ori(double xa0, double ya0, double za0) ;  

	/** Sets a new value to the angle \c rot_phi 
	 *  between the \e x --axis and the absolute frame \e X --axis
	 */
	void set_rot_phi(double rot_phi_i) ; 

    // Miscellaneous
    // -------------
    	
    public:    	
	/// Comparison operator (egality)
	virtual bool operator==(const Base_vect& ) const ; 

	/// Change the basis in which the components of a tensor are expressed
	virtual void change_basis(Tenseur& ) const ; 

	/** Returns a number to identify the sub-classe of \c Base_vect the
	 *  object belongs to. 
	 */
	virtual int identify() const ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


};

}
#endif
