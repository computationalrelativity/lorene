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
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
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

// Standard C++
class ostream ; 

// Headers C
#include <stdio.h>

// Lorene classes
class Tenseur ; 

		    //-----------------------------------//
		    //	    class Base_vect (base class) //
		    //-----------------------------------//

/**
 * Vectorial bases (triads) with respect to which the tensorial components are
 * defined. 
 *
 * @version #$Id$#
 */
class Base_vect {
    
    // Data : 
    // -----
    protected:
	char name[100] ;		/// Name of the basis


    // Constructors - Destructor
    // -------------------------
	
    protected:
	Base_vect() ;			    /// Standard constructor

	/// Standard constructor with name
	explicit Base_vect(const char* name_i) ; 

	Base_vect(const Base_vect& ) ;	/// Copy constructor 
	 
    protected:
	/** Constructor from a file.
	 *  This constructor is protected because any {\tt Base\_vect} 
	 *  construction from a file must be done via the function 
	 *  {\tt Base\_vect::bvect\_from\_file}. 
	 */
	explicit Base_vect(FILE* ) ; 

    public:
	virtual ~Base_vect() ;			/// Destructor

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
	const char* get_name() const ;	/// Returns the basis name
	
	/** Returns a number to identify the sub-classe of {\tt Base\_vect} the
	 *  object belongs to. 
	 */
	virtual int identify() const = 0 ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	/// Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Base_vect& ) ;	

    protected: 
	virtual ostream& operator>>(ostream &) const = 0 ;    /// Operator >>

    // Miscellaneous
    // -------------
    
    public:
	/** Construction of a vectorial basis from a file 
	 * (see {\tt sauve(FILE* )}).
	 */
	static Base_vect* bvect_from_file(FILE* ) ; 

	/// Comparison operator (egality)
	virtual bool operator==(const Base_vect& ) const = 0 ;  

	/// Comparison operator (difference)
	bool operator!=(const Base_vect& ) const ;  
	
	/// Change the basis in which the components of a tensor are expressed
	virtual void change_basis(Tenseur& ) const = 0 ; 

};


		    //-----------------------------------//
		    //	    class Base_vect_cart	 //
		    //-----------------------------------//


/**
 * Cartesian vectorial bases (triads). 
 *
 * @version #$Id$#
 */
class Base_vect_cart : public Base_vect {
    
    // Data : 
    // -----
    private:
	/// Angle between the $x$--axis and the absolute frame $X$--axis
	double rot_phi ;	

	/**
	 * Indicator of alignment with respect to the absolute frame: \\
	 *   {\tt align = 1} : basis aligned with the absolute frame 
	 *			(${\tt rot\_phi = 0}$) \\
	 *   {\tt align = -1} : basis anti-aligned with the absolute frame 
	 *			(${\tt rot\_phi} = \pi$) \\
	 *   {\tt align = 0} : general case 
	 */
	int align ; 

    // Constructors - Destructor
    // -------------------------
	
    public:
	explicit Base_vect_cart(double rot_phi_i) ;   /// Standard constructor

	/// Standard constructor with name
	Base_vect_cart(double rot_phi_i, const char* name_i) ; 

 	Base_vect_cart(const Base_vect_cart& ) ;    /// Copy constructor
	 
    protected:
	/** Constructor from a file.
	 *  This constructor is protected because any {\tt Base\_vect\_cart} 
	 *  construction from a file must be done via the function 
	 *  {\tt Base\_vect::bvect\_from\_file}. 
	 */
	explicit Base_vect_cart(FILE* ) ; 

	/// The construction function from a file
	friend Base_vect* Base_vect::bvect_from_file(FILE* ) ; 

    public:
	virtual ~Base_vect_cart() ;			/// Destructor

    // Mutators
    // --------

    public:
	/// Assignment to another {\tt Base\_vect\_cart}
	void operator=(const Base_vect_cart& ) ;

	/** Sets a new value to the angle {\tt rot\_phi} 
	 *  between the $x$--axis and the absolute frame $X$--axis
	 */
	void set_rot_phi(double rot_phi_i) ; 

    private: 
	/// Computes {\tt align} from the value of {\tt rot\_phi}
	void set_align() ; 
	
    // Miscellaneous
    // -------------

    public:    	
	/// Comparison operator (egality)
	virtual bool operator==(const Base_vect& ) const ; 

	/// Change the basis in which the components of a tensor are expressed
	virtual void change_basis(Tenseur& ) const ; 

	/** Returns a number to identify the sub-classe of {\tt Base\_vect} the
	 *  object belongs to. 
	 */
	virtual int identify() const ; 

	/** Returns the
	 *  indicator of alignment with respect to the absolute frame.
	 * 
	 *  @return 
	 *   1 : basis aligned with the absolute frame 
	 *			(${\tt rot\_phi = 0}$) \\
	 *   -1 : basis anti-aligned with the absolute frame 
	 *			(${\tt rot\_phi} = \pi$) \\
	 *   0 : general case 
	 */
	int get_align() const {return align;} ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	/// Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};

		    //-----------------------------------//
		    //	    class Base_vect_spher	 //
		    //-----------------------------------//


/**
 * Spherical orthonormal vectorial bases (triads).
 *
 * @version #$Id$#
 */
class Base_vect_spher : public Base_vect {
    
    // Data : 
    // -----
    private:
	double ori_x ;		/// Absolute coordinate $X$ of the origin
	double ori_y ;		/// Absolute coordinate $Y$ of the origin
	double ori_z ;		/// Absolute coordinate $Z$ of the origin
    
	/// Angle between the $x$--axis and the absolute frame $X$--axis
	double rot_phi ;	

	

    // Constructors - Destructor
    // -------------------------
	
    public:
	Base_vect_spher(double xa0, double ya0, double za0,  
			double rot_phi_i) ;   /// Standard constructor

	/// Standard constructor with name
	Base_vect_spher(double xa0, double ya0, double za0, 
			double rot_phi_i, const char* name_i) ; 

	Base_vect_spher(const Base_vect_spher& ) ;	/// Copy constructor
	 
    protected:
	/** Constructor from a file.
	 *  This constructor is protected because any {\tt Base\_vect\_spher} 
	 *  construction from a file must be done via the function 
	 *  {\tt Base\_vect::bvect\_from\_file}. 
	 */
	explicit Base_vect_spher(FILE* ) ; 

	/// The construction function from a file
	friend Base_vect* Base_vect::bvect_from_file(FILE* ) ; 

    public:
	virtual ~Base_vect_spher() ;			/// Destructor

    // Mutators
    // --------
    public:
	/// Assignment to another {\tt Base\_vect\_spher}
	void operator=(const Base_vect_spher& ) ;

    public:
	/// Sets a new origin
	void set_ori(double xa0, double ya0, double za0) ;  

	/** Sets a new value to the angle {\tt rot\_phi} 
	 *  between the $x$--axis and the absolute frame $X$--axis
	 */
	void set_rot_phi(double rot_phi_i) ; 

    // Miscellaneous
    // -------------
    	
    public:    	
	/// Comparison operator (egality)
	virtual bool operator==(const Base_vect& ) const ; 

	/// Change the basis in which the components of a tensor are expressed
	virtual void change_basis(Tenseur& ) const ; 

	/** Returns a number to identify the sub-classe of {\tt Base\_vect} the
	 *  object belongs to. 
	 */
	virtual int identify() const ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	/// Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};

#endif
