/*
 *  Definition of Lorene class XXX
 *
 */

/*
 *   Copyright (c) year  your_name
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

#ifndef __XXX_H_ 
#define __XXX_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 * template files
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <>

// C headers
#include <>

// Lorene headers
#include ""

/**
 * Extended description of the class for Doc++ documentation
 * 
 * @version #$Id$#
 */
class XXX {

    // Data : 
    // -----
    protected:

    // Derived data : 
    // ------------
    protected:
	mutable ?? p_?? ;   /// Comment for Doc++

    // Constructors - Destructor
    // -------------------------
    public:
	XXX(??) ;			/// Standard constructor
	XXX(const XXX& ) ;		/// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	XXX(FILE* ) ;    		

	virtual ~XXX() ;			/// Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Logical destructor
	virtual void del_t() ;	
	    
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	virtual void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another XXX
	void operator=(const XXX&) ;	
	
    // Accessors
    // ---------
    public:

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const XXX& ) ;	



};

#endif
