/*
 *  Definition of Lorene class Dim_tbl
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


#ifndef __DIM_TBL_H_ 
#define __DIM_TBL_H_ 

/**
 * Storage of array dimensions.
 * This class is designed for internal purposes related to the class
 * {\tt Tbl}, namely the storage of the {\tt Tbl} dimensions. 
 *
 * @version #$Id$#
 * @see	tbl.h
 */


/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.7  1999/11/23  12:16:24  eric
 * Modif commentaires (dimension 0 autorisee).
 *
 * Revision 2.6  1999/10/01  10:17:10  eric
 * Amelioration des commentaires.
 *
 * Revision 2.5  1999/09/30  12:49:09  eric
 * Constructeur a 1 parametre rendu explicit.
 * Amelioration des commentaires
 *
 * Revision 2.4  1999/09/24  14:22:36  eric
 * Declaration de methodes const
 * Amelioration commentaires.
 *
 * Revision 2.3  1999/09/22  11:38:59  eric
 * *** empty log message ***
 *
 * Revision 2.2  1999/09/22  11:24:39  eric
 * Amelioration commentaires
 *
 * Revision 2.1  1999/09/16  16:23:52  eric
 * Doc++
 *
 * Revision 2.0  1999/02/15  10:41:51  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

#include <stdio.h>

class Dim_tbl {
    public:
	int ndim ;	/// Number of dimensions of the {\tt Tbl}: can be 1, 2 or 3.
	int* dim ;	/// Array of dimensions (size: {\tt ndim}). 

	/**
	 * Total size of the array {\tt Tbl::t}. 
	 * \begin{itemize}
	 * \item[] {\tt taille = dim[0]} if {\tt ndim = 1}
	 * \item[] {\tt taille = dim[0]*dim[1]} if {\tt ndim = 2}
	 * \item[] {\tt taille = dim[0]*dim[1]*dim[2]} if {\tt ndim = 3}
	 * \end{itemize}
	 */
	int taille ;	
	
    public:
	// Constructeurs
	/**
	 * 1D constructor 
	 * @param size0  [input] Number of elements of the array {\tt Tbl::t}.
	 *		  Will be assigned to {\tt dim[0]}.  
	 *		  The size 0 is allowed for the 1D constructor but
	 *		  not for the 2D or 3D ones. 
	 */
	explicit Dim_tbl(int size0) ; 
	
	/**
	 * 2D constructor
	 * @param size1  [input] Defines the range [0, size1-1] of the outermost
	 *		  index in the storage of
	 *		  the array {\tt Tbl::t}. 
	 *		  Will be assigned to {\tt dim[1]}.
	 * @param size0  [input] Defines the range [0, size0-1]  of the 
	 *		  innermost index
	 *		  in the storage of
	 *		  the array {\tt Tbl::t}. 
	 *		  Will be assigned to {\tt dim[0]}.
	 */ 
	Dim_tbl(int size1, int size0) ;	    

	/**
	 * 3D constructor
	 * @param size2  [input] Defines the range [0, size2-1]  of the 
	 *		  outermost index in the storage of
	 *		  the array {\tt Tbl::t}. 
	 *		  Will be assigned to {\tt dim[2]}.
	 * @param size1  [input] Defines the range [0, size1-1] of the 
	 *		  intermediate index in the storage of
	 *		  the array {\tt Tbl::t}. 
	 *		  Will be assigned to {\tt dim[1]}.
	 * @param size0  [input] Defines the range [0, size0-1] of the 
	 *		  innermost index in the storage of
	 *		  the array {\tt Tbl::t}. 
	 *		  Will be assigned to {\tt dim[0]}.
	 */ 
	Dim_tbl(int size2, int size1, int size0) ; 

	Dim_tbl(const Dim_tbl & ) ; /// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )} )	
	explicit Dim_tbl(FILE* ) ;   	    

	~Dim_tbl() ;	    	    /// Destructor

	void operator=(const Dim_tbl &) ;   	/// Assignment
    
	void sauve(FILE* ) const ;			/// Save in a file

	bool operator==(const Dim_tbl &) const ;  /// Comparison operator

    	friend ostream& operator<<(ostream& , const Dim_tbl &) ;  /// Display
	
};

#endif
