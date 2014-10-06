/*
 * Methods of class Dim_tbl
 *
 *  (see file dim_tbl.h for documentation)
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


char dim_tbl[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/06 15:13:18  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2006/09/26 07:21:07  p_grandclement
 * Minor change in the indices
 *
 * Revision 1.4  2006/09/25 10:01:50  p_grandclement
 * Addition of N-dimensional Tbl
 *
 * Revision 1.3  2002/10/16 14:37:13  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/12/04 21:27:54  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1  2001/11/23 09:37:51  e_gourgoulhon
 * dim_tbl.C now in directory Tbl
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  1999/11/23  12:16:49  eric
 * Dimension 0 autorisee dans le constructeur 1D.
 *
 * Revision 2.4  1999/09/24  14:24:09  eric
 * Declaration de methodes const.
 *
 * Revision 2.3  1999/09/22  11:24:53  eric
 * Correction erreur ecriture/lecture fichier de taille (double->int).
 * Initialisation de ndim.
 *
 * Revision 2.2  1999/09/16  16:24:08  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/03/01  14:56:17  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// Headers C
#include <cassert>

// Headers Lorene
#include "dim_tbl.h"
#include "utilitaires.h"

			//---------------//
			// Constructeurs //
			//---------------//

// 1D constructor
Dim_tbl::Dim_tbl(int i) : ndim(1) {
    assert(i >= 0) ;			// The dimension 0 is allowed
    dim = new int[ndim] ;
    dim[0] = i ;
    taille = i ;
}
// 2D constructor
Dim_tbl::Dim_tbl(int j, int i) : ndim(2) {
    assert(j > 0) ;
    assert(i > 0) ;
    dim = new int[ndim] ;
    dim[0] = i ; dim[1] = j ;
    taille = i * j ;
}
// 3D constructor
Dim_tbl::Dim_tbl(int k, int j, int i) : ndim(3) {
    assert(k > 0) ;
    assert(j > 0) ;
    assert(i > 0) ;
    dim = new int[ndim] ;
    dim[0] = i ; dim[1] = j ; dim[2] = k ;
    taille = i * j * k ;
}
	
// N-dimensional constructor
Dim_tbl::Dim_tbl(int n, int* sizes) : ndim(n) {
    for (int i=0 ; i<ndim ; i++)
        assert(sizes[i] > 0) ;
    dim = new int[ndim] ;
    taille = 1 ;
    for (int i=0 ; i<ndim ; i++) {
    	dim[i] = sizes[ndim-i-1] ;
	taille *= sizes[i] ;
    }
}

// Copy
Dim_tbl::Dim_tbl(const Dim_tbl & titi) : ndim(titi.ndim) {
    dim = new int[ndim] ;
    for (int i=0 ; i<ndim ; i++) {
    	dim[i] = titi.dim[i] ;
    }
    taille = titi.taille ;
}
	
// From a file
Dim_tbl::Dim_tbl(FILE* fd) {
    fread_be(&ndim, sizeof(int), 1, fd) ;		// ndim
    dim = new int[ndim] ;
    fread_be(dim, sizeof(int), ndim, fd) ;		// dim[]
    taille = dim[0] ;
    for (int i=1; i<ndim; i++) {
	taille *= dim[i] ; 
    }
}

			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Dim_tbl::~Dim_tbl() {
    delete [] dim ;
}

			//-------------//
			// Affectation //
			//-------------//

// From Dim_tbl
void Dim_tbl::operator=(const Dim_tbl & titi) {
    ndim = titi.ndim ;
    delete [] dim ;
    dim = new int[ndim] ;
    for (int i=0 ; i<ndim ; i++) {
    	dim[i] = titi.dim[i] ;
    }
    taille = titi.taille ;
}
    
			//------------//
			// Sauvegarde //
			//------------//

// Save in a file
void Dim_tbl::sauve(FILE* fd) const {
    fwrite_be(&ndim, sizeof(int), 1, fd) ;		    // ndim
    fwrite_be(dim, sizeof(int), ndim, fd) ;	    // dim[]
}
    
			//------------//
			// Impression //
			//------------//

// Operateurs <<
ostream& operator<<(ostream& o, const Dim_tbl & titi) {
    o << titi.ndim << " dimension(s):" ;
    for (int i=0 ; i<titi.ndim ; i++) {
    	o << " " << titi.dim[i] ;
    }
    return o ;
}


		    //---------------------//
	    	    // Operateurs logiques //
		    //---------------------//

bool Dim_tbl::operator==(const Dim_tbl & ti) const {
    
    // Peut-etre faux ?
    if (ndim != ti.ndim) return false ;
    for (int i=0 ; i<ndim ; i++) {
    	if (dim[i] != ti.dim[i]) return false ;
    }
    
    // Non ! juste
    return true ;
}
