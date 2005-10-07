/*
 * Methods of class Mg3d
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


char mg3d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.11  2005/10/07 08:47:21  j_novak
 * Addition of the pointer g_non_axi on a grid, with at least 5 points in the
 * theta direction and 4 in the phi one (for tensor rotations).
 *
 * Revision 1.10  2004/07/06 13:36:29  j_novak
 * Added methods for desaliased product (operator |) only in r direction.
 *
 * Revision 1.9  2003/10/27 16:21:54  e_gourgoulhon
 * Treated the special case nz=1 in the simplified constructor.
 *
 * Revision 1.8  2003/06/20 14:50:15  f_limousin
 * Add the operator==
 *
 * Revision 1.7  2003/06/18 08:45:27  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
 * Revision 1.6  2002/10/16 14:36:42  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.5  2002/05/07 07:36:03  e_gourgoulhon
 * Compatibilty with xlC compiler on IBM SP2:
 *    suppressed the parentheses around argument of instruction new:
 * 	e.g.   t = new (Tbl *[nzone])  -->   t = new Tbl*[nzone]
 *
 * Revision 1.4  2001/12/12 09:23:46  e_gourgoulhon
 * Parameter compact added to the simplified constructor of class Mg3d
 *
 * Revision 1.3  2001/12/11 06:48:30  e_gourgoulhon
 * Addition of the simplified constructor
 *
 * Revision 1.2  2001/12/04 21:27:54  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.10  2001/05/26  14:50:46  eric
 * *** empty log message ***
 *
 * Revision 2.9  2001/05/26  13:25:59  eric
 * Ajout du membre g_twice (grille double pour le desaliasing)
 * Modif de la declaration de g_angu (pointeur mutable)
 *   g_twice et g_angu ne sont calcules que si necessaire (cad si
 *   on appelle la fonction get_twice() ou get_angu()).
 *
 * Revision 2.8  2000/03/22  13:38:51  eric
 * Remplacement des iendl par endl dans <<
 *
 * Revision 2.7  1999/10/12  15:04:29  eric
 * *** empty log message ***
 *
 * Revision 2.6  1999/10/12  15:03:30  eric
 * *** empty log message ***
 *
 * Revision 2.5  1999/09/30  14:58:16  eric
 * Operator!= declare const
 *
 * Revision 2.4  1999/09/30  14:12:04  eric
 * sauve declaree const.
 *
 * Revision 2.3  1999/09/30  12:52:52  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.2  1999/03/01  14:35:21  eric
 * Modif affichage (operator<<)
 *
 *
 * $Header$
 *
 */

// C Headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "grilles.h"
#include "type_parite.h"
#include "utilitaires.h"

		//--------------//
		// Multi-grille //
		//--------------//


//=============================================================================
//    General constructor
//=============================================================================


Mg3d::Mg3d(int nz,
    int nbr[], int typr[], int nbt[], int typt, int nbp[], int typp)
    : nzone(nz), type_t(typt), type_p(typp)
{

    // Type d'echantillonnage dans chaque zone
    type_r = new int[nz];
    for (int i=0 ; i<nz ; i++) {
	type_r[i] = typr[i];
    }

    // Nombre de points
    nr = new int[nz];
    nt = new int[nz];
    np = new int[nz];
    for (int i=0 ; i<nz ; i++) {
	nr[i] = nbr[i] ;
	nt[i] = nbt[i] ;
	np[i] = nbp[i] ;
    }

    // Les grilles
    // -----------
    g = new Grille3d*[nz] ;

    for (int i=0; i<nz; i++) {

    //... Echantillonnage selon le type demande:
    switch (type_p) {
//----------------------------------------------------------------------------
	case NONSYM :    // echantillonnage en phi sur [0, 2 pi[
//----------------------------------------------------------------------------
	    switch (type_t) {
		case SYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_feq(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_req(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_ieq(nbr[i], nbt[i], nbp[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = SYM et type_p = NONSYM

		case NONSYM :    // echantillonnage en theta sur [0,pi]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_f(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_r(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_i(nbr[i], nbt[i], nbp[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = NONSYM et type_p = NONSYM

		default :
		    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			 << " et type_t = " << type_t <<
			    " n'est pas prevu !" << endl ;
		    abort() ;
		
	    }  // fin des differents types d'echantillonnage en theta
	    break ;	// fin du cas type_p = NONSYM

//----------------------------------------------------------------------------
	case SYM :    // echantillonnage en phi sur [0,pi[
//----------------------------------------------------------------------------
	    switch (type_t) {
		case SYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_fs(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_rs(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_is(nbr[i], nbt[i], nbp[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = SYM et type_p = SYM

		case NONSYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_f2p(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_r2p(nbr[i], nbt[i], nbp[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_i2p(nbr[i], nbt[i], nbp[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = NONSYM et type_p = SYM

		default :
		    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			 << " et type_t = " << type_t <<
			    " n'est pas prevu !" << endl ;
		    abort() ;
		
	    }  // fin des differents types d'echantillonnage en theta
	    break ;	// fin du cas type_p = SYM

//----------------------------------------------------------------------------
	default :
//----------------------------------------------------------------------------
	    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
		 << " n'est pas prevu !" << endl ;
	    abort() ;
    }	// fin des differents types d'echantillonnage en phi
}   // fin de la boucle sur les zones

    // Pointers on derived grids initiated to 0x0:
    // -------------------------------------------

    set_deriv_0x0() ;


}

//=============================================================================
//    Simplified constructor
//=============================================================================

Mg3d::Mg3d(int nz, int nbr, int nbt, int nbp, int typt, int typp, 
	   bool compact)
    		: nzone(nz),
    		  type_t(typt),
    		  type_p(typp)   {

    // Type of r sampling in each domain: 
    type_r = new int[nz];
    type_r[0] = RARE ;
    for (int l=1; l<nz-1; l++) {
	type_r[l] = FIN ;
    }
    if (nz > 1) {
		type_r[nz-1] = compact ? UNSURR : FIN ;
	}

    // Same number of points in all domains:
    nr = new int[nz];
    nt = new int[nz];
    np = new int[nz];
    for (int l=0 ; l<nz ; l++) {
		nr[l] = nbr ;
		nt[l] = nbt ;
		np[l] = nbp ;
    }

    // Les grilles
    // -----------
    g = new Grille3d*[nz] ;

    for (int i=0; i<nz; i++) {

    //... Echantillonnage selon le type demande:
    switch (type_p) {
//----------------------------------------------------------------------------
	case NONSYM :    // echantillonnage en phi sur [0, 2 pi[
//----------------------------------------------------------------------------
	    switch (type_t) {
		case SYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_feq(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_req(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_ieq(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = SYM et type_p = NONSYM

		case NONSYM :    // echantillonnage en theta sur [0,pi]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_f(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_r(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_i(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = NONSYM et type_p = NONSYM

		default :
		    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			 << " et type_t = " << type_t <<
			    " n'est pas prevu !" << endl ;
		    abort() ;
		
	    }  // fin des differents types d'echantillonnage en theta
	    break ;	// fin du cas type_p = NONSYM

//----------------------------------------------------------------------------
	case SYM :    // echantillonnage en phi sur [0,pi[
//----------------------------------------------------------------------------
	    switch (type_t) {
		case SYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_fs(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_rs(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_is(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = SYM et type_p = SYM

		case NONSYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_f2p(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_r2p(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_i2p(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = NONSYM et type_p = SYM

		default :
		    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			 << " et type_t = " << type_t <<
			    " n'est pas prevu !" << endl ;
		    abort() ;
		
	    }  // fin des differents types d'echantillonnage en theta
	    break ;	// fin du cas type_p = SYM

//----------------------------------------------------------------------------
	default :
//----------------------------------------------------------------------------
	    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
		 << " n'est pas prevu !" << endl ;
	    abort() ;
    }	// fin des differents types d'echantillonnage en phi
}   // fin de la boucle sur les zones

    // Pointers on derived grids initiated to 0x0:
    // -------------------------------------------

    set_deriv_0x0() ;
	
}



//=============================================================================
//    Constructor from a file
//=============================================================================

/*
 * Construction a partir d'un fichier.
 * Cette facon de faire est abominable. Cependant je ne vois pas comment
 * faire autrement... j.a.
 */
Mg3d::Mg3d(FILE* fd)
{
    // Lecture sur le fichier
    fread_be(&nzone, sizeof(int), 1, fd) ;		// nzone
    nr = new int[nzone] ;
    fread_be(nr, sizeof(int), nzone, fd) ;		// nr
    nt = new int[nzone] ;
    fread_be(nt, sizeof(int), nzone, fd) ;		// nt
    np = new int[nzone] ;
    fread_be(np, sizeof(int), nzone, fd) ;		// np
    type_r = new int[nzone] ;
    fread_be(type_r, sizeof(int), nzone, fd) ;	// type_r
    fread_be(&type_t, sizeof(int), 1, fd) ;	// type_t
    fread_be(&type_p, sizeof(int), 1, fd) ;	// type_p

    // Les grilles
    // -----------

    g = new Grille3d*[nzone] ;

    for (int i=0; i<nzone; i++) {

    //... Echantillonnage selon le type demande:
    switch (type_p) {
//----------------------------------------------------------------------------
	case NONSYM :    // echantillonnage en phi sur [0, 2 pi[
//----------------------------------------------------------------------------
	    switch (type_t) {
		case SYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_feq(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_req(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_ieq(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = SYM et type_p = NONSYM

		case NONSYM :    // echantillonnage en theta sur [0,pi]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_f(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_r(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_i(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = NONSYM et type_p = NONSYM

		default :
		    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			 << " et type_t = " << type_t <<
			    " n'est pas prevu !" << endl ;
		    abort() ;
		
	    }  // fin des differents types d'echantillonnage en theta
	    break ;	// fin du cas type_p = NONSYM

//----------------------------------------------------------------------------
	case SYM :    // echantillonnage en phi sur [0,pi[
//----------------------------------------------------------------------------
	    switch (type_t) {
		case SYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_fs(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_rs(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_is(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = SYM et type_p = SYM

		case NONSYM :    // echantillonnage en theta sur [0,pi/2]
		    switch (type_r[i]) {
			case FIN :    // echantillonnage fin
			    g[i] = new
			     Grille3d_f2p(nr[i], nt[i], np[i]) ;
			    break ;
			case RARE :    // echantillonnage rarefie
			    g[i] = new
			     Grille3d_r2p(nr[i], nt[i], np[i]) ;
			    break ;
			case UNSURR :    // echantillonnage fin en 1/r
			    g[i] = new
			     Grille3d_i2p(nr[i], nt[i], np[i]) ;
			    break ;

			default :
			cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			     << " , type_t = " << type_t
			     << " et type_r = " << type_r[i] << endl ;
			cout << "n'est pas prevu !" << endl ;
			abort() ;
			
		    }	// fin des differents types d'echantillonnage en r
		    break ;  // fin du cas type_t = NONSYM et type_p = SYM

		default :
		    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
			 << " et type_t = " << type_t <<
			    " n'est pas prevu !" << endl ;
		    abort() ;
		
	    }  // fin des differents types d'echantillonnage en theta
	    break ;	// fin du cas type_p = SYM

//----------------------------------------------------------------------------
	default :
//----------------------------------------------------------------------------
	    cout << "Mg3d::Mg3d : Le cas type_p = " << type_p
		 << " n'est pas prevu !" << endl ;
	    abort() ;
    }	// fin des differents types d'echantillonnage en phi
}   // fin de la boucle sur les zones

    // Pointers on derived grids initiated to 0x0:
    // -------------------------------------------

    set_deriv_0x0() ;

}


// Destructeur
// -----------
Mg3d::~Mg3d() {

    del_deriv() ;   // Deletes the derived quantities

    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    for (int i=0 ; i<nzone ; i++) {
	delete g[i] ;
    }
    delete [] g ;

}

//==================================================================
//  Write in a file
//==================================================================

void Mg3d::sauve(FILE* fd) const {	
	    fwrite_be(&nzone, sizeof(int), 1, fd) ;	// nzone
	    fwrite_be(nr, sizeof(int), nzone, fd) ;	// nr
	    fwrite_be(nt, sizeof(int), nzone, fd) ;	// nt
	    fwrite_be(np, sizeof(int), nzone, fd) ;	// np
	    fwrite_be(type_r, sizeof(int), nzone, fd) ;	// type_r
	    fwrite_be(&type_t, sizeof(int), 1, fd) ;	// type_t
	    fwrite_be(&type_p, sizeof(int), 1, fd) ;	// type_p
}




		//--------------------------//
		// Surcharge des operateurs //
		//--------------------------//

// Operateur <<
ostream& operator<<(ostream& o, const Mg3d& g) {
    char* tr[3] ;
    tr[FIN] = "FIN" ; tr[RARE] = "RARE" ; tr[UNSURR] = "UNSURR" ;
    char* tang[2] ;
    tang[NONSYM] = "NONSYM" ; tang[SYM] = "SYM" ;
    o << "Number of domains: " << g.nzone << endl ;
    for (int i=0 ; i< g.nzone ; i++) {
	o << "  Domain #" << i << ": "
	<< "nr = " << g.nr[i] << ", " << tr[g.type_r[i]] << "; "
	<< "nt = " << g.nt[i] << ", " << tang[g.type_t] << "; "
	<< "np = " << g.np[i] << ", " << tang[g.type_p] << endl ;
    }
    o << endl ;
    return o ;
}

// Operateur !=
bool Mg3d::operator!=(const Mg3d & titi) const {

    if (nzone != titi.nzone) return true ;   // C'est vrai que c'est faux...

    for (int i=0 ; i<nzone ; i++) {
	if (nr[i] != titi.nr[i]) return true ;
	if (nt[i] != titi.nt[i]) return true ;
	if (np[i] != titi.np[i]) return true ;

	if (type_r[i] != titi.type_r[i]) return true ;
    }

    if (type_t != titi.type_t) return true ;
    if (type_p != titi.type_p) return true ;

    // C'est faux que c'est vrai...
    return false ;
}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Mg3d::del_deriv() const {

    if (g_angu != 0x0) delete g_angu ;
    if (g_radial != 0x0) delete g_radial ;
    if (g_twice != 0x0) delete g_twice ;
    if (g_plus_half != 0x0) delete g_plus_half ;
    if (g_non_axi != 0x0) delete g_non_axi ;

    set_deriv_0x0() ;

}

void Mg3d::set_deriv_0x0() const {

    g_angu = 0x0 ;
    g_radial = 0x0 ;
    g_twice = 0x0 ;
    g_plus_half = 0x0 ;
    g_non_axi = 0x0 ;
}


			    //--------------//
			    // Angular grid //
			    //--------------//

const Mg3d* Mg3d::get_angu() const {

    if (g_angu == 0x0) {	  // The construction is required

	int* nbr_angu = new int[nzone] ;
	for (int i=0 ; i<nzone ; i++) {
	    nbr_angu[i] = 1 ;
	}
	g_angu = new Mg3d(nzone, nbr_angu, type_r, nt, type_t, np, type_p) ;
	delete [] nbr_angu ;
    }

    return g_angu ;

}
	
			    //--------------//
			    //  Radial grid //
			    //--------------//

const Mg3d* Mg3d::get_radial() const {

    if (g_radial == 0x0) {	  // The construction is required

	int* nbr_radial = new int[nzone] ;
	for (int i=0 ; i<nzone ; i++) {
	    nbr_radial[i] = 1 ;
	}
	g_radial = new Mg3d(nzone, nr, type_r, nbr_radial, SYM, nbr_radial, 
			  SYM) ;
	delete [] nbr_radial ;
    }

    return g_radial ;

}
	
		  //--------------------------------------//
		  // Grid with twice the number of points //
		  //--------------------------------------//

const Mg3d* Mg3d::get_twice() const {

    if (g_twice == 0x0) {	  // The construction is required

	int* nbr = new int[nzone] ;
	int* nbt = new int[nzone] ;
	int* nbp = new int[nzone] ;

	for (int l=0; l<nzone; l++) {
	    if (nr[l] == 1) {
		nbr[l] = 1 ;
	    }
	    else {
		nbr[l] = 2*nr[l] - 1 ;
	    }

	    if (nt[l] == 1) {
		nbt[l] = 1 ;
	    }
	    else {
		nbt[l] = 2*nt[l] - 1 ;
	    }
	
	    if (np[l] == 1) {
		nbp[l] = 1 ;
	    }
	    else {
		nbp[l] = 2*np[l] ;
	    }
	}

	g_twice = new Mg3d(nzone, nbr, type_r, nbt, type_t, nbp, type_p) ;

	delete [] nbr ;
	delete [] nbt ;
	delete [] nbp ;

    }

    return g_twice ;

}
	

		  //--------------------------------------//
		  // Grid with 50% additional points in r //
		  //--------------------------------------//

const Mg3d* Mg3d::plus_half() const {

  if (g_plus_half == 0x0) {	  // The construction is required

    int* nbr = new int[nzone] ;
    
    for (int l=0; l<nzone; l++) {
      if (nr[l] == 1) 
	nbr[l] = 1 ;
      else 
	nbr[l] = (3*nr[l])/2  ;
    }
    
    g_plus_half = new Mg3d(nzone, nbr, type_r, nt, type_t, np, type_p) ;

    delete [] nbr ;


  }

  return g_plus_half ;

}
	
		  //----------------------------------------------//
		  // Grid for rotations (5/4 points in theta/phi) //
		  //----------------------------------------------//

const Mg3d* Mg3d::get_non_axi() const {

  if (g_non_axi == 0x0) {	  // The construction is required

    int* nbt = new int[nzone] ;
    int* nbp = new int[nzone] ;
    
    for (int l=0; l<nzone; l++) {
      if (nt[l] < 5)   
	  nbt[l] = 5 ;
      else
	  nbt[l] = nt[l] ;
      if (np[l] < 4)
	  nbp[l] = 4 ;
      else
	  nbp[l] = np[l] ;
    }
    
    g_non_axi = new Mg3d(nzone, nr, type_r, nbt, type_t, nbp, type_p) ;

    delete [] nbt ;
    delete [] nbp ;


  }

  return g_non_axi ;

}
	

bool Mg3d::operator==(const Mg3d& mgi) const {
  
  bool resu = true ;

  if (mgi.get_nzone() != nzone) {
    resu = false ;
  }
  else {
    for (int i=0; i<nzone; i++) {
      if (mgi.get_nr(i) != nr[i]) resu = false ;
      if (mgi.get_np(i) != np[i]) resu = false ;
      if (mgi.get_nt(i) != nt[i]) resu = false ;
      if (mgi.get_type_r(i) != type_r[i]) resu =false ;
    }
  }
  
  if (mgi.get_type_t() != type_t) resu = false ;
  if (mgi.get_type_p() != type_p) resu = false ;

  return resu ;

}
