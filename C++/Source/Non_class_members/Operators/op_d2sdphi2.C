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


 

/*
 * Ensemble des routines de base de derivation seconde par rapport a Phi
 * (Utilisation interne)
 * 
 *	void _d2sdphi2_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:07  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:24  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2006/03/10 12:45:38  j_novak
 * Use of C++-style cast.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/10/04  10:14:52  eric
 * Ajout d' abort() dans le cas non prevu.
 * ./
 *
 * Revision 2.3  2000/02/24  16:40:19  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.2  1999/11/15  16:37:33  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.1  1999/03/01  15:05:37  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/23  11:40:34  hyc
 * *** empty log message ***
 *
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tbl.h"

// Routine pour les cas non prevus
//--------------------------------
namespace Lorene {
void _d2sdphi2_pas_prevu(Tbl* , int & b) {
    cout << "Unknown phi basis in Mtbl_cf::d2sdp2() !" << endl ;
    cout << " basis: " << hex << b << endl ;
    abort() ; 
}

// cas P_COSSIN
//-------------
void _d2sdphi2_p_cossin(Tbl* tb, int & )
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cx = 0 ;
    static int np_pre =0 ;

    // Test sur np pour initialisation eventuelle
    if (np > np_pre) {
	np_pre = np ;
	cx = reinterpret_cast<double*>(realloc(cx, (np+2) * sizeof(double))) ;
	for (int i=0 ; i<np+2 ; i++) {
	    cx[i] = - (i/2) * (i/2) ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    // k = 0
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cx[0] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cx[k] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

}
