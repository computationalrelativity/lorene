/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char op_sxm1_C[] = "$Header$" ;

/* 
 * Ensemble des routines de base pour l'operateur 1/(x-1) Id
 * (Utilisation interne)
 * 
 *  Prototype :
 *  ---------
 *	void _sxm1_XXXX(Tbl * tb, int& base)
 * 
 *  Entree/Sortie :
 *  -------------
 *	tb	pointeur sur le Tbl d'entree-sortie
 *
 *  Entree :
 *  ------
 *	base	base de travail
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:26  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  1999/11/15  16:39:47  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.1  1999/04/26  16:02:21  phil
 * changement commentaire (!)
 *
 * Revision 2.0  1999/04/26  15:54:52  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//Lorene
#include "tbl.h"

// Prototypage
#include "proto.h"

		//-----------------------------------
		// Routine qui ne fait rien	   --
		//-----------------------------------

namespace Lorene {
void _sxm1_identite(Tbl * , int& ) {
    return ; 
}

			//-------------------------
			// cas R_CHEB et R_CHEBU --
			//-------------------------

void _sxm1_cheb(Tbl* tb, int&)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    int ntnr = nt*nr ; 
    
    int k, j ; 
    for (k=0 ; k<np+1 ; k++) {
	if (k==1) continue ;	// On ne traite pas le coefficient de sin(0*phi) 
	for (j=0 ; j<nt ; j++) {
	
	    double* cf = tb->t + k*ntnr + j*nr ;
	    
	    sxm1_1d_cheb(nr, cf) ;	// division par (x-1)

	    
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    
    // base de developpement
    // inchangee
}
}
