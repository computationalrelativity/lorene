/*
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


char op_mult_xm1_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.2  1999/11/15  16:38:48  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.1  1999/04/26  16:19:38  phil
 * changement commentaire
 *
 * Revision 2.0  1999/04/26  16:16:24  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

/* 
 * Ensemble des routines de base pour l'operateur (x-1) Id
 * (Utilisation interne)
 * 
 *  Prototype :
 *  ---------
 *	void _mult_xm1_XXXX(Tbl * tb, int& base)
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
 
 
 //Lorene
#include "tbl.h"
#include "type_parite.h"

// Prototypage
#include "proto.h"


		//------------------------------------
		// Routine qui ne fait rien	   ---
		//------------------------------------

void _mult_xm1_identite(Tbl* , int& ) {
    return ; 
}

			//-------------------------
			// cas R_CHEB et R_CHEBU --
			//-------------------------

void _mult_xm1_cheb(Tbl *tb, int& )
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
    
    double* trav = new double[nr] ;
    
    int k, j, i ; 
    for (k=0 ; k<np+1 ; k++) {
	if (k==1) continue ;	// On ne traite pas le coefficient de sin(0*phi) 
	for (j=0 ; j<nt ; j++) {
	
	    double* cf = tb->t + k*ntnr + j*nr ;
	    
	    mult_xm1_1d_cheb(nr, cf, trav) ;	// multiplication par (x-1)

	    for (i=0; i<nr; i++) {
		cf[i] = trav[i] ; 
	    }
	    
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
   delete [] trav ; 
    
    // base de developpement
    // inchangee
}

