/*
 *   Copyright (c) 1999-2001 Jerome Novak
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


char op_mult_xp1_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/12/11 15:42:23  jl_cornou
 * Premiere version des fonctions liees aux polynomes de Jacobi(0,2)
 *
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/09/07  12:49:53  phil
 * *** empty log message ***
 *
 * Revision 1.2  2000/02/24  16:42:18  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 1.1  1999/11/16  13:37:41  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

/* 
 * Ensemble des routines de base de multiplication par x+1
 * (Utilisation interne)
 * 
 *	void _mult_x_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */
 
 // Fichier includes
#include "tbl.h"


		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _mult_xp1_pas_prevu(Tbl * tb, int& base) {
    cout << "mult_xp1 pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

			//-------------
			// Identite ---
			//-------------

void _mult_xp1_identite(Tbl* , int& ) {
    return ;
}

			//---------------
			// cas R_JACO02 -
			//---------------

void _mult_xp1_r_jaco02(Tbl* tb, int& base)
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
    
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants

    int borne_phi = np + 1 ; 
    if (np == 1) {
	borne_phi = 1 ; 
    }
    
    for (int k=0 ; k< borne_phi ; k++)
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	}
	else {
	for (int j=0 ; j<nt ; j++) {

    		xco[0] = 1.5*xci[0] + 0.3*xci[1] ;
    		for (int i = 1 ; i < nr-1 ; i++) {
		xco[i] = i*(i+2)/double((i+1)*(2*i+1))*xci[i-1] + (i*i+3*i+3)/double((i+1)*(i+2))*xci[i] + (i+1)*(i+3)/double((i+2)*(2*i+5))*xci[i+1] ;
    		}
    		xco[nr-1] = (nr*nr-1)/double((nr)*(2*nr-1))*xci[nr-2] + (1+1/double((nr)*(nr+1)))*xci[nr-1] ;
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangée

}
