/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


char op_sxpun_C[] = "$Header$" ;

/* 
 * Ensemble des routines de base de division par rapport a x+1
 * (Utilisation interne)
 * 
 *	void _sx_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */
 
 /*
 * $Id$
 * $Log$
 * Revision 1.4  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.3  2007/12/21 13:59:02  j_novak
 * Suppression of call to pow(-1, something).
 *
 * Revision 1.2  2007/12/20 09:11:08  jl_cornou
 * Correction of an error in op_sxpun about Jacobi(0,2) polynomials
 *
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
 * Revision 2.5  2000/09/07  12:50:48  phil
 * *** empty log message ***
 *
 * Revision 2.4  2000/02/24  16:42:59  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.3  1999/11/15  16:39:17  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.2  1999/10/18  13:40:11  eric
 * Suppression de la routine _sx_r_chebu car doublon avec _sxm1_cheb.
 *
 * Revision 2.1  1999/09/13  11:35:42  phil
 * gestion des cas k==1 et k==np
 *
 * Revision 2.0  1999/04/26  14:59:42  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

 // Fichier includes
#include "tbl.h"
#include <math.h>

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _sxpun_pas_prevu(Tbl * tb, int& base) {
    cout << "sxpun pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

			//-------------
			// Identite ---
			//-------------

void _sxpun_identite(Tbl* , int& ) {
    return ;
}

			//--------------
			// cas R_JACO02-
			//--------------

void _sxpun_r_jaco02(Tbl* tb, int&)
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

		xco[nr-1] = 0 ;
    		double somme ;
    		for (int i = 0 ; i < nr-1 ; i++) {
			somme = 0 ;
			for (int m = i+1 ; m < nr ; m++) {
			int signe = ((m-1-i)%2 == 0 ? 1 : -1) ; 
			somme += signe*((m+1)*(m+2)/double((i+1)*(i+2))-(i+1)*(i+2)/double((m+1)*(m+2)))*xci[m] ;
			}
		xco[i] = (2*i+3)/double(4)*somme ;
		}	
    
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

