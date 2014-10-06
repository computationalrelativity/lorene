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


char op_sxp12_C[] = "$Header$" ;

/* 
 * Ensemble des routines de base de division par rapport a (x+1)²
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
 * Revision 1.3  2014/10/06 15:16:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
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
#include <cmath>

void sxp12_1d(int , double **, int )	;

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _sxp12_pas_prevu(Tbl * tb, int& base) {
    cout << "sxp12 pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

			//-------------
			// Identite ---
			//-------------

void _sxp12_identite(Tbl* , int& ) {
    return ;
}

			//--------------
			// cas R_JACO02-
			//--------------

void _sxp12_r_jaco02(Tbl* tb, int& )
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

		double** tb1 = new double*[nr] ;
		for (int i = 0 ; i<nr ; i++ ) {
			*tb1[i]=xci[i] ;
			}
		sxp12_1d(nr,tb1,R_JACO02 >> TRA_R) ;
		for (int i = 0 ; i<nr ; i++ ) {
			xco[i] = *tb1[i] ;
		}
		delete [] (*tb1) ;		
	
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

