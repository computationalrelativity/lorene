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


char som_r_C[] = "$Header$" ;

/*
 * Ensemble des routine pour la sommation directe en r
 * 
 *   SYNOPSYS:
 *     double som_r_XX
 *	(double* ti, int nr, int nt, int np, double x, double* trtp)
 *
 *     x est l'argument du polynome de Chebychev: x in [0, 1] ou x in [-1, 1]
 * 
 *   ATTENTION: np est suppose etre le nombre de points (sans le +2)
 *		on suppose que trtp tient compte de ca.
 */


/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.4  2000/03/06  09:34:21  eric
 * Suppression des #include inutiles.
 *
 * Revision 2.3  1999/04/28  12:11:15  phil
 * changements de sommations
 *
 * Revision 2.2  1999/04/26  14:26:31  phil
 * remplacement des malloc en new
 *
 * Revision 2.1  1999/04/13  14:44:06  phil
 * ajout de som_r_chebi
 *
 * Revision 2.0  1999/04/12  15:42:33  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/04/12  15:40:25  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C++
#include <iostream.h>

// Headers C
#include <stdlib.h>

			//-------------------
			//- Cas Non-Prevu ---
			//-------------------

void som_r_pas_prevu
    (double* ti, const int nr, const int nt, const int np, const double x, double* trtp) {
	cout << "Sommation en r sur une base non prevue" << endl ;
	cout << ti << " " << nr << " " << nt << " " << np << " " 
	    << x << " " << trtp << endl ;
	abort () ;
    }

			//----------------
			//- Cas R_CHEB ---
			//----------------

void som_r_cheb
    (double* ti, const int nr, const int nt, const int np, const double x, 
    double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Chebyshev au point x demande
    double* cheb = new double [nr] ;
    cheb[0] = 1. ;
    cheb[1] = x ;
    for (i=2; i<nr; i++) {
	cheb[i] = 2*x* cheb[i-1] - cheb[i-2] ;	    
    }
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[i] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }
    
    }	// fin du cas np > 1 

    // Menage
    delete [] cheb ;

}


			//-----------------
			//- Cas R_CHEBP ---
			//-----------------

void som_r_chebp
    (double* ti, const int nr, const int nt, const int np, const double x, 
    double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Chebyshev au point x demande
    double* cheb = new double [nr] ;
    cheb[0] = 1. ;
    double t2im1 = x ;
    for (i=1; i<nr; i++) {
	cheb[i] = 2*x* t2im1 - cheb[i-1] ;
	t2im1 = 2*x* cheb[i] - t2im1 ;	    
    }
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[i] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }

    }	// fin du cas np > 1 
    // Menage
    delete [] cheb ;

}


			//-----------------
			//- Cas R_CHEBI ---
			//-----------------

void som_r_chebi
    (double* ti, const int nr, const int nt, const int np, const double x, 
    double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Chebyshev au point x demande
    double* cheb = new double [nr] ;
    cheb[0] = x ;
    double t2im1 = 1. ;
    for (i=1; i<nr; i++) {
	t2im1 = 2*x* cheb[i-1] - t2im1 ;
	cheb[i] = 2*x* t2im1 - cheb[i-1] ;	    
    }
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[i] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }

    }	// fin du cas np > 1 
    // Menage
   delete [] cheb ;

}
			//-----------------
			//- Cas R_CHEBU ---
			//-----------------

void som_r_chebu
    (double* ti, const int nr, const int nt, const int np, const double x, double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Chebyshev au point x demande
    double* cheb = new double [nr] ;
    cheb[0] = 1. ;
    cheb[1] = x ;
    for (i=2; i<nr; i++) {
	cheb[i] = 2*x* cheb[i-1] - cheb[i-2] ;	    
    }
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[i] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }

    }	// fin du cas np > 1 
    // Menage
    delete [] cheb ;
}

			//----------------------
			//  Cas R_CHEBPIM_P ---
			//---------------------

void som_r_chebpim_p
    (double* ti, const int nr, const int nt, const int np, const double x, double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Chebyshev au point x demande
    double* cheb = new double [2*nr] ;
    cheb[0] = 1. ;
    cheb[1] = x ;
    for (i=2 ; i<2*nr ; i++) {
	cheb[i] = 2*x* cheb[i-1] - cheb[i-2] ;	    
    }
    
    // Sommation pour le premier phi, k=0
    int m = 0;
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[2*i] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2 ;		    // parite: 0 <-> 1
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[2*i + m] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }
    
    }	// fin du cas np > 1 
    // Menage
   delete [] cheb ;
    
}

			//----------------------
			//- Cas R_CHEBPIM_I ----
			//----------------------

void som_r_chebpim_i
    (double* ti, const int nr, const int nt, const int np, const double x, double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Chebyshev au point x demande
    double* cheb = new double [2*nr] ;
    cheb[0] = 1. ;
    cheb[1] = x ;
    for (i=2 ; i<2*nr ; i++) {
	cheb[i] = 2*x* cheb[i-1] - cheb[i-2] ;	    
    }
    
    // Sommation pour le premier phi, k=0
    int m = 0;
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[2*i+1] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2  ;		    // parite: 0 <-> 1
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[2*i + 1 - m] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }
    
    }	// fin du cas np > 1 
    // Menage
    delete [] cheb ;
    
}

