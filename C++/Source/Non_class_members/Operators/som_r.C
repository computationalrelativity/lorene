/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
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
 * Revision 1.11  2014/10/13 08:53:27  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:16:07  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2013/06/13 14:18:18  j_novak
 * Inclusion of new bases R_LEG, R_LEGP and R_LEGI.
 *
 * Revision 1.8  2008/08/27 08:50:10  jl_cornou
 * Added Jacobi(0,2) polynomials
 *
 * Revision 1.7  2007/12/11 15:28:18  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.6  2006/05/17 13:15:19  j_novak
 * Added a test for the pure angular grid case (nr=1), in the shell.
 *
 * Revision 1.5  2005/02/18 13:14:17  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.4  2004/11/23 15:16:02  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.3  2002/10/16 14:36:59  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/05/05 16:20:40  e_gourgoulhon
 * Error message (in unknwon basis case) in English
 * Added the basis T_COSSIN_SP
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
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
// Headers C
#include <cstdlib>

#include "headcpp.h"
#include "proto.h"

namespace Lorene {
			//-------------------
			//- Cas Non-Prevu ---
			//-------------------

void som_r_pas_prevu
    (double*, const int, const int, const int, const double, double*) {
	cout << "Mtbl_cf::val_point: r basis not implemented yet !"
	     << endl ;
	abort() ;
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
    if (nr > 1) cheb[1] = x ;
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

			//----------------------
			//  Cas R_CHEBPI_P ---
			//---------------------

void som_r_chebpi_p
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
    for (j=0 ; j<nt ; j++) {
      int l = j%2;
      if(l==0){
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	  *po += (*pi) * cheb[2*i] ;
	  pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
      } else {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[2*i+1] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
      }	
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	for (j=0 ; j<nt ; j++) {
	  int l = j% 2 ;		    // parite: 0 <-> 1
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[2*i + l] ;
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
			//- Cas R_CHEBPI_I ----
			//----------------------

void som_r_chebpi_i
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
    for (j=0 ; j<nt ; j++) {
      int l = j%2 ;
      if(l==1){
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	  *po += (*pi) * cheb[2*i] ;
	  pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
      } else {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * cheb[2*i+1] ;
	    pi++ ;  // R suivant
	}
	po++ ;	    // Theta suivant
      }	
    }
    
    if (np > 1) {	

    // On saute le deuxieme phi (sin(0)), k=1
    pi += nr*nt ;
    po += nt ;
    
    // Sommation sur les phi suivants (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
      for (j=0 ; j<nt ; j++) {
	  int l = j % 2  ;		    // parite: 0 <-> 1
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[2*i + 1 - l] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }
    
    }	// fin du cas np > 1 
    // Menage
    delete [] cheb ;
    
}
			//----------------
			//- Cas R_LEG  ---
			//----------------

void som_r_leg
    (double* ti, const int nr, const int nt, const int np, const double x, 
    double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Legendre au point x demande
    double* leg = new double [nr] ;
    leg[0] = 1. ;
    if (nr > 1) leg[1] = x ;
    for (i=2; i<nr; i++) {
      leg[i] = (double(2*i-1)*x* leg[i-1] - double(i-1)*leg[i-2]) / double(i) ; 
    }
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * leg[i] ;
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
		*po += (*pi) * leg[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }
    
    }	// fin du cas np > 1 

    // Menage
    delete [] leg ;

}


			//-----------------
			//- Cas R_LEGP ---
			//-----------------

void som_r_legp
    (double* ti, const int nr, const int nt, const int np, const double x, 
    double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Legendre au point x demande
    double* leg = new double [nr] ;
    leg[0] = 1. ;
    double p2im1 = x ;
    for (i=1; i<nr; i++) {
      leg[i] = (double(4*i-1)*x* p2im1 - double(2*i-1)*leg[i-1]) / double(2*i) ;
      p2im1 = (double(4*i+1)*x* leg[i] - double(2*i)*p2im1) / double(2*i+1) ;	    
    }
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * leg[i] ;
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
		*po += (*pi) * leg[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }

    }	// fin du cas np > 1 
    // Menage
    delete [] leg ;

}


			//-----------------
			//- Cas R_LEGI ---
			//-----------------

void som_r_legi
    (double* ti, const int nr, const int nt, const int np, const double x, 
    double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Legendre au point x demande
    double* leg = new double [nr] ;
    leg[0] = x ;
    double p2im1 = 1. ;
    for (i=1; i<nr; i++) {
      p2im1 = (double(4*i-1)*x* leg[i-1] - double(2*i-1)*p2im1) / double(2*i) ;
      leg[i] = (double(4*i+1)*x* p2im1 - double(2*i)*leg[i-1]) / double(2*i+1) ;
    }
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * leg[i] ;
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
		*po += (*pi) * leg[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }

    }	// fin du cas np > 1 
    // Menage
   delete [] leg ;

}

			//----------------
			//- Cas R_JACO02 -
			//----------------

void som_r_jaco02
    (double* ti, const int nr, const int nt, const int np, const double x, 
    double* trtp) {

// Variables diverses
int i, j, k ;
double* pi = ti ;	    // pointeur courant sur l'entree
double* po = trtp ;	    // pointeur courant sur la sortie
    
    // Valeurs des polynomes de Jacobi(0,2) au point x demande
    double* jaco = jacobi(nr-1,x) ;
    
    // Sommation pour le premier phi, k=0
    for (j=0 ; j<nt ; j++) {
	*po = 0 ;
	// Sommation sur r
	for (i=0 ; i<nr ; i++) {
	    *po += (*pi) * jaco[i] ;
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
		*po += (*pi) * jaco[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
    }
    
    }	// fin du cas np > 1 

    // Menage
    delete [] jaco ;

}
}
