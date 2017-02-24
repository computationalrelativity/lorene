/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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
 * Ensemble des routine pour la sommation directe en r, theta et phi dans
 * le cas symetrie equatoriale (plan z=0) + antisymetrie par rapport au plan y=0
 * 
 *   SYNOPSYS:
 *     double som_r_XX_asymy
 *	(double* ti, int nr, int nt, int np, double x, double* trtp)
 *
 *     x est l'argument du polynome de Chebychev: x in [0, 1] ou x in [-1, 1]
 *
 *
 *     double som_tet_XX_asymy
 *	(double* ti, int nt, int np, double tet, double* to)
 * 
 *     double som_phi_XX_asymy
 *	(double* ti, int np, double phi, double* xo)
 *
 *   ATTENTION: np est suppose etre le nombre de points (sans le +2)
 *		on suppose que trtp tient compte de ca.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.5  2017/02/24 15:34:59  j_novak
 * Removal of spurious comments
 *
 * Revision 1.4  2016/12/05 16:18:08  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:26  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:16:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/03/06  10:27:45  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers C
#include <cmath>

namespace Lorene {

//****************************************************************************
//			    Sommation en r
//****************************************************************************

			//////////////////
			//  Cas R_CHEB  //
			//////////////////

void som_r_cheb_asymy(double* ti, const int nr, const int nt, const int np, 
		      const double x, double* trtp) {

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
    
    // On saute les 3 premiers coef. en phi, qui sont respectivemment
    //  cos(0 phi), sin(0 phi) et cos(phi)
    pi += 3*nr*nt ;
    po += 3*nt ;
        
    // Sommation sur les phi suivants (pour k=3,...,np)
    //  en sautant les cosinus (d'ou le k+=2)
    for (k=3 ; k<np+1 ; k+=2) {
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
	// On saute le cos(m*phi) : 
	pi += nr*nt ;
	po += nt ;
    }
    
    // Menage
    delete [] cheb ;
}



			///////////////////
			//  Cas R_CHEBU  //
			///////////////////

void som_r_chebu_asymy(double* ti, const int nr, const int nt, const int np, 
		       const double x, double* trtp) {

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
    
    // On saute les 3 premiers coef. en phi, qui sont respectivemment
    //  cos(0 phi), sin(0 phi) et cos(phi)
    pi += 3*nr*nt ;
    po += 3*nt ;
    
    // Sommation sur les phi suivants (pour k=3,...,np)
    //  en sautant les cosinus (d'ou le k+=2)
    for (k=3 ; k<np+1 ; k+=2) {
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[i] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
	// On saute le cos(m*phi) : 
	pi += nr*nt ;
	po += nt ;
    }

    // Menage
    delete [] cheb ;
    
}



			///////////////////////
			//  Cas R_CHEBPIM_P  //
			///////////////////////

void som_r_chebpim_p_asymy(double* ti, const int nr, const int nt, 
			   const int np, const double x, double* trtp) {

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
    

    // On saute les 3 premiers coef. en phi, qui sont respectivemment
    //  cos(0 phi), sin(0 phi) et cos(phi)
    pi += 3*nr*nt ;
    po += 3*nt ;
    
    // Sommation sur les phi suivants (pour k=3,...,np)
    //  en sautant les cosinus (d'ou le k+=2)
    for (k=3 ; k<np+1 ; k+=2) {
	int m = (k/2) % 2 ;		    // parite: 0 <-> 1
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[2*i + m] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
	// On saute le cos(m*phi) : 
	pi += nr*nt ;
	po += nt ;
    }
    

    // Menage
   delete [] cheb ;
    
}

			///////////////////////
			//  Cas R_CHEBPIM_I  //
			///////////////////////

void som_r_chebpim_i_asymy(double* ti, const int nr, const int nt, 
			   const int np, const double x, double* trtp) {

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
    
    // On saute les 3 premiers coef. en phi, qui sont respectivemment
    //  cos(0 phi), sin(0 phi) et cos(phi)
    pi += 3*nr*nt ;
    po += 3*nt ;
        
    // Sommation sur les phi suivants (pour k=3,...,np)
    //  en sautant les cosinus (d'ou le k+=2)
    for (k=3 ; k<np+1 ; k+=2) {
	int m = (k/2) % 2 ;		    // parite: 0 <-> 1
	for (j=0 ; j<nt ; j++) {
	    // Sommation sur r
	    *po = 0 ;
	    for (i=0 ; i<nr ; i++) {
		*po += (*pi) * cheb[2*i + 1 - m] ;
		pi++ ;	// R suivant
	    }
	    po++ ;	// Theta suivant
	}
	// On saute le cos(m*phi) : 
	pi += nr*nt ;
	po += nt ;
    }
    
    // Menage
    delete [] cheb ;
    
}



//****************************************************************************
//			    Sommation en theta
//****************************************************************************

			///////////////////////
			//  Cas T_COSSIN_CP  //
			///////////////////////

void som_tet_cossin_cp_asymy(double* ti, const int nt, const int np,
			     const double tet, double* to) {
    
// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* cossin = new double [2*nt] ;
    for (j=0 ; j<2*nt ; j +=2) {
	cossin[j] = cos(j * tet) ;
	cossin[j+1] = sin((j+1) * tet) ;
    }
    
    // On saute les 3 premiers coef. en phi, qui sont respectivemment
    //  cos(0 phi), sin(0 phi) et cos(phi)
    pi += 3*nt ;
    po += 3 ;

    // Sommation sur le reste des phi (pour k=3,...,np), suivant parite de m
    for (k=3 ; k<np+1 ; k+=2) {
	int m = (k/2) % 2 ;	    // parite: 0 <-> 1
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cossin[2*j + m] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant

	// On saute le cos(m*phi) : 
	pi += nt ;
	po++ ;

    }

    // Menage
    delete [] cossin ;
}

			///////////////////////
			//  Cas T_COSSIN_CI  //
			///////////////////////

void som_tet_cossin_ci_asymy(double* ti, const int nt, const int np,
			     const double tet, double* to) {
    
// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* cossin = new double [2*nt] ;
    for (j=0 ; j<2*nt ; j +=2) {
	cossin[j] = cos((j+1) * tet) ;
	cossin[j+1] = sin(j * tet) ;
    }
    
    // On saute les 3 premiers coef. en phi, qui sont respectivemment
    //  cos(0 phi), sin(0 phi) et cos(phi)
    pi += 3*nt ;
    po += 3 ;

    // Sommation sur le reste des phi (pour k=3,...,np), suivant parite de m
    for (k=3 ; k<np+1 ; k+=2) {
	int m = (k/2) % 2 ;	    // parite: 0 <-> 1
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cossin[2*j + m] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant

	// On saute le cos(m*phi) : 
	pi += nt ;
	po++ ;

    }

    // Menage
    delete [] cossin ;
}


//****************************************************************************
//			    Sommation en phi
//****************************************************************************

void som_phi_cossin_asymy(double* ti, const int np, const double phi, 
			  double* xo) {
    
    *xo = 0 ; 

    // Sommation sur les sinus seulement
    for (int k=3 ; k<np+1 ; k +=2 ) {
	int m = k/2 ;
	*xo += ti[k] * sin(m * phi) ;
    }

}

}
