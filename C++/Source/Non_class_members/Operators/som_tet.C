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


char som_tet_C[] = "$Header$" ;

/*
 * Ensemble des routine pour la sommation directe en theta
 * 
 *   SYNOPSYS:
 *     double som_tet_XX
 *	(double* ti, int nt, int np, double tet, double* to)
 * 
 *   ATTENTION: np est le vrai nombre de points (sans le +2)
 *		on suppose que les tableaux tiennent compte de ca.
 */


/*
 * $Id$
 * $Log$
 * Revision 1.4  2002/05/11 12:36:54  e_gourgoulhon
 * Added basis T_COSSIN_SI.
 *
 * Revision 1.3  2002/05/05 16:20:40  e_gourgoulhon
 * Error message (in unknwon basis case) in English
 * Added the basis T_COSSIN_SP
 *
 * Revision 1.2  2002/05/01 07:41:05  e_gourgoulhon
 * Correction of an ERROR in som_tet_sin_p :
 *    sin(2*(j+1) * tet) --> sin(2*j * tet)
 * idem in som_tet_sin:
 *    sin( (j+1) * tet) --> sin(j * tet)
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  2000/09/08  16:26:32  eric
 * Ajout de la base T_SIN_I.
 *
 * Revision 2.4  2000/09/06  14:00:01  eric
 * Ajout de la base T_COS_I.
 *
 * Revision 2.3  2000/03/06  09:34:47  eric
 * Suppression des #include inutiles.
 *
 * Revision 2.2  1999/04/28  12:27:52  phil
 * Correction tailles des tableaux
 *
 * Revision 2.1  1999/04/26  14:31:17  phil
 * remplacement des malloc en new
 *
 * Revision 2.0  1999/04/12  15:42:03  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/04/12  15:41:20  phil
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
#include <math.h>


			//--------------------
			//- Cas Non-Prevu ---
			//------------------

void som_tet_pas_prevu
    (double* , const int, const int, const double, double*) {
	cout << "Mtbl_cf::val_point: theta basis not implemented yet ! "
	     << endl ;
	abort () ;
}

			//----------------
			//  Cas T_COS ---
			//----------------

void som_tet_cos
    (double* ti, const int nt, const int np,
    const double tet, double* to) {

// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* cosinus = new double [nt] ;
    for (j=0 ; j<nt ; j++) {
	cosinus[j] = cos(j * tet) ;
    }
    
    // Sommation sur le premier phi, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * cosinus[j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cosinus[j] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }

    }	// fin du cas np > 1 

    // Menage
    delete [] cosinus ;
    
}

			//-------------------
			//   Cas T_COS_P ---
			//------------------

void som_tet_cos_p
    (double* ti, const int nt, const int np,
    const double tet, double* to) {
    
// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* cosinus = new double [nt] ;
    for (j=0 ; j<nt ; j++) {
	cosinus[j] = cos(2*j * tet) ;
    }
    
    // Sommation sur le premier phi, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * cosinus[j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cosinus[j] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }

    }	// fin du cas np > 1 
    // Menage
   delete [] cosinus ;
    
}

			//----------------------
			//- Cas T_COS_I ---
			//---------------------

void som_tet_cos_i
    (double* ti, const int nt, const int np,
    const double tet, double* to) {
    
// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* cosinus = new double [nt] ;
    for (j=0 ; j<nt-1 ; j++) {
	cosinus[j] = cos((2*j+1) * tet) ;
    }
    cosinus[nt-1] = 0 ; 
    
    // Sommation sur le premier phi, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * cosinus[j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cosinus[j] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1 

    // Menage
    delete [] cosinus ;
    
}





			//---------------
			//  Cas T_SIN ---
			//---------------

void som_tet_sin
    (double* ti, const int nt, const int np,
    const double tet, double* to) {
    
// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* sinus = new double [nt] ;
    for (j=0 ; j<nt ; j++) {
	sinus[j] = sin(j * tet) ;
    }
    
    // Sommation sur le premier phi, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * sinus[j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * sinus[j] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1 

    // Menage
    delete [] sinus ;
    
}

			//------------------
			//  Cas T_SIN_P ---
			//-----------------

void som_tet_sin_p
    (double* ti, const int nt, const int np,
    const double tet, double* to) {
    
// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* sinus = new double [nt] ;
    for (j=0 ; j<nt-1 ; j++) {
	sinus[j] = sin(2*j * tet) ;
    }
    sinus[nt-1] = 0 ;

    // Sommation sur le premier phi, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * sinus[j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * sinus[j] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1 

    // Menage
   delete [] sinus ;
    
}

			//------------------
			//  Cas T_SIN_I ---
			//-----------------

void som_tet_sin_i
    (double* ti, const int nt, const int np,
    const double tet, double* to) {
    
// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* sinus = new double [nt] ;
    for (j=0 ; j<nt-1 ; j++) {
	sinus[j] = sin( (2*j+1) * tet) ;
    }
    sinus[nt-1] = 0 ; 
    
    // Sommation sur le premier phi, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * sinus[j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np)
    for (k=2 ; k<np+1 ; k++) {
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * sinus[j] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1 

    // Menage
   delete [] sinus ;
    
}


			//---------------------
			//  Cas T_COSSIN_CP ---
			//---------------------

void som_tet_cossin_cp
    (double* ti, const int nt, const int np,
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
    
    // Sommation sur le premier phi -> cosinus, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * cossin[2*j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np), suivant parite de m
    for (k=2 ; k<np+1 ; k++) {
	int m = (k/2) % 2 ;	    // parite: 0 <-> 1
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cossin[2*j + m] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1 

    // Menage
    delete [] cossin ;
    
}


			//----------------------
			//- Cas T_COSSIN_CI ---
			//---------------------

void som_tet_cossin_ci
    (double* ti, const int nt, const int np,
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
    
    // Sommation sur le premier phi -> cosinus, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * cossin[2*j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant
    
    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;
    
    // Sommation sur le reste des phi (pour k=2,...,np), suivant parite de m
    for (k=2 ; k<np+1 ; k++) {
	int m = (k/2) % 2 ;	    // parite: 0 <-> 1
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cossin[2*j + m] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1 

    // Menage
    delete [] cossin ;
    
}


			//---------------------
			//  Cas T_COSSIN_SP ---
			//---------------------

void som_tet_cossin_sp
    (double* ti, const int nt, const int np,
    const double tet, double* to) {

// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* cossin = new double [2*nt] ;
    for (j=0 ; j<2*nt ; j +=2) {
	cossin[j] = sin(j * tet) ;
	cossin[j+1] = cos((j+1) * tet) ;
    }

    // Sommation sur le premier phi -> cosinus, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * cossin[2*j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant

    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;

    // Sommation sur le reste des phi (pour k=2,...,np), suivant parite de m
    for (k=2 ; k<np+1 ; k++) {
	int m = (k/2) % 2 ;	    // parite: 0 <-> 1
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cossin[2*j + m] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1

    // Menage
    delete [] cossin ;

}


			//---------------------
			//  Cas T_COSSIN_SI ---
			//---------------------

void som_tet_cossin_si
    (double* ti, const int nt, const int np,
    const double tet, double* to) {

// Variables diverses
int j, k ;
double* pi = ti ;	    // Pointeur courant sur l'entree
double* po = to ;	    // Pointeur courant sur la sortie

    // Initialisation des tables trigo
    double* cossin = new double [2*nt] ;
    for (j=0 ; j<2*nt ; j +=2) {
	cossin[j] = sin((j+1) * tet) ;
	cossin[j+1] = cos(j * tet) ;
    }

    // Sommation sur le premier phi -> cosinus, k=0
    *po = 0 ;
    for (j=0 ; j<nt ; j++) {
	*po += (*pi) * cossin[2*j] ;
	pi++ ;	// Theta suivant
    }
    po++ ;	// Phi suivant

    if (np > 1) {	

    // On saute le phi suivant (sin(0)), k=1
    pi += nt ;
    po++ ;

    // Sommation sur le reste des phi (pour k=2,...,np), suivant parite de m
    for (k=2 ; k<np+1 ; k++) {
	int m = (k/2) % 2 ;	    // parite: 0 <-> 1
	(*po) = 0 ;
	for (j=0 ; j<nt ; j++) {
	    *po += (*pi) * cossin[2*j + m] ;
	    pi++ ;  // Theta suivant
	}
	po++ ;	    // Phi suivant
    }
    }	// fin du cas np > 1

    // Menage
    delete [] cossin ;

}



