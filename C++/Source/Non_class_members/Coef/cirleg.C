/*
 *   Copyright (c) 2013 Jerome Novak
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


char cirleg_C[] = "$Header$" ;


/*
 * Transformation de Legendre inverse sur le troisieme indice 
 * (indice correspondant a r) d'un tableau 3-D.
 *
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:11  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2013/06/13 14:17:47  j_novak
 * Implementation of Legendre inverse coefficient transform.
 *
 * Revision 1.1  2013/06/06 15:31:33  j_novak
 * Functions to compute Legendre coefficients (not fully tested yet).
 *
 *
 * $Header $
 *
 */

// headers du C
#include <cassert>
#include <cstdlib>

//Lorene prototypes
#include "tbl.h"
#include "proto.h"
#include "utilitaires.h"


namespace Lorene {
//*****************************************************************************

void cirleg(const int* deg, const int* dimc, double* cf, const int* dimf,
		double* ff)

{

// Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en r :    
    int nr = deg[2] ;
    
// Tests de dimension:
    if (nr > n3c) {
	cout << "cirleg: nr > n3c : nr = " << nr << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nr > n3f) {
	cout << "cirleg: nr > n3f : nr = " << nr << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1c > n1f) {
	cout << "cirleg: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2c > n2f) {
	cout << "cirleg: n2c > n2f : n2c = " << n2c << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }

// boucle sur phi et theta

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1c = 1, on effectue la boucle une fois seulement.
 *    si n1c > 1, on va jusqu'a j = n1c-2 en sautant j = 1 (les coefficients
 *	j=n1c-1 et j=0 ne sont pas consideres car nuls). 
 */
    int borne_phi = ( n1c > 1 ) ? n1c-1 : 1 ;
    
    double* colloc = new double[nr] ;
    legendre_collocation_points(nr, colloc) ;
    
    for (int j=0; j< borne_phi; j++) {
      
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)

	for (int k=0; k<n2c; k++) {

	    int i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat

	    for (int i = 0; i<nr; i++) {
	      double x0 = colloc[i] ;
	      double Pi = 1. ;
	      double Pip1 = x0 ;
	      double som = cf0[0] + cf0[1]*x0 ;
	      for (int h=2; h<nr; h++) {
		double Pip2 = (2. - 1./double(h))*x0*Pip1
		  - (1. - 1./double(h))*Pi ;
		som += cf0[h]*Pip2 ;
		Pi = Pip1 ;
		Pip1 = Pip2 ;
	      }
	      ff0[i] = som ;
	    }
	} 	// fin de la boucle sur theta 
    }	// fin de la boucle sur phi
    delete [] colloc ;
}

//*****************************************************************************

void cirlegp(const int* deg, const int* dimc, double* cf, 
		    const int* dimf, double* ff)

{

// Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en r :    
    int nr = deg[2] ;
    
// Tests de dimension:
    if (nr > n3c) {
	cout << "cirlegp: nr > n3c : nr = " << nr << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nr > n3f) {
	cout << "cirlegp: nr > n3f : nr = " << nr << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1c > n1f) {
	cout << "cirlegp: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2c > n2f) {
	cout << "cirlegp: n2c > n2f : n2c = " << n2c << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }

// Nombre de points
    int nm1 = nr - 1;
    int dnm1 = 2*nr - 1 ;

// boucle sur phi et theta

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1c = 1, on effectue la boucle une fois seulement.
 *    si n1c > 1, on va jusqu'a j = n1c-2 en sautant j = 1 (les coefficients
 *	j=n1c-1 et j=0 ne sont pas consideres car nuls). 
 */
    int borne_phi = ( n1c > 1 ) ? n1c-1 : 1 ;

    double* colloc = new double[dnm1] ;
    legendre_collocation_points(dnm1, colloc) ;

    for (int j=0; j< borne_phi; j++) {
    
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)

	for (int k=0; k<n2c; k++) {

	    int i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat

	    for (int i = 0; i<nr; i++) {
	      double x0 = colloc[nm1+i] ;
	      double Pi = 1. ;
	      double Pip1 = x0 ;
	      double som = cf0[0] ;
	      for (int h=2; h<dnm1; h++) {
		double Pip2 = (2. - 1./double(h))*x0*Pip1
		  - (1. - 1./double(h))*Pi ;
		if (h%2 == 0) som += cf0[h/2]*Pip2 ;
		Pi = Pip1 ;
		Pip1 = Pip2 ;
	      }
	      ff0[i] = som ;
	    }

	} 	// fin de la boucle sur theta 
   }	// fin de la boucle sur phi

    delete [] colloc ;

}

//*****************************************************************************

void cirlegi(const int* deg, const int* dimc, double* cf, 
		    const int* dimf, double* ff)

{

// Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en r :    
    int nr = deg[2] ;
    
// Tests de dimension:
    if (nr > n3c) {
	cout << "cirlegi: nr > n3c : nr = " << nr << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nr > n3f) {
	cout << "cirlegi: nr > n3f : nr = " << nr << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1c > n1f) {
	cout << "cirlegi: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2c > n2f) {
	cout << "cirlegi: n2c > n2f : n2c = " << n2c << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }

// Nombre de points
    int nm1 = nr - 1;
    int dnm1 = 2*nr - 1 ;

// boucle sur phi et theta

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1c = 1, on effectue la boucle une fois seulement.
 *    si n1c > 1, on va jusqu'a j = n1c-2 en sautant j = 1 (les coefficients
 *	j=n1c-1 et j=0 ne sont pas consideres car nuls). 
 */
    int borne_phi = ( n1c > 1 ) ? n1c-1 : 1 ;

    double* colloc = new double[dnm1] ;
    legendre_collocation_points(dnm1, colloc) ;

    for (int j=0; j< borne_phi; j++) {
    
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)

	for (int k=0; k<n2c; k++) {

	    int i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat

	    for (int i = 0; i<nr; i++) {
	      double x0 = colloc[nm1+i] ;
	      double Pi = 1. ;
	      double Pip1 = x0 ;
	      double som = cf0[0]*x0 ;
	      for (int h=2; h<dnm1; h++) {
		double Pip2 = (2. - 1./double(h))*x0*Pip1
		  - (1. - 1./double(h))*Pi ;
		if (h%2 == 1) som += cf0[h/2]*Pip2 ;
		Pi = Pip1 ;
		Pip1 = Pip2 ;
	      }
	      ff0[i] = som ;
	    }

	} 	// fin de la boucle sur theta 
   }	// fin de la boucle sur phi

    delete [] colloc ;

}

}
