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


char cfrleg_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.1  2013/06/05 15:08:13  j_novak
 * Initial revision. Not ready yet...
 *
 *
 *
 * $Header$
 *
 */


// headers du C
#include <cstdlib>
#include <cassert>

//Lorene prototypes
#include "tbl.h"
#include "utilitaires.h"

void get_legendre_data(int, Tbl*&, Tbl*& ) ;

//*****************************************************************************

void cfrleg(const int* deg, const int* dimf, double* ff, const int* dimc, 
	    double* cf)

{
  // Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en r :    
    int nr = deg[2] ;
    int nm1 = nr - 1 ;

    Tbl* Pni = 0x0 ;
    Tbl* wn = 0x0 ;
    get_legendre_data(nr, Pni, wn) ;
    assert( (Pni != 0x0) && (wn != 0x0) ) ;

    // boucle sur phi et theta

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1f = 1, on effectue la boucle une fois seulement.
 *    si n1f > 1, on va jusqu'a j = n1f-2 en sautant j = 1 (les coefficients
 *	j=n1f-1 et j=0 ne sont pas consideres car nuls). 
 */
    int borne_phi = ( n1f > 1 ) ? n1f-1 : 1 ;

    for (int j=0; j< borne_phi; j++) {
    
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)

	for (int k=0; k<n2f; k++) {

	    int i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau resultat

	    for (int ii=0; ii<nr; ii++) {
	      cf0[ii] = 0. ;
	      for (int jj = 0; jj<nr; jj++) 
		cf0[ii] += ff0[jj] * (*wn)(jj) * (*Pni)(ii, jj) ;
	      cf0[ii] /= double(2) / double(2*ii+1) ;
	    }
	    cf0[nm1] /= double(nr+nm1) / double(nm1) ;

	} 	// fin de la boucle sur theta 
    }	// fin de la boucle sur phi
    
    
}

void cfrlegp(const int*, const int*, double*, const int*, double*)

{
  c_est_pas_fait("cfrlegp") ;
}

void cfrlegi(const int*, const int*, double*, const int*, double*)

{
  c_est_pas_fait("cfrlegi") ;
}
