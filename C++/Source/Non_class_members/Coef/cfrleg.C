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
 * Revision 1.4  2013/06/13 14:17:47  j_novak
 * Implementation of Legendre inverse coefficient transform.
 *
 * Revision 1.3  2013/06/07 14:44:34  j_novak
 * Coefficient computation for even Legendre basis.
 *
 * Revision 1.2  2013/06/06 15:31:32  j_novak
 * Functions to compute Legendre coefficients (not fully tested yet).
 *
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
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en r :    
    int nr = deg[2] ;
    int nm1 = nr - 1 ;

// Tests de dimension:
    if (nr > n3f) {
        cout << "cfrleg: nr > n3f : nr = " << nr << " ,  n3f = " 
        << n3f << endl ;
        abort () ;
        exit(-1) ;
    }
    if (nr > n3c) {
        cout << "cfrleg: nr > n3c : nr = " << nr << " ,  n3c = " 
        << n3c << endl ;
        abort () ;
        exit(-1) ;
    }
    if (n1f > n1c) {
        cout << "cfrleg: n1f > n1c : n1f = " << n1f << " ,  n1c = " 
        << n1c << endl ;
        abort () ;
        exit(-1) ;
    }
    if (n2f > n2c) {
        cout << "cfrleg: n2f > n2c : n2f = " << n2f << " ,  n2c = " 
        << n2c << endl ;
        abort () ;
        exit(-1) ;
    }

    Tbl* Pni = 0x0 ;
    Tbl* wn = 0x0 ;
    get_legendre_data(nr, Pni, wn) ;
    assert( (Pni != 0x0) && (wn != 0x0) ) ;
    double* cf_tmp = new double[nr] ;

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
	      cf_tmp[ii] = 0. ;
	      for (int jj = 0; jj<nr; jj++) 
		cf_tmp[ii] += ff0[jj] * (*wn)(jj) * (*Pni)(ii, jj) ;
	      cf_tmp[ii] /= double(2) / double(2*ii+1) ;
	    }
	    cf_tmp[nm1] /= double(nr+nm1) / double(nm1) ;
	    for (int i=0; i<nr; i++)
	      cf0[i] = cf_tmp[i] ;

	} 	// fin de la boucle sur theta 
    }	// fin de la boucle sur phi
    
    delete [] cf_tmp ;
}

void cfrlegp(const int* deg, const int* dimf, double* ff, const int* dimc,
	     double* cf)
  
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
  if (nr > n3f) {
    cout << "cfrlegp: nr > n3f : nr = " << nr << " ,  n3f = " 
	 << n3f << endl ;
    abort () ;
    exit(-1) ;
  }
  if (nr > n3c) {
    cout << "cfrlegp: nr > n3c : nr = " << nr << " ,  n3c = " 
	 << n3c << endl ;
    abort () ;
    exit(-1) ;
  }
  if (n1f > n1c) {
    cout << "cfrlegp: n1f > n1c : n1f = " << n1f << " ,  n1c = " 
	 << n1c << endl ;
    abort () ;
    exit(-1) ;
  }
  if (n2f > n2c) {
    cout << "cfrlegp: n2f > n2c : n2f = " << n2f << " ,  n2c = " 
	 << n2c << endl ;
    abort () ;
    exit(-1) ;
  }
  
  // Nombre de points:
  int nm1 = nr - 1;
  int dnm1 = 2*nr - 1 ;
  
  Tbl* Pni = 0x0 ;
  Tbl* wn = 0x0 ;
  get_legendre_data(dnm1, Pni, wn) ;
  double* cf_tmp = new double[nr] ;
  
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
	cf_tmp[ii] = 0.5*ff0[0]*(*Pni)(2*ii, nm1) 
	  * (*wn)(nm1) ;
	for (int jj=1; jj<nr; jj++) {
	  cf_tmp[ii] += ff0[jj]* (*wn)(nm1+jj) * (*Pni)(2*ii, nm1+jj) ;
	}
	cf_tmp[ii] *= double(4*ii+1) ;
      }
      cf_tmp[nm1] /= double(4*nm1+1) / double(2*nm1) ;
      for (int i=0; i<nr; i++)
	cf0[i] = cf_tmp[i] ;
      
    } 	// fin de la boucle sur theta 
  }	// fin de la boucle sur phi

  delete [] cf_tmp ;
  
}


void cfrlegi(const int* deg, const int* dimf, double* ff, const int* dimc,
	     double* cf)
  
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
  if (nr > n3f) {
    cout << "cfrlegi: nr > n3f : nr = " << nr << " ,  n3f = " 
	 << n3f << endl ;
    abort () ;
    exit(-1) ;
  }
  if (nr > n3c) {
    cout << "cfrlegi: nr > n3c : nr = " << nr << " ,  n3c = " 
	 << n3c << endl ;
    abort () ;
    exit(-1) ;
  }
  if (n1f > n1c) {
    cout << "cfrlegi: n1f > n1c : n1f = " << n1f << " ,  n1c = " 
	 << n1c << endl ;
    abort () ;
    exit(-1) ;
  }
  if (n2f > n2c) {
    cout << "cfrlegi: n2f > n2c : n2f = " << n2f << " ,  n2c = " 
	 << n2c << endl ;
    abort () ;
    exit(-1) ;
  }
  
  // Nombre de points:
  int nm1 = nr - 1;
  int dnm1 = 2*nr - 1 ;
  
  Tbl* Pni = 0x0 ;
  Tbl* wn = 0x0 ;
  get_legendre_data(dnm1, Pni, wn) ;
  double* cf_tmp = new double[nr] ;
  double* gam_tmp = new double[nr] ;
  
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
      
      for (int ii=0; ii<nr-1; ii++) {
	cf_tmp[ii] = 0.5*ff0[0]*(*Pni)(2*ii+1, nm1) 
	  * (*wn)(nm1) ;
	gam_tmp[ii] = 0. ;
	for (int jj=1; jj<nr; jj++) {
	  cf_tmp[ii] += ff0[jj]* (*wn)(nm1+jj) * (*Pni)(2*ii+1, nm1+jj) ;
	  gam_tmp[ii] += (*Pni)(2*ii+1, nm1+jj) * (*Pni)(2*ii+1, nm1+jj) * (*wn)(nm1+jj) ;
	}
	cf_tmp[ii] *= double(4*ii+3) ;
      }
      cf_tmp[nm1] = 0. ;
      for (int i=0; i<nr; i++)
	cf0[i] = cf_tmp[i] ;
      
    } 	// fin de la boucle sur theta 
  }	// fin de la boucle sur phi
  
  delete [] cf_tmp ;

}
