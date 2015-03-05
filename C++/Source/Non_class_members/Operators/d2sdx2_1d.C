/*
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


char d2sdx2_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2015/03/05 08:49:31  j_novak
 * Implemented operators with Legendre bases.
 *
 * Revision 1.6  2014/10/13 08:53:23  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:16:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2007/12/11 15:28:18  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.3  2002/10/16 15:05:54  j_novak
 * *** empty log message ***
 *
 * Revision 1.2  2002/10/16 14:36:58  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  1999/10/11  15:18:14  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/10/11  14:27:06  phil
 * initialisation des variables statiques
 *
 * Revision 2.2  1999/09/03  14:15:56  phil
 * Correction termes 0 (/2)
 *
 * Revision 2.1  1999/07/08  09:53:13  phil
 * correction gestion memoire
 *
 * Revision 2.0  1999/07/07  10:15:26  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

#include <cstdlib>
#include "type_parite.h"
#include "headcpp.h"
#include "proto.h"

/*
 * Routine appliquant l'operateur d2sdx2.
 * 
 * Entree : tb contient les coefficients du developpement 
 *	    int nr : nombre de points en r. 
 *
 * Sortie : tb contient d2sdx2
 * 
 */

		//----------------------------------
		// Routine pour les cas non prevus --
		//----------------------------------

namespace Lorene {
void _d2sdx2_1d_pas_prevu(int nr, double* tb, double *xo) {
    cout << "d2sdx2 pas prevu..." << endl ;
    cout << "Nombre de points : " << nr << endl ;
    cout << "Valeurs : " << tb << "  " << xo <<endl ;
    abort() ;
    exit(-1) ;
}

			//----------------
			// cas R_CHEBU ---
			//----------------

void _d2sdx2_1d_r_chebu(int nr,  double* tb, double *xo)
{
    // Variables statiques
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  (i+2)*(i+2)*(i+2) ;
	    cx2[i] =  (i+2) ;
	    cx3[i] =  i*i ;
	    }
	}

    double som1, som2 ;
	    
    xo[nr-1] = 0 ;
    som1 = (nr-1)*(nr-1)*(nr-1) * tb[nr-1] ;
    som2 = (nr-1) * tb[nr-1] ;
    xo[nr-3] = som1 - (nr-3)*(nr-3)*som2 ;
    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
	som1 += cx1[i] * tb[i+2] ;
	som2 += cx2[i] * tb[i+2] ;
	xo[i] = som1 - cx3[i] * som2 ;
    }	// Fin de la premiere boucle sur r
	    
    xo[nr-2] = 0 ;
    som1 = (nr-2)*(nr-2)*(nr-2) * tb[nr-2] ;
    som2 = (nr-2) * tb[nr-2] ;
    xo[nr-4] = som1 - (nr-4)*(nr-4)*som2 ;
    for (int i = nr-6 ; i >= 0 ; i -= 2 ) {
	som1 += cx1[i] * tb[i+2] ;
	som2 += cx2[i] * tb[i+2] ;
	xo[i] = som1 - cx3[i] * som2 ;
    }	    // Fin de la deuxieme boucle sur r
	    xo[0] *= 0.5 ;
	
}

			//---------------
			// cas R_CHEB ---
			//---------------

void _d2sdx2_1d_r_cheb(int nr, double* tb, double *xo)
{
  
    // Variables statiques
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  (i+2)*(i+2)*(i+2) ;
	    cx2[i] =  (i+2) ;
	    cx3[i] =  i*i ;
	    }
	}

    double som1, som2 ;
	    
    xo[nr-1] = 0 ;
    som1 = (nr-1)*(nr-1)*(nr-1) * tb[nr-1] ;
    som2 = (nr-1) * tb[nr-1] ;
    xo[nr-3] = som1 - (nr-3)*(nr-3)*som2 ;
    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
	som1 += cx1[i] * tb[i+2] ;
	som2 += cx2[i] * tb[i+2] ;
	xo[i] = som1 - cx3[i] * som2 ;
    }	// Fin de la premiere boucle sur r
    xo[nr-2] = 0 ;
    som1 = (nr-2)*(nr-2)*(nr-2) * tb[nr-2] ;
    som2 = (nr-2) * tb[nr-2] ;
    xo[nr-4] = som1 - (nr-4)*(nr-4)*som2 ;
    for (int i = nr-6 ; i >= 0 ; i -= 2 ) {
	som1 += cx1[i] * tb[i+2] ;
	som2 += cx2[i] * tb[i+2] ;
	xo[i] = som1 - cx3[i] * som2 ;
    }	// Fin de la deuxieme boucle sur r
    xo[0] *= .5 ;
	   
}

			  //----------------
			 // cas R_JACO02 --
			//----------------

void _d2sdx2_1d_r_jaco02(int nr, double* tb, double *xo)
{
	dsdx_1d(nr, &tb, R_JACO02 >> TRA_R) ;
	dsdx_1d(nr, &tb, R_JACO02 >> TRA_R) ;
	for (int i = 0 ;  i<nr ; i++) {
		xo[i] = tb[i] ;
	}
}


			//-----------------
			// cas R_CHEBP ---
			//----------------

void _d2sdx2_1d_r_chebp(int nr, double* tb, double *xo)
{
    // Variables statiques
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  8*(i+1)*(i+1)*(i+1) ;
	    cx2[i] =  2*(i+1) ;
	    cx3[i] =  4*i*i ;
	    }
	}
    
    double som1, som2 ;

    xo[nr-1] = 0 ; 
    som1 = 8*(nr-1)*(nr-1)*(nr-1) * tb[nr-1] ;
    som2 = 2*(nr-1) * tb[nr-1] ;
    xo[nr-2] = som1 - 4*(nr-2)*(nr-2)*som2 ;
  
    for (int i = nr-3 ; i >= 0 ; i-- ) {
	som1 += cx1[i] * tb[i+1] ;
	som2 += cx2[i] * tb[i+1] ;
	xo[i] = som1 - cx3[i] * som2 ;
    }	// Fin de la boucle sur r
    xo[0] *= .5 ;
}

			//---------------
			// cas R_CHEBI --
			//---------------

void _d2sdx2_1d_r_chebi(int nr, double* tb, double *xo)
{
    // Variables statiques
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  (2*i+3)*(2*i+3)*(2*i+3) ;
	    cx2[i] =  (2*i+3) ;
	    cx3[i] =  (2*i+1)*(2*i+1) ;
	    }
	}
  
    // pt. sur le tableau de double resultat
    double som1, som2 ;
	    
    xo[nr-1] = 0 ;
    som1 = (2*nr-1)*(2*nr-1)*(2*nr-1) * tb[nr-1] ;
    som2 = (2*nr-1) * tb[nr-1] ;
    xo[nr-2] = som1 - (2*nr-3)*(2*nr-3)*som2 ;
    for (int i = nr-3 ; i >= 0 ; i-- ) {
	som1 += cx1[i] * tb[i+1] ;
	som2 += cx2[i] * tb[i+1] ;
	xo[i] = som1 - cx3[i] * som2 ;
    }	// Fin de la boucle su r

}

			//--------------
			// cas R_LEG ---
			//--------------

void _d2sdx2_1d_r_leg(int nr, double* tb, double *xo)
{
  
    // Variables statiques
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;
	for (int i=0 ; i<nr ; i++) {
	  cx1[i] =  (i+2)*(i+3) ;
	  cx2[i] =  i*(i+1) ;
	  cx3[i] =  double(i) + 0.5 ;
	    }
	}

    double som1, som2 ;
	    
    xo[nr-1] = 0 ;
    som1 = (nr-1)* nr * tb[nr-1] ;
    som2 = tb[nr-1] ;
    if (nr > 2) 
      xo[nr-3] = (double(nr) - 2.5) * (som1 - (nr-3)*(nr-2)*som2) ;
    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
	som1 += cx1[i] * tb[i+2] ;
	som2 += tb[i+2] ;
	xo[i] = cx3[i]*(som1 - cx2[i] * som2) ;
    }	// Fin de la premiere boucle sur r
    if (nr > 1) xo[nr-2] = 0 ;
    if (nr > 3) {
      som1 = (nr-2)*(nr-1) * tb[nr-2] ;
      som2 = tb[nr-2] ;
      xo[nr-4] = (double(nr) - 3.5) * (som1 - (nr-4)*(nr-3)*som2) ;
    }
    for (int i = nr-6 ; i >= 0 ; i -= 2 ) {
	som1 += cx1[i] * tb[i+2] ;
	som2 += tb[i+2] ;
	xo[i] = cx3[i]*(som1 - cx2[i] * som2) ;
    }	// Fin de la deuxieme boucle sur r
	   
}

			//---------------
			// cas R_LEGP ---
			//---------------

void _d2sdx2_1d_r_legp(int nr, double* tb, double *xo)
{
    // Variables statiques
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
      nr_pre = nr ;
      if (cx1 != 0x0) delete [] cx1 ;
      if (cx2 != 0x0) delete [] cx2 ;
      if (cx3 != 0x0) delete [] cx3 ;
      cx1 = new double [nr] ;
      cx2 = new double [nr] ;
      cx3 = new double [nr] ;
      for (int i=0 ; i<nr ; i++) {
	cx1[i] =  (2*i+2)*(2*i+3) ;
	cx2[i] =  2*i*2*(i+1) ;
	cx3[i] =  double(2*i)+ 0.5 ;
      }
    }
    
    double som1, som2 ;

    xo[nr-1] = 0 ; 
    som1 = (2*nr-2)*(2*nr-1)* tb[nr-1] ;
    som2 = tb[nr-1] ;
    if (nr > 1) 
      xo[nr-2] = (double(2*nr) - 1.5)*(som1 - 2*(nr-2)*(2*nr-1)*som2) ;
    for (int i = nr-3 ; i >= 0 ; i-- ) {
      som1 += cx1[i] * tb[i+1] ;
      som2 += tb[i+1] ;
      xo[i] = cx3[i]*(som1 - cx2[i]*som2) ;
    }	// Fin de la boucle sur r
}

			//---------------
			// cas R_LEGI --
			//---------------

void _d2sdx2_1d_r_legi(int nr, double* tb, double *xo)
{
    // Variables statiques
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  (2*i+3)*(2*i+4) ;
	    cx2[i] =  (2*i+1)*(2*i+2) ;
	    cx3[i] =  double(2*i) + 1.5 ;
	    }
	}
  
    // pt. sur le tableau de double resultat
    double som1, som2 ;
	    
    xo[nr-1] = 0 ;
    som1 = (2*nr-1)*(2*nr) * tb[nr-1] ;
    som2 = tb[nr-1] ;
    if (nr > 1)
      xo[nr-2] = (double(nr) - 1.5)*(som1 - (2*nr-3)*(2*nr-2)*som2) ;
    for (int i = nr-3 ; i >= 0 ; i-- ) {
	som1 += cx1[i] * tb[i+1] ;
	som2 += tb[i+1] ;
	xo[i] = cx3[i]*(som1 - cx2[i] * som2) ;
    }	// Fin de la boucle sur r
}


		// ---------------------
		// La routine a appeler
		//----------------------
		
		
void d2sdx2_1d(int nr, double** tb, int base_r)
{

		// Routines de derivation
    static void (*d2sdx2_1d[MAX_BASE])(int, double*, double *) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    d2sdx2_1d[i] = _d2sdx2_1d_pas_prevu ;
	}
		// Les routines existantes
	d2sdx2_1d[R_CHEB >> TRA_R] = _d2sdx2_1d_r_cheb ;
	d2sdx2_1d[R_CHEBU >> TRA_R] = _d2sdx2_1d_r_chebu ;
	d2sdx2_1d[R_CHEBP >> TRA_R] = _d2sdx2_1d_r_chebp ;
	d2sdx2_1d[R_CHEBI >> TRA_R] = _d2sdx2_1d_r_chebi ;
	d2sdx2_1d[R_LEG >> TRA_R] = _d2sdx2_1d_r_leg ;
	d2sdx2_1d[R_LEGP >> TRA_R] = _d2sdx2_1d_r_legp ;
	d2sdx2_1d[R_LEGI >> TRA_R] = _d2sdx2_1d_r_legi ;
	d2sdx2_1d[R_JACO02 >> TRA_R] = _d2sdx2_1d_r_jaco02 ;
    }
    
    double *result = new double[nr] ;
    
    d2sdx2_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}
}
