/*
 * Method of the class Base_val to get the values of the phi basis functions
 *  at the phi collocation points.
 *
 * (see file base_val.h for the documentation)
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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

char base_val_phi_funct_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2013/04/25 15:46:05  j_novak
 * Added special treatment in the case np = 1, for type_p = NONSYM.
 *
 * Revision 1.4  2012/01/17 14:44:35  j_penner
 * Modified phi variables to only use 16 integers in arrays
 *
 * Revision 1.3  2006/05/30 13:06:12  n_vasset
 *   Implemented function P_COSSIN_I in base_val_phi_funct.C
 *
 * Revision 1.2  2002/10/16 14:36:30  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  1999/12/29  10:49:35  eric
 * Methode const.
 *
 * Revision 1.1  1999/12/28  12:58:29  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <math.h>


// Headers Lorene
#include "base_val.h"
#include "type_parite.h"
#include "tbl.h"

// Local prototypes
void phi_funct_pas_prevu(int, double*) ;
void phi_funct_cossin(int, double*) ;
void phi_funct_cossin_p(int, double*) ;
void phi_funct_cossin_i(int, double*) ;

//************************************************************************
//  user interface : method Base_val::phi_functions
//************************************************************************

const Tbl& Base_val::phi_functions(int l, int np) const {
    
    const int nmax = 20 ;	    // maximum number of couples (base_p, np) 
    static int nb_done = 0 ;	    // number of Tbl already computed
    static int base_p_done[nmax] ;  // phi bases already treated
    static int np_done[nmax] ;	    // number of points already treated
    static Tbl* tab[nmax] ;	    // result for couples (base_p, np)
    
    static void(*vbasecol[MAX_BASE_2])(int, double*) ;  // computation routines

    static int premier_appel = 1 ;

    // Initializations at first call
    // -----------------------------
    if (premier_appel == 1) {

	premier_appel = 0 ;

	for (int i=0 ; i<MAX_BASE_2 ; i++) {
	    vbasecol[i] = phi_funct_pas_prevu ;
	}

	vbasecol[P_COSSIN >> TRA_P] = phi_funct_cossin ;
	vbasecol[P_COSSIN_P >> TRA_P] = phi_funct_cossin_p ;
        vbasecol[P_COSSIN_I >> TRA_P] = phi_funct_cossin_p ;

    }

    // Computation 
    // -----------

    int base_p = ( b[l] & MSQ_P ) >> TRA_P ;

    // Has this couple (base_p, np) been previously considered ?
    // ---------------------------------------------------------
    int index = -1 ; 
    for (int i=0; i<nb_done; i++) {
	if ( (base_p_done[i] == base_p) && (np_done[i] == np) ) {
	    index = i ; 
	}
    }
    
    // If not, a new computation must be performed 
    // -------------------------------------------
    if (index == -1) {
	if ( nb_done >= nmax ) {
	    cout << "Base_val::phi_functions :  nb_done >= nmax ! " << endl ; 
	    abort() ; 
	}
	
	index = nb_done ; 

	tab[index] = new Tbl( np+1, np ) ; 
	(tab[index])->set_etat_qcq() ; 

	vbasecol[base_p](np, (tab[index])->t ) ; 
	
	base_p_done[index] = base_p ; 
	np_done[index] = np ;
	nb_done++ ;
	
    }  // end of the case where the computation had to be done


    return *(tab[index]) ;
    
}


//************************************************************************
//  computational subroutines
//************************************************************************

//====================================
//  Unknown case
//====================================

void phi_funct_pas_prevu(int, double*) {
    
    cout << "Base_val::phi_functions : phi basis not implemented !" 
	 << endl ; 
    abort() ; 
    
}

//==============================================
//  Basis P_COSSIN
//==============================================

void phi_funct_cossin(int np, double* ff) {

    double xx = 2.*M_PI / double(np) ;

    if (np == 1) {
	ff[0] = 1. ; // cos (0 * phi)
	ff[1] = 0. ; // sin (0 * phi)
    }
    else {
      for (int i = 0; i < np-1 ; i+=2 ) {
	int m = i/2 ;
	for (int k = 0; k < np ; k++ ) {
	  double phi = xx*k ;
	  ff[np*i + k] = cos(m * phi) ;	
	  ff[np*(i+1) + k] = sin(m * phi) ;	
	}
      }

      for (int k = 0; k < np ; k++ ) {
	double phi = xx*k ;
	ff[np*np + k] = cos(np/2 * phi) ;	
      }
    }

}

//==============================================
//  Basis P_COSSIN_P
//==============================================

void phi_funct_cossin_p(int np, double* ff) {
    
    double xx = M_PI/double(np) ;
    
    for (int i = 0; i < np+1 ; i+=2 ) {
	for (int k = 0; k < np ; k++ ) {
	    double phi = xx*k ;
	    ff[np*i+ k] = cos(i * phi);	
	}
    }

    for (int i = 1; i < np ; i+=2 ) {
	for (int k = 0; k < np ; k++ ) {
	    double phi = xx*k ;
	    ff[np*i+ k] = sin((i-1) * phi);	
	}
    }

    
}

//==============================================
//  Basis P_COSSIN_I
//==============================================

void phi_funct_cossin_i(int np, double* ff) {
    
    double xx = M_PI/double(np) ;
    
    for (int i = 0; i < np+1 ; i+=2 ) {
	for (int k = 0; k < np ; k++ ) {
	    double phi = xx*k ;
	    ff[np*i+ k] = sin(i * phi);	
	}
    }

    for (int i = 1; i < np ; i+=2 ) {
	for (int k = 0; k < np ; k++ ) {
	    double phi = xx*k ;
	    ff[np*i+ k] = cos((i-1) * phi);	
	}
    }

    
}

