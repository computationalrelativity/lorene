/*
 * Method of the class Base_val to get the values of the theta basis functions
 *  at the theta collocation points.
 *
 * (see file base_val.h for the documentation)
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2001 Jerome Novak
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


char base_val_theta_funct_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2006/05/30 08:30:12  n_vasset
 * Implementation of sine-like bases (T_SIN_P, T_SIN_I, T_COSSIN_SI, etc...).
 *
 * Revision 1.4  2005/05/27 14:54:59  j_novak
 * Added new bases T_COSSIN_CI and T_COS_I
 *
 * Revision 1.3  2004/11/23 15:08:01  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.2  2002/10/16 14:36:30  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2001/11/16  09:28:35  novak
 * The case nt=1 is treated separately
 *
 * Revision 1.2  1999/12/29 10:49:47  eric
 * Mehtode const.
 *
 * Revision 1.1  1999/12/28  12:58:35  eric
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
void theta_funct_pas_prevu(int, double*) ;
void theta_funct_cos_p(int, double*) ;
void theta_funct_cos_i(int, double*) ;
void theta_funct_sin_p(int, double*) ;
void theta_funct_sin_i(int, double*) ;
void theta_funct_cossin_cp(int, double*) ;
void theta_funct_cossin_ci(int, double*) ;
void theta_funct_cossin_sp(int, double*) ;
void theta_funct_cossin_si(int, double*) ;
void theta_funct_cossin_c(int, double*) ;
void theta_funct_cossin_s(int, double*) ;

//************************************************************************
//  user interface : method Base_val::theta_functions
//************************************************************************

const Tbl& Base_val::theta_functions(int l, int nt) const {
    
    const int nmax = 20 ;	    // maximum number of couples (base_t, nt) 
    static int nb_done = 0 ;	    // number of Tbl already computed
    static int base_t_done[nmax] ;  // theta bases already treated
    static int nt_done[nmax] ;	    // number of points already treated
    static Tbl* tab[nmax] ;	    // result for couples (base_t, nt)
    
    static void(*vbasecol[MAX_BASE])(int, double*) ;  // computation routines
    static int dim2[MAX_BASE] ;			 // dim(2) of the Tbl's

    static int premier_appel = 1 ;

    // Initializations at first call
    // -----------------------------
    if (premier_appel == 1) {

	premier_appel = 0 ;

	for (int i=0 ; i<MAX_BASE ; i++) {
	    vbasecol[i] = theta_funct_pas_prevu ;
	    dim2[i] = 0 ; 
	}

	vbasecol[T_COS_P >> TRA_T] = theta_funct_cos_p ;
	vbasecol[T_COS_I >> TRA_T] = theta_funct_cos_i ;
	vbasecol[T_SIN_I >> TRA_T] = theta_funct_sin_i ;
	vbasecol[T_SIN_P >> TRA_T] = theta_funct_sin_p ;
	vbasecol[T_COSSIN_CP >> TRA_T] = theta_funct_cossin_cp ;
	vbasecol[T_COSSIN_CI >> TRA_T] = theta_funct_cossin_ci ;
	vbasecol[T_COSSIN_SP >> TRA_T] = theta_funct_cossin_sp ;
	vbasecol[T_COSSIN_SI >> TRA_T] = theta_funct_cossin_si ;
	vbasecol[T_COSSIN_C >> TRA_T] = theta_funct_cossin_c ;
	vbasecol[T_COSSIN_S >> TRA_T] = theta_funct_cossin_s ;

	dim2[T_COS_P >> TRA_T] = 1 ;
	dim2[T_COS_I >> TRA_T] = 1 ;
	dim2[T_SIN_P >> TRA_T] = 1 ;
	dim2[T_SIN_I >> TRA_T] = 1 ;
	dim2[T_COSSIN_CP >> TRA_T] = 2 ;
	dim2[T_COSSIN_CI >> TRA_T] = 2 ;
	dim2[T_COSSIN_SP >> TRA_T] = 2 ;
	dim2[T_COSSIN_SI >> TRA_T] = 2 ;

    }

    // Computation 
    // -----------

    int base_t = ( b[l] & MSQ_T ) >> TRA_T ;

    // Has this couple (base_t, nt) been previously considered ?
    // ---------------------------------------------------------
    int index = -1 ; 
    for (int i=0; i<nb_done; i++) {
	if ( (base_t_done[i] == base_t) && (nt_done[i] == nt) ) {
	    index = i ; 
	}
    }
    
    // If not, a new computation must be performed 
    // -------------------------------------------
    if (index == -1) {
	if ( nb_done >= nmax ) {
	    cout << "Base_val::theta_functions :  nb_done >= nmax ! " << endl ; 
	    abort() ; 
	}
	
	index = nb_done ; 

	tab[index] = new Tbl( dim2[base_t], nt, nt ) ; 
	(tab[index])->set_etat_qcq() ; 

	vbasecol[base_t](nt, (tab[index])->t ) ; 
	
	base_t_done[index] = base_t ; 
	nt_done[index] = nt ;
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

void theta_funct_pas_prevu(int, double*) {
    
    cout << "Base_val::theta_functions : theta basis not implemented !" 
	 << endl ; 
    abort() ; 
    
}

//==============================================
//  Basis cos(2*j* theta)  T_COS_P
//==============================================

void theta_funct_cos_p(int nt, double* ff) {

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = cos(2*i * theta);	
	}
    }

}

//==============================================
//  Basis cos((2*j+1)* theta)  T_COS_I
//==============================================

void theta_funct_cos_i(int nt, double* ff) {

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = cos((2*i+1) * theta);	
	}
    }

}

//==============================================
//  Basis sin(2*j* theta)  T_SIN_P
//==============================================

void theta_funct_sin_p(int nt, double* ff) {

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = sin(2*i * theta);	
	}
    }

}

//==============================================
//  Basis sin((2*j+1)* theta)  T_SIN_I
//==============================================

void theta_funct_sin_i(int nt, double* ff) {

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = sin((2*i+1) * theta);	
	}
    }

}

//===========================================================
//  Basis cos(2*j* theta)/sin((2*j+1) theta)   T_COSSIN_CP
//===========================================================

void theta_funct_cossin_cp(int nt, double* ff) {
    
    int nt2 = nt*nt ;

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    // cos(2i theta_j)
    // ---------------
    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = cos(2*i * theta);	
	}
    }
    
    // sin((2i+1) theta_j)
    // -------------------

    for (int i = 0; i < nt-1 ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt2+nt*i+ j] = sin((2*i+1) * theta);	
	}
    }
    
    for (int j = 0; j < nt ; j++ ) {
	ff[nt2+nt*(nt-1) + j] = 0 ;	
    }
    
    
}

//===========================================================
//  Basis cos((2*j+1)* theta)/sin(2*j*theta)   T_COSSIN_CI
//===========================================================

void theta_funct_cossin_ci(int nt, double* ff) {
    
    int nt2 = nt*nt ;

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    // cos((2i+1) theta_j)
    // ---------------
    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = cos((2*i+1) * theta);	
	}
    }
    
    // sin(2i theta_j)
    // -------------------

    for (int i = 0; i < nt-1 ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt2+nt*i+ j] = sin(2*i * theta);	
	}
    }
    
    for (int j = 0; j < nt ; j++ ) {
	ff[nt2+nt*(nt-1) + j] = 0 ;	
    }
    
    
}


//===========================================================
//  Basis sin(2*j* theta)/cos((2*j+1) theta)   T_COSSIN_SP
//===========================================================

void theta_funct_cossin_sp(int nt, double* ff) {
    
    int nt2 = nt*nt ;

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    // sin(2i theta_j)
    // ---------------
    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = sin(2*i * theta);	
	}
    }
    
    // cos((2i+1) theta_j)
    // -------------------

    for (int i = 0; i < nt-1 ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt2+nt*i+ j] = cos((2*i+1) * theta);	
	}
    }
    
    for (int j = 0; j < nt ; j++ ) {
	ff[nt2+nt*(nt-1) + j] = 0 ;	
    }
    
    
}

//===========================================================
//  Basis sin((2*j+1)* theta)/cos(2*j*theta)   T_COSSIN_SI
//===========================================================

void theta_funct_cossin_si(int nt, double* ff) {
    
    int nt2 = nt*nt ;

    double xx = ( nt > 1 ? M_PI / double(2*(nt-1)) : 0.) ;

    // sin((2i+1) theta_j)
    // ---------------
    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = sin((2*i+1) * theta);	
	}
    }
    
    // cos(2i theta_j)
    // -------------------

    for (int i = 0; i < nt-1 ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt2+nt*i+ j] = cos(2*i * theta);	
	}
    }
    
    for (int j = 0; j < nt ; j++ ) {
	ff[nt2+nt*(nt-1) + j] = 0 ;	
    }
    
    
}

//===========================================================
//  Basis cos(j* theta)/sin(j*theta)   T_COSSIN_C
//===========================================================

void theta_funct_cossin_c(int nt, double* ff) {
    
    int nt2 = nt*nt ;

    double xx = ( nt > 1 ? M_PI / double(nt-1) : 0.) ;

    // cos(i theta_j)
    // ---------------
    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = cos(i * theta);	
	}
    }
    
    // sin(i theta_j)
    // -------------------

    for (int i = 0; i < nt-1 ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt2+nt*i+ j] = sin(i * theta);	
	}
    }
    
    for (int j = 0; j < nt ; j++ ) {
	ff[nt2+nt*(nt-1) + j] = 0 ;	
    }
    
    
}

//===========================================================
//  Basis sin(j* theta)/cos(j*theta)   T_COSSIN_S
//===========================================================

void theta_funct_cossin_s(int nt, double* ff) {
    
    int nt2 = nt*nt ;

    double xx = ( nt > 1 ? M_PI / double(nt-1) : 0.) ;

    // sin(i theta_j)
    // ---------------
    for (int i = 0; i < nt-1 ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt*i+ j] = sin(i * theta);	
	}
    }

    for (int j = 0; j < nt ; j++ ) {
	ff[nt*(nt-1) + j] = 0 ;	
    }
    
    // cos(i theta_j)
    // -------------------

    for (int i = 0; i < nt ; i++ ) {
	for (int j = 0; j < nt ; j++ ) {
	    double theta = xx*j ;
	    ff[nt2+nt*i+ j] = cos(i * theta);	
	}
    }
           
}

