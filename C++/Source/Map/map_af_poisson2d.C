/*
 * Method of the class Map_af for the resolution of the 2-D Poisson
 *  equation.
 *
 * (see file map.h for documentation).
 */

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


char map_af_poisson2d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/09/09 13:54:20  e_gourgoulhon
 *
 * Change of name of the Fortran subroutines
 * 	poisson2d -> poiss2d
 * 	poisson2di -> poiss2di
 * to avoid any clash with Map::poisson2d and Map::poisson2di
 *
 * Revision 1.2  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/10/11  15:15:26  eric
 * 1ere version operationnelle.
 *
 * Revision 2.0  2000/10/09  13:47:10  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"
#include "param.h"
#include "proto_f77.h"

//*****************************************************************************

void Map_af::poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
		       Param& par, Cmp& uu) const {
    
    assert(source_mat.get_etat() != ETATNONDEF) ; 
    assert(source_quad.get_etat() != ETATNONDEF) ; 
    assert(source_mat.get_mp()->get_mg() == mg) ; 
    assert(source_quad.get_mp()->get_mg() == mg) ; 
    assert(uu.get_mp()->get_mg() == mg) ; 

    assert( source_quad.check_dzpuis(4) ) ; 
    
    double& lambda = par.get_double_mod(0) ; 

    // Special case of a vanishing source 
    // ----------------------------------

    if (    (source_mat.get_etat() == ETATZERO) 
	 && (source_quad.get_etat() == ETATZERO) ) {
	
	uu = 0 ; 
	lambda = 1 ; 
	return ; 
    }

    // ---------------------------------------------------------------------
    // Preparation of the parameters for the Fortran subroutines F77_poisson2d
    //  and F77_poisson2di
    // ---------------------------------------------------------------------
    
    int nz = mg->get_nzone() ;
    int np1 = 1 ;		// Axisymmetry enforced
    int nt = mg->get_nt(0) ; 
    int nt2 = 2*nt - 1 ;	// Number of points for the theta sampling
				//  in [0,Pi], instead of [0,Pi/2]
    
    // Array NDL
    // ---------
    int* ndl = new int[nz+4] ; 
    ndl[0] = nz ; 
    for (int l=0; l<nz; l++) {
	ndl[1+l] = mg->get_nr(l) ; 
    }
    ndl[1+nz] = nt2 ; 
    ndl[2+nz] = np1 ; 
    ndl[3+nz] = nz ; 
    
    // Array INDD
    // ----------
    int* indd = new int[nz] ; 
    for (int l=0; l<nz; l++) {
	switch ( mg->get_type_r(l) ) {
	    case RARE : {
		indd[l] = 0 ; 
		break ; 
	    }
	    case FIN : {
		indd[l] = 1 ; 
		break ; 
	    }
	    case UNSURR : {
		indd[l] = 2 ; 
		break ; 
	    }
	    default : {
		cout << "Map_af::poisson2d: unknown type_r !" << endl ; 
		abort() ; 
		break ; 
	    }
	}
    }
    
    // Parameters NDR, NDT, NDP and NDZ
    // --------------------------------
    int nrmax = 0 ;
    for (int l=0; l<nz ; l++) {
	nrmax = ( ndl[1+l] > nrmax ) ? ndl[1+l] : nrmax ; 
    }
    int ndr = nrmax + 5 ;   // Le +5 est impose par les routines de resolution
			    // de l'equation de Poisson du type gr2p3s_
    int ndt = nt2 + 2 ;
    int ndp = np1 + 2 ;
    int ndz = nz ; 
    
    // Array ERRE
    // ----------
    
    double* erre = new double [ndz*ndr] ; 
    
    for (int l=0; l<nz; l++) {
	for (int i=0; i<ndl[1+l]; i++) {
	    double xr = mg->get_grille3d(l)->x[i] ;
	    erre[ ndr*l + i ] = alpha[l] * xr + beta[l]  ;
	}
    }

    // Arrays containing the data
    // --------------------------
    
    int ndrt = ndr*ndt ; 
    int ndrtp = ndr*ndt*ndp ; 
    int taille = ndrtp*ndz ;

    double* tsou_m = new double[ taille ] ;    
    double* tsou_q = new double[ taille ] ;    
    double* tuu = new double[ taille ] ;    
    
    // Initialisation to zero :
    for (int i=0; i<taille; i++) {
	tsou_m[i] = 0 ; 
	tsou_q[i] = 0 ; 
	tuu[i] = 0 ; 
    }
    
    // Copy of source_mat into tsou_m
    // ------------------------------
    const Valeur& va_m = source_mat.va ; 
    assert(va_m.get_etat() == ETATQCQ) ; 
    va_m.coef_i() ; 
    const Mtbl* s_m = va_m.c ; 
    assert(s_m->get_etat() == ETATQCQ) ; 
    
    Base_val base_s = va_m.base ; 

    for (int l=0; l<nz; l++) {
	int nr = mg->get_nr(l) ; 
	int nrt = nr*nt ;  
	if (s_m->t[l]->get_etat() == ETATZERO) {
	    for (int k=0; k<np1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			tsou_m[ndrtp*l + ndrt*k + ndr*j + i] = 0 ;
			// point symetrique par rapport au plan theta = pi/2 :
			tsou_m[ndrtp*l + ndrt*k + ndr*(nt2-1-j) + i] = 0 ;			
		    }
		}
	    }
	    
	}
	else {
	    assert( s_m->t[l]->get_etat() == ETATQCQ ) ; 
	    for (int k=0; k<np1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			double xx = s_m->t[l]->t[nrt*k + nr*j + i] ;
			tsou_m[ndrtp*l + ndrt*k + ndr*j + i] = xx ;
			// point symetrique par rapport au plan theta = pi/2 :
			tsou_m[ndrtp*l + ndrt*k + ndr*(nt2-1-j) + i] = xx ;			
		    }
		}
	    }
	}  // End of case etat != ETATZERO
    }  

    // Copy of source_quad into tsou_q
    // -------------------------------

    if (source_quad.get_etat() != ETATZERO) {

	const Valeur& va_q = source_quad.va ; 
	assert(va_q.get_etat() == ETATQCQ) ; 
	va_q.coef_i() ; 
	const Mtbl* s_q = va_q.c ;  

	assert( va_q.base == base_s ) ; 

	assert(s_q->get_etat() == ETATQCQ) ; 

	for (int l=0; l<nz; l++) {
	    int nr = mg->get_nr(l) ; 
	    int nrt = nr*nt ;  
	    if (s_q->t[l]->get_etat() == ETATZERO) {
		for (int k=0; k<np1; k++) {
		    for (int j=0; j<nt; j++) {
			for (int i=0; i<nr; i++) {
			tsou_q[ndrtp*l + ndrt*k + ndr*j + i] = 0 ;
			// point symetrique par rapport au plan theta = pi/2 :
			tsou_q[ndrtp*l + ndrt*k + ndr*(nt2-1-j) + i] = 0 ;			
			}
		    }
		}
	    
	    }
	    else {
		assert( s_q->t[l]->get_etat() == ETATQCQ ) ; 
		for (int k=0; k<np1; k++) {
		    for (int j=0; j<nt; j++) {
			for (int i=0; i<nr; i++) {
			double xx = s_q->t[l]->t[nrt*k + nr*j + i] ;
			tsou_q[ndrtp*l + ndrt*k + ndr*j + i] = xx ;
			// point symetrique par rapport au plan theta = pi/2 :
			tsou_q[ndrtp*l + ndrt*k + ndr*(nt2-1-j) + i] = xx ;			
			}
		    }
		}
	    }  // End of case s_q->t[l]->etat != ETATZERO
	}
	
    } // End of case source_quad.etat != ETATZERO  

    //-----------------------------------------------------------
    //  Call of the Fortran subroutine poisson2d_ or poisson2di_
    //-----------------------------------------------------------

    int base_t = (va_m.base).get_base_t(0) ; 

    Base_val base_uu(nz) ;  // Output spectral bases

    switch (base_t) {

	case T_COS_P : {
	
	    double lambda0 ; 
	    
	    F77_poiss2d(ndl, &ndr, &ndt, &ndp, indd, erre, tsou_m, tsou_q, 
		       &lambda0, tuu) ;
		      
	    base_uu = base_s ;		// output bases = input bases
		       
	    lambda = lambda0 ; 
	    break ; 	    
	}
	
	case T_SIN_I : {

	    double* tsou = new double[taille] ; 
	    for (int i=0; i<taille; i++) {
		tsou[i] = tsou_m[i] + tsou_q[i] ; 
	    }

	    F77_poiss2di(ndl, &ndr, &ndt, &ndp, indd, erre, tsou, tuu) ;
		
	    base_uu = base_s ;		// output bases = input bases

	    lambda = double(1) ; 
	    
	    delete [] tsou ; 
	    break ; 	    
	}
	
	default : {
	    cout << "Map_af::poisson2d : unkown theta basis !" << endl ;
	    cout << "  basis : " << hex << base_t << endl ; 
	    abort() ; 
	    break ; 
	}
    }
    
    //-------------------------------
    // Copy of tuu into uu
    //-------------------------------

    uu.allocate_all() ; 
    (uu.va).set_etat_c_qcq() ;	// To make sure that the coefficients are
				//  deleted
    
    for (int l=0; l<nz; l++) {
	int nr = mg->get_nr(l) ; 
	for (int k=0; k<mg->get_np(l); k++) {
	    for (int j=0; j<nt; j++) {
		for (int i=0; i<nr; i++) {
		    uu.set(l, k, j, i) = tuu[ndrtp*l + ndr*j + i]  ;
		}
	    }
	}
    }

    (uu.va).set_base( base_uu ) ;	// Bases for spectral expansions
    
    uu.set_dzpuis(0) ; 
   
    // Cleaning
    // --------
    
    delete [] ndl ; 
    delete [] indd ; 
    delete [] erre ; 
    delete [] tsou_m ; 
    delete [] tsou_q ; 
    delete [] tuu ; 
    
    
}


