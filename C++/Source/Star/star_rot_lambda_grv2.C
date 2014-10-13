/*
 * Method Star_rot::lambda_grv2.
 *
 * (see file star_rot.h for documentation)
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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


char star_rot_lambda_grv2_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:39  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2013/06/05 15:10:43  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "star_rot.h"
#include "proto_f77.h"

namespace Lorene {
double Star_rot::lambda_grv2(const Scalar& sou_m, const Scalar& sou_q) {

	const Map_radial* mprad = dynamic_cast<const Map_radial*>( &sou_m.get_mp() ) ;
	
	if (mprad == 0x0) {
		cout << "Star_rot::lambda_grv2: the mapping of sou_m does not"
			 << endl << " belong to the class Map_radial !" << endl ;
		abort() ;
	} 	

	assert( &sou_q.get_mp() == mprad ) ;
	
	sou_q.check_dzpuis(4) ;
	
	const Mg3d* mg = mprad->get_mg() ;
	int nz = mg->get_nzone() ;
		
	// Construction of a Map_af which coincides with *mp on the equator
    // ----------------------------------------------------------------

    double theta0 = M_PI / 2 ;	    // Equator
    double phi0 = 0 ;

    Map_af mpaff(*mprad) ;

    for (int l=0 ; l<nz ; l++) {
		double rmax = mprad->val_r(l, double(1), theta0, phi0) ;
		switch ( mg->get_type_r(l) ) {
	    	case RARE:	{
				double rmin = mprad->val_r(l, double(0), theta0, phi0) ;
				mpaff.set_alpha(rmax - rmin, l) ;
				mpaff.set_beta(rmin, l) ;
				break ;
	    	}
	
	    	case FIN:	{
				double rmin = mprad->val_r(l, double(-1), theta0, phi0) ;
				mpaff.set_alpha( double(.5) * (rmax - rmin), l ) ;
				mpaff.set_beta( double(.5) * (rmax + rmin), l) ;
				break ;
	    	}

	
	    	case UNSURR: {
				double rmin = mprad->val_r(l, double(-1), theta0, phi0) ;
				double umax = double(1) / rmin ;
				double umin = double(1) / rmax ;
				mpaff.set_alpha( double(.5) * (umin - umax),  l) ;
				mpaff.set_beta( double(.5) * (umin + umax), l) ;
				break ;
	    	}
	
	    	default: {
				cout << "Star_rot::lambda_grv2: unknown type_r ! " << endl ;
				abort () ;
				break ;
	    	}
	
		}
    }


	// Reduced Jacobian of
	// the transformation  (r,theta,phi) <-> (dzeta,theta',phi')
	// ------------------------------------------------------------
	
	Mtbl jac = 1 / ( (mprad->xsr) * (mprad->dxdr) ) ;	
								// R/x dR/dx in the nucleus
								// R dR/dx   in the shells
								// - U/(x-1) dU/dx in the ZEC						
	for (int l=0; l<nz; l++) {
		switch ( mg->get_type_r(l) ) {
	    	case RARE:	{
	    		double a1 = mpaff.get_alpha()[l] ;
				*(jac.t[l]) =  *(jac.t[l]) / (a1*a1) ;
				break ;
	    	}
	
	    	case FIN:	{
				double a1 = mpaff.get_alpha()[l] ;
				double b1 = mpaff.get_beta()[l] ;
				assert( jac.t[l]->get_etat() == ETATQCQ ) ;
				double* tjac = jac.t[l]->t ;
				double* const xi = mg->get_grille3d(l)->x ;
				for (int k=0; k<mg->get_np(l); k++) {
					for (int j=0; j<mg->get_nt(l); j++) {
						for (int i=0; i<mg->get_nr(l); i++) {
							*tjac = *tjac /
									(a1 * (a1 * xi[i] + b1) ) ;
							tjac++ ; 	
						}
					}
				}				
				
				break ;
	    	}

	
	    	case UNSURR: {
	    		double a1 = mpaff.get_alpha()[l] ;
				*(jac.t[l]) = - *(jac.t[l]) / (a1*a1) ;
				break ;
	    	}
	
	    	default: {
				cout << "Star_rot::lambda_grv2: unknown type_r ! " << endl ;
				abort () ;
				break ;
	    	}
	
		}
	
	}


	// Multiplication of the sources by the reduced Jacobian:
	// -----------------------------------------------------
		
	Mtbl s_m(mg) ;
	if ( sou_m.get_etat() == ETATZERO ) {
		s_m = 0 ;
	}
	else{
		assert(sou_m.get_spectral_va().get_etat() == ETATQCQ) ;	
		sou_m.get_spectral_va().coef_i() ;	
		s_m = *(sou_m.get_spectral_va().c) ;
    }
		
	Mtbl s_q(mg) ;
	if ( sou_q.get_etat() == ETATZERO ) {
		s_q = 0 ;
	}
	else{
		assert(sou_q.get_spectral_va().get_etat() == ETATQCQ) ;	
		sou_q.get_spectral_va().coef_i() ;	
		s_q = *(sou_q.get_spectral_va().c) ;
    }
			
	s_m *= jac ;
	s_q *= jac ;
		
	
	// Preparations for the call to the Fortran subroutine
	// ---------------------------------------------------								
	
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

	// Parameters NDR, NDT, NDP
    // ------------------------
    int nrmax = 0 ;
    for (int l=0; l<nz ; l++) {
		nrmax = ( ndl[1+l] > nrmax ) ? ndl[1+l] : nrmax ;
    }
    int ndr = nrmax + 5 ;
    int ndt = nt2 + 2 ;
    int ndp = np1 + 2 ;

    // Array ERRE
    // ----------

    double* erre = new double [nz*ndr] ;

    for (int l=0; l<nz; l++) {
		double a1 = mpaff.get_alpha()[l] ;
		double b1 = mpaff.get_beta()[l] ;
		for (int i=0; i<ndl[1+l]; i++) {
	    	double xi = mg->get_grille3d(l)->x[i] ;
	    	erre[ ndr*l + i ] = a1 * xi + b1 ;
		}
    }

    // Arrays containing the data
    // --------------------------

    int ndrt = ndr*ndt ;
    int ndrtp = ndr*ndt*ndp ;
    int taille = ndrtp*nz ;

    double* tsou_m = new double[ taille ] ;
    double* tsou_q = new double[ taille ] ;

    // Initialisation to zero :
    for (int i=0; i<taille; i++) {
		tsou_m[i] = 0 ;
		tsou_q[i] = 0 ;
    }

    // Copy of s_m into tsou_m
    // -----------------------

    for (int l=0; l<nz; l++) {
	   for (int k=0; k<np1; k++) {
			for (int j=0; j<nt; j++) {
		 		for (int i=0; i<mg->get_nr(l); i++) {
					double xx = s_m(l, k, j, i) ;
					tsou_m[ndrtp*l + ndrt*k + ndr*j + i] = xx ;
					// point symetrique par rapport au plan theta = pi/2 :
					tsou_m[ndrtp*l + ndrt*k + ndr*(nt2-1-j) + i] = xx ;			
		   		}
			}
	  	}
    }

    // Copy of s_q into tsou_q
    // -----------------------

    for (int l=0; l<nz; l++) {
	   for (int k=0; k<np1; k++) {
			for (int j=0; j<nt; j++) {
		 		for (int i=0; i<mg->get_nr(l); i++) {
					double xx = s_q(l, k, j, i) ;
					tsou_q[ndrtp*l + ndrt*k + ndr*j + i] = xx ;
					// point symetrique par rapport au plan theta = pi/2 :
					tsou_q[ndrtp*l + ndrt*k + ndr*(nt2-1-j) + i] = xx ;			
		   		}
			}
	  	}
    }

	
    // Computation of the integrals
    // ----------------------------

    double int_m, int_q ;
    F77_integrale2d(ndl, &ndr, &ndt, &ndp, erre, tsou_m, &int_m) ;
    F77_integrale2d(ndl, &ndr, &ndt, &ndp, erre, tsou_q, &int_q) ;

    // Cleaning
    // --------

    delete [] ndl ;
    delete [] erre ;
    delete [] tsou_m ;
    delete [] tsou_q ;

    // Computation of lambda
    // ---------------------

    double lambda ;
    if ( int_q != double(0) ) {
		lambda = - int_m / int_q ;
    }
    else{
		lambda = 0 ;
    }
	
    return lambda ;
	
}
}
