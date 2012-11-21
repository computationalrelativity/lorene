/*
 *  Methods of the class Star_QI to compute global quantities
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char star_QI_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2012/11/21 14:54:13  c_some
 * First version
 *
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>
#include <cstdlib>

// Lorene headers
#include "compobj.h"
#include "unites.h"
#include "proto_f77.h"

			//------------------------//
			//	Gravitational mass    //
			//------------------------//

double Star_QI::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
		Scalar s_euler = stress_euler.trace(gamma) ; 
	
		//## alternative:
		// assert(*(stress_euler.get_triad()) == mp.get_bvect_spher()) ; 
		// Scalar s_euler = ( stress_euler(1,1) + stress_euler(2,2) ) / a_car 
		// 	+ stress_euler(3,3) / b_car ; 
			
		// Cf. Eq. (4.18) of arXiv:1003.5015v2 with (E+p) U = B * mom_euler(3)
		
	    Scalar source = nn * (ener_euler + s_euler) 
				+ 2 * b_car * mom_euler(3) 
				    * tnphi  ; 
	    source = a_car * bbb * source ;
	    source.std_spectral_base() ; 

	    p_mass_g = new double( source.integrale() ) ;

    }
    
    return *p_mass_g ; 

} 
		

			//------------------------//
			//	Angular momentum      //
			//------------------------//

double Star_QI::angu_mom() const {

    if (p_angu_mom == 0x0) {    // a new computation is required
	
		// Cf. Eq. (4.39) of arXiv:1003.5015v2 with (E+p) U = B * mom_euler(3) 
		
		assert(*(mom_euler.get_triad()) == mp.get_bvect_spher()) ; 
	
		Scalar dens = mom_euler(3) ; 

		dens.mult_r() ;			//  Multiplication by
		dens.set_spectral_va() = (dens.get_spectral_va()).mult_st() ;	//    r sin(theta)

		dens = a_car * b_car * bbb * dens ; 

		p_angu_mom = new double( dens.integrale() ) ;

    }
    
    return *p_angu_mom ; 

}

			//--------------------//
			//	     GRV2	      //
			//--------------------//

double Star_QI::grv2() const {

      using namespace Unites ;	
      if (p_grv2 == 0x0) {    // a new computation is required
	
		assert( *(stress_euler.get_triad()) == mp.get_bvect_spher() ) ; 
		Scalar sou_m =  2 * qpig * a_car * b_car * stress_euler(3,3) ;
	
		Vector dlogn = logn.derive_cov( mp.flat_met_spher() ) ; 
        Scalar sou_q =  1.5 * ak_car
	  - dlogn(1)*dlogn(1) - dlogn(2)*dlogn(2) - dlogn(3)*dlogn(3) ;	
	
		p_grv2 = new double( double(1) - lambda_grv2(sou_m, sou_q) ) ; 
	
      }
    
      return *p_grv2 ; 
      
}


			//--------------------//
			//	     GRV3	      //
			//--------------------//

double Star_QI::grv3(ostream* ost) const {

  using namespace Unites ;	
  
  if (p_grv3 == 0x0) {    // a new computation is required

    Scalar source(mp) ; 
    
    // Gravitational term [cf. Eq. (43) of Gourgoulhon & Bonazzola
    // ------------------	    Class. Quantum Grav. 11, 443 (1994)]
    
    Vector dlogn = logn.derive_cov( mp.flat_met_spher() ) ;

      Scalar alpha = dzeta - logn ; 
      Scalar bet = log( bbb ) ; 
      bet.std_spectral_base() ; 
      
      Vector dalpha = alpha.derive_cov( mp.flat_met_spher() ) ;
      Vector dbet = bet.derive_cov( mp.flat_met_spher() ) ;

      source = 0.75 * ak_car 
	- dlogn(1)*dlogn(1) - dlogn(2)*dlogn(2) - dlogn(3)*dlogn(3) 
	+ 0.5 * ( dalpha(1)*dbet(1) + dalpha(2)*dbet(2) + dalpha(3)*dbet(3) ) ; 
      
      Scalar aa = alpha - 0.5 * bet ; 
      Scalar daadt = aa.srdsdt() ;	// 1/r d/dth
	    
      // What follows is valid only for a mapping of class Map_radial : 
      const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ; 
      if (mpr == 0x0) {
	cout << "Star_rot::grv3: the mapping does not belong"
	     << " to the class Map_radial !" << endl ; 
	abort() ; 
      }
      
      // Computation of 1/tan(theta) * 1/r daa/dtheta
      if (daadt.get_etat() == ETATQCQ) {
	Valeur& vdaadt = daadt.set_spectral_va() ; 
	vdaadt = vdaadt.ssint() ;	// division by sin(theta)
	vdaadt = vdaadt.mult_ct() ;	// multiplication by cos(theta)
      }
      
      Scalar temp = aa.dsdr() + daadt ; 
      temp = ( bbb - a_car/bbb ) * temp ; 
      temp.std_spectral_base() ; 
      
      // Division by r 
      Valeur& vtemp = temp.set_spectral_va() ; 
      vtemp = vtemp.sx() ;    // division by xi in the nucleus
      // Id in the shells
      // division by xi-1 in the ZEC
      vtemp = (mpr->xsr) * vtemp ; // multiplication by xi/r in the nucleus
      //		  by 1/r in the shells
      //		  by r(xi-1) in the ZEC
      
      // In the ZEC, a multiplication by r has been performed instead
      //   of the division: 			
      temp.set_dzpuis( temp.get_dzpuis() + 2 ) ;  
      
      source = bbb * source + 0.5 * temp ; 
          
    source.std_spectral_base() ; 
    
    double int_grav = source.integrale() ; 
    
    // Matter term
    // -----------
    
	Scalar s_euler = stress_euler.trace(gamma) ; 
	
		//## alternative:
		// assert(*(stress_euler.get_triad()) == mp.get_bvect_spher()) ; 
		// Scalar s_euler = ( stress_euler(1,1) + stress_euler(2,2) ) / a_car 
		// 	+ stress_euler(3,3) / b_car ; 
			
    source  = qpig * a_car * bbb * s_euler ;
     
    source.std_spectral_base() ; 

    double int_mat = source.integrale() ; 
    
    // Virial error
    // ------------
    if (ost != 0x0) {
      *ost << "Star_rot::grv3 : gravitational term : " << int_grav 
	   << endl ;
      *ost << "Star_rot::grv3 : matter term :        " << int_mat 
	   << endl ;
    }
    
    p_grv3 = new double( (int_grav + int_mat) / int_mat ) ; 	 
    
  }
  
  return *p_grv3 ; 
  
}

			//----------------------------//
			//     Quadrupole moment      //
			//----------------------------//

double Star_QI::mom_quad() const {

  using namespace Unites ;
    if (p_mom_quad == 0x0) {    // a new computation is required
	
	// Source for of the Poisson equation for nu
	// -----------------------------------------

	Scalar source(mp) ; 
	
	Scalar s_euler = stress_euler.trace(gamma) ; 
	
		//## alternative:
		// assert(*(stress_euler.get_triad()) == mp.get_bvect_spher()) ; 
		// Scalar s_euler = ( stress_euler(1,1) + stress_euler(2,2) ) / a_car 
		// 	+ stress_euler(3,3) / b_car ; 

	    Scalar bet = log(bbb) ; 
	    bet.std_spectral_base() ; 
	    
	    Vector dlogn = logn.derive_cov( mp.flat_met_spher() ) ;
	    Vector dlogn_bet = dlogn + bet.derive_cov( mp.flat_met_spher() ) ;

	    source =  qpig * a_car *( ener_euler + s_euler ) + ak_car 
	      - dlogn(1)*dlogn_bet(1) - dlogn(2)*dlogn_bet(2) - dlogn(3)*dlogn_bet(3)   ; 
	
	source.std_spectral_base() ;

	// Multiplication by -r^2 P_2(cos(theta))
	//  [cf Eq.(7) of Salgado et al. Astron. Astrophys. 291, 155 (1994) ]
	// ------------------------------------------------------------------
	
	// Multiplication by r^2 : 
	// ----------------------
	source.mult_r() ; 
	source.mult_r() ; 
	if (source.check_dzpuis(2)) {
	    source.inc_dzpuis(2) ; 
	}
		
	// Muliplication by cos^2(theta) :
	// -----------------------------
	Scalar temp = source ; 
	
	// What follows is valid only for a mapping of class Map_radial : 
	assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ; 
		
	if (temp.get_etat() == ETATQCQ) {
	    Valeur& vtemp = temp.set_spectral_va() ; 
	    vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
	    vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
	}
	
	// Muliplication by -P_2(cos(theta)) :
	// ----------------------------------
	source = 0.5 * source - 1.5 * temp ; 
	
	// Final result
	// ------------

	p_mom_quad = new double( source.integrale() / qpig ) ; 	 

    }
    
    return *p_mom_quad ; 

}


// Function Star_QI::lambda_grv2

double Star_QI::lambda_grv2(const Scalar& sou_m, const Scalar& sou_q) {

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

	    	case FINJAC:	{
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

	    	case FINJAC:	{
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

