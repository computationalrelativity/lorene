/*
 * Functions Star_rot_Dirac_diff::fait_omega_field
 *	     Star_rot_Dirac_diff::fait_prim_field
 *
 *    (see file star_rot_dirac_diff.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Motoyuki Saijo
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


char strot_dirac_diff_faitomeg_C[] = "$Header$" ;

/*
 * $Header$
 *
 */


// Headers Lorene
#include "star_rot_dirac_diff.h"
#include "utilitaires.h"
#include "param.h"

namespace Lorene {
double strot_dirac_diff_fzero(double omeg, const Param& par) ;


		    //-----------------------------------//
		    //		fait_omega_field         //
		    //-----------------------------------//

void Star_rot_Dirac_diff::fait_omega_field(double omeg_min, double omeg_max,
				   double precis, int nitermax) {

    const Mg3d& mg = *(mp.get_mg()) ; 
    int nz = mg.get_nzone() ; 
    int nzm1 = nz - 1 ; 

    Scalar brst2 = gamma.cov()(3,3) ;
    brst2.annule(nzm1, nzm1) ; 
            brst2.mult_rsint() ;
            brst2.mult_rsint() ;
    Scalar nnn2 = qqq * qqq / psi4 ; 
    Scalar betp = beta(3) ;
           betp.div_rsint() ;

    Param par ;
    par.add_scalar(nnn2, 0) ;
    par.add_scalar(brst2, 1) ;
    par.add_scalar(betp, 2) ;

    int l, k, j, i ;
    par.add_int(l, 0) ;
    par.add_int(k, 1) ;
    par.add_int(j, 2) ;
    par.add_int(i, 3) ;

    par.add_star(*this) ;

    // Loop on the grid points
    // -----------------------

    int niter ;
    
    bool prev_zero = (omega_field.get_etat() == ETATZERO) ;
    
    omega_field.allocate_all() ;

    // cout << "Star_rot_Dirac_diff::fait_omega_field: niter : " << '\n' ; 
    for (l=0; l<nzet+1; l++) {
    	Tbl& tom = omega_field.set_domain(l) ;
	for (k=0; k<mg.get_np(l); k++) {
	    for (j=0; j<mg.get_nt(l); j++) {
		for (i=0; i<mg.get_nr(l); i++) {

		    double& omeg =  tom.set(k, j, i) ; 

		    double omeg_min1, omeg_max1 ; 
		    if ( prev_zero || omeg == double(0)) {
			omeg_min1 = omeg_min ; 
			omeg_max1 = omeg_max ; 
		    } 
		    else{
			omeg_min1 = 0.8 * omeg ; 
			omeg_max1 = 1.2 * omeg ; 
		    } 
		
		    omeg = zerosec(strot_dirac_diff_fzero, par, omeg_min1,
			 	   omeg_max1, precis, nitermax, niter) ;
			
		    // cout << "  " << niter ; 

		}
	    }
	}
	// cout << '\n' ; 
    }

    // omega_field is set to 0 far from the star:
    // ---------------------------------------------------
    for (l=nzet+1; l<nz; l++) {
	omega_field.set_domain(l) = 0 ;
    }

    omega_field.std_spectral_base() ;
    
    // Min and Max of Omega
    // --------------------
    
    omega_min = min( omega_field.domain(0) ) ; 
    for (l=1; l<nzet; l++) {
	double xx = min( omega_field.domain(l) ) ; 
	omega_min = (xx < omega_min) ? xx : omega_min ; 
    }
    
    omega_max = max( max( omega_field ) ) ;     
    
    // Update of prim_field
    // --------------------
    fait_prim_field() ; 
    
}

		    //----------------------------------------//
		    //  	strot_dirac_diff_fzero        //
		    //----------------------------------------//

double strot_dirac_diff_fzero(double omeg, const Param& par) {

	const Scalar& nnn2 =  par.get_scalar(0) ;
	const Scalar& brst2 = par.get_scalar(1) ;
    	const Scalar& betp =  par.get_scalar(2) ;
    	int l = par.get_int(0) ;
    	int k = par.get_int(1) ;
   	int j = par.get_int(2) ;
   	int i = par.get_int(3) ;
   	
   	const Star_rot_Dirac_diff* star =
	  dynamic_cast<const Star_rot_Dirac_diff*>(&par.get_star()) ;
   		
     	double fom = star->funct_omega(omeg) ;
     	double omnp = omeg + betp.val_grid_point(l, k, j, i) ;
    	
     	return fom - brst2.val_grid_point(l, k, j, i) * omnp /
     	       ( nnn2.val_grid_point(l, k, j, i) - 
                 brst2.val_grid_point(l, k, j, i) * omnp * omnp ) ;

}


		    //-----------------------------------//
		    //		fait_prim_field          //
		    //-----------------------------------//


void Star_rot_Dirac_diff::fait_prim_field() {
   
    const Mg3d& mg = *(mp.get_mg()) ; 
    int nz = mg.get_nzone() ; 

    // Loop on the grid points in the vicinity of the star
    // ---------------------------------------------------

    prim_field.allocate_all() ;

    for (int l=0; l<nzet+1; l++) {
    	Tbl& tprim = prim_field.set_domain(l) ;
	for (int k=0; k<mg.get_np(l); k++) {
	    for (int j=0; j<mg.get_nt(l); j++) {
		for (int i=0; i<mg.get_nr(l); i++) {
		
			tprim.set(k, j, i) = 
			    primfrot( omega_field.val_grid_point(l, k, j, i), 
                                      par_frot ) ; 
			
		}
	    }
	}
    }

    // prim_field is set to 0 far from the star:
    // -----------------------------------------
    for (int l=nzet+1; l<nz; l++) {
	prim_field.set_domain(l) = 0 ;
    }

    prim_field.std_spectral_base() ;
    
    
}
}
