/*
 * Functions Et_rot_diff::fait_omega_field
 *	     Et_rot_diff::fait_prim_field
 *
 *    (see file et_rot_diff.h for documentation).
 *
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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


char et_rot_diff_faitomeg_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/05/14 20:07:00  e_gourgoulhon
 * Suppressed the outputs (cout)
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2001/10/25  09:43:18  eric
 * omega_min est determine dans l'etoile seulement (l<nzet).
 *
 * Revision 1.2  2001/10/25  09:21:29  eric
 * Recherche de Omega dans un intervalle de 20% autour de la valeur precedente.
 *
 * Revision 1.1  2001/10/19  08:18:23  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Headers Lorene
#include "et_rot_diff.h"
#include "utilitaires.h"
#include "param.h"

double et_rot_diff_fzero(double omeg, const Param& par) ;


		    //-----------------------------------//
		    //		fait_omega_field         //
		    //-----------------------------------//

void Et_rot_diff::fait_omega_field(double omeg_min, double omeg_max,
				   double precis, int nitermax) {
    
    const Mg3d& mg = *(mp.get_mg()) ; 
    int nz = mg.get_nzone() ; 
    int nzm1 = nz - 1 ; 

    // Computation of B^2 r^2 sin^2(theta):
    // -----------------------------------
    
    Cmp brst2 = bbb() ;
    brst2.annule(nzm1) ; 
    brst2.std_base_scal() ;
    brst2.mult_rsint() ;	    //  Multiplication by r sin(theta)
    brst2 = brst2*brst2 ;
    
    Cmp nnn2 = nnn() * nnn() ; 

    Param par ; 
    par.add_cmp(nnn2, 0) ; 
    par.add_cmp(brst2, 1) ; 
    par.add_cmp(nphi(), 2) ;

    int l, k, j, i ;
    par.add_int(l, 0) ;
    par.add_int(k, 1) ;
    par.add_int(j, 2) ;
    par.add_int(i, 3) ;

    par.add_etoile(*this) ;

    // Loop on the grid points
    // -----------------------

    int niter ;
    
    bool prev_zero = (omega_field.get_etat() == ETATZERO) ;
    
    omega_field.allocate_all() ;

    // cout << "Et_rot_diff::fait_omega_field: niter : " << endl ; 
    for (l=0; l<nzet+1; l++) {
    	Tbl& tom = omega_field.set().set(l) ;
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
		
		    omeg = zerosec(et_rot_diff_fzero, par, omeg_min1,
			 	   omeg_max1, precis, nitermax, niter) ;
			
		    // cout << "  " << niter ; 

		}
	    }
	}
	// cout << endl ; 
    }

    // omega_field is set to 0 far from the star:
    // ---------------------------------------------------
    for (l=nzet+1; l<nz; l++) {
	omega_field.set().set(l) = 0 ;
    }

    omega_field.set_std_base() ;
    
    // Min and Max of Omega
    // --------------------
    
    omega_min = min( omega_field()(0) ) ; 
    for (l=1; l<nzet; l++) {
	double xx = min( omega_field()(l) ) ; 
	omega_min = (xx < omega_min) ? xx : omega_min ; 
    }
    
    omega_max = max( max( omega_field() ) ) ;     
    
    // Update of prim_field
    // --------------------
    fait_prim_field() ; 
    
}

		    //-----------------------------------//
		    //		et_rot_diff_fzero        //
		    //-----------------------------------//

double et_rot_diff_fzero(double omeg, const Param& par) {

	const Cmp& nnn2 =  par.get_cmp(0) ;
	const Cmp& brst2 = par.get_cmp(1) ;
    	const Cmp& nphi =  par.get_cmp(2) ;
    	int l = par.get_int(0) ;
    	int k = par.get_int(1) ;
   	int j = par.get_int(2) ;
   	int i = par.get_int(3) ;
   	
   	const Et_rot_diff* et =
   		dynamic_cast<const Et_rot_diff*>(&par.get_etoile()) ;
   		
     	double fom = et->funct_omega(omeg) ;
     	double omnp = omeg - nphi(l, k, j, i) ;
    	
     	return fom - brst2(l, k, j, i) * omnp /
     	       ( nnn2(l, k, j, i) - brst2(l, k, j, i) * omnp*omnp ) ;

}


		    //-----------------------------------//
		    //		fait_prim_field          //
		    //-----------------------------------//


void Et_rot_diff::fait_prim_field() {
   
    const Mg3d& mg = *(mp.get_mg()) ; 
    int nz = mg.get_nzone() ; 

    // Loop on the grid points in the vicinity of the star
    // ---------------------------------------------------

    prim_field.allocate_all() ;

    for (int l=0; l<nzet+1; l++) {
    	Tbl& tprim = prim_field.set().set(l) ;
	for (int k=0; k<mg.get_np(l); k++) {
	    for (int j=0; j<mg.get_nt(l); j++) {
		for (int i=0; i<mg.get_nr(l); i++) {
		
			tprim.set(k, j, i) = 
			    primfrot( omega_field()(l, k, j, i), par_frot ) ; 
			
		}
	    }
	}
    }

    // prim_field is set to 0 far from the star:
    // -----------------------------------------
    for (int l=nzet+1; l<nz; l++) {
	prim_field.set().set(l) = 0 ;
    }

    prim_field.set_std_base() ;
    
    
}
