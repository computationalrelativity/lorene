/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose-Luis Jaramillo
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


char binhor_equations_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.18  2007/08/22 16:12:33  f_limousin
 * Changed the name of the function computing \tilde{\gamma}_{ij}
 *
 * Revision 1.17  2007/04/13 15:28:55  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 * Revision 1.16  2006/08/01 14:37:19  f_limousin
 * New version
 *
 * Revision 1.15  2006/06/29 08:51:00  f_limousin
 * *** empty log message ***
 *
 * Revision 1.14  2006/05/24 16:56:37  f_limousin
 * Many small modifs.
 *
 * Revision 1.13  2005/11/15 14:04:00  f_limousin
 * Minor change to control the resolution of the equation for psi.
 *
 * Revision 1.12  2005/10/23 16:39:43  f_limousin
 * Simplification of the equation in the case of a conformally
 * flat metric and maximal slicing
 *
 * Revision 1.11  2005/09/13 18:33:15  f_limousin
 * New function vv_bound_cart_bin(double) for computing binaries with
 * berlin condition for the shift vector.
 * Suppress all the symy and asymy in the importations.
 *
 * Revision 1.10  2005/07/11 08:21:57  f_limousin
 * Implementation of a new boundary condition for the lapse in the binary
 * case : boundary_nn_Dir_lapl().
 *
 * Revision 1.9  2005/06/09 16:12:04  f_limousin
 * Implementation of amg_mom_adm().
 *
 * Revision 1.8  2005/04/29 14:02:44  f_limousin
 * Important changes : manage the dependances between quantities (for
 * instance psi and psi4). New function write_global(ost).
 *
 * Revision 1.7  2005/03/10 17:21:52  f_limousin
 * Add the Berlin boundary condition for the shift.
 * Some changes to avoid warnings.
 *
 * Revision 1.6  2005/03/03 10:29:02  f_limousin
 * Delete omega in the parameters of the function boundary_beta_z().
 *
 * Revision 1.5  2005/02/24 17:25:23  f_limousin
 * The boundary conditions for psi, N and beta are now parameters in
 * par_init.d and par_coal.d.
 *
 * Revision 1.4  2005/02/11 18:21:38  f_limousin
 * Dirichlet_binaire and neumann_binaire are taking Scalars as arguments
 * instead of Cmps.
 *
 * Revision 1.3  2005/02/07 10:46:28  f_limousin
 * Many changes !! The sources are written differently to minimize the
 * numerical error, the boundary conditions are changed...
 *
 * Revision 1.2  2004/12/31 15:41:26  f_limousin
 * Correction of an error
 *
 * Revision 1.1  2004/12/29 16:11:34  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tensor.h"
#include "tenseur.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
//#include "graphique.h"

// Resolution for the lapse
// ------------------------

void Bin_hor::solve_lapse (double precision, double relax, int bound_nn,
			   double lim_nn) {
    
    assert ((relax >0) && (relax<=1)) ;
    
    cout << "-----------------------------------------------" << endl ;
    cout << "Resolution LAPSE" << endl ;
    
    Scalar lapse_un_old (hole1.n_auto) ;
    Scalar lapse_deux_old (hole2.n_auto) ;

    Sym_tensor taa_un = hole1.aa.up_down(hole1.tgam) ;         
    Scalar aa_quad_un = contract(taa_un, 0, 1, hole1.aa_auto, 0, 1) ; 

    Sym_tensor taa_deux = hole2.aa.up_down(hole2.tgam) ;         
    Scalar aa_quad_deux = contract(taa_deux, 0, 1, hole2.aa_auto, 0, 1) ; 
       
    Tensor hdirac1 (contract((hole1.hh).derive_cov(hole1.ff),0,2)) ;
    Tensor hdirac2 (contract((hole2.hh).derive_cov(hole2.ff),0,2)) ;

    // Source 1
    // --------
 
    Scalar source_un (hole1.mp) ;
    
    // Conformally flat
    /*
    source_un = hole1.get_psi4()*hole1.nn*aa_quad_un
      -2.*contract(hole1.dn, 0, hole1.psi_auto
		    .derive_con(hole1.ff), 0)/hole1.psi ;
    */
    
    source_un = hole1.get_psi4()*( hole1.nn*( aa_quad_un + 0.3333333333333333 * 
					    hole1.trK*hole1.trK*hole1.decouple)
			       - hole1.trK_point*hole1.decouple )     

       -2.*contract(contract(hole1.hh, 0, hole1.n_auto
	       .derive_cov(hole1.ff), 0), 0, hole1.dpsi, 0)/hole1.psi

	-2.*contract(hole1.dn, 0, hole1.psi_auto
		    .derive_con(hole1.ff), 0)/hole1.psi ;

    - contract(hdirac1, 0, hole1.n_auto.derive_cov(hole1.ff), 0) ; 


    Scalar tmp_un (hole1.mp) ;
    
    tmp_un = hole1.get_psi4()* contract(hole1.beta_auto, 0, hole1.trK.
				    derive_cov(hole1.ff), 0) 
	- contract( hole1.hh, 0, 1, hole1.n_auto.derive_cov(hole1.ff)
		    .derive_cov(hole1.ff), 0, 1 ) ;
        
    tmp_un.inc_dzpuis() ; // dzpuis: 3 -> 4

    source_un += tmp_un ;


    // Source 2
    // ---------

    Scalar source_deux (hole2.mp) ;

    // Conformally flat
    /*
    source_deux = hole2.get_psi4()*hole2.nn*aa_quad_deux
      -2.*contract(hole2.dn, 0, hole2.psi_auto
		    .derive_con(hole2.ff), 0)/hole2.psi ;
    */
    
    source_deux = hole2.get_psi4()*( hole2.nn*( aa_quad_deux + 0.3333333333333333
					  * hole2.trK*hole2.trK*hole2.decouple)
			       - hole2.trK_point*hole2.decouple ) 

	-2.*contract(contract(hole2.hh, 0, hole2.n_auto
	       .derive_cov(hole2.ff), 0), 0, hole2.dpsi, 0)/hole2.psi

	-2.*contract(hole2.dn, 0, hole2.psi_auto
		     .derive_con(hole2.ff), 0)/hole2.psi ;

     - contract(hdirac2, 0, hole2.n_auto.derive_cov(hole2.ff), 0) ; 

    Scalar tmp_deux (hole2.mp) ;
    
    tmp_deux = hole2.get_psi4()* contract(hole2.beta_auto, 0, hole2.trK.
				    derive_cov(hole2.ff), 0) 
	- contract( hole2.hh, 0, 1, hole2.n_auto.derive_cov(hole2.ff)
		    .derive_cov(hole2.ff), 0, 1 ) ;
        
    tmp_deux.inc_dzpuis() ; // dzpuis: 3 -> 4
  	   
    source_deux += tmp_deux ;
    
    cout << "source lapse" << endl << norme(source_un) << endl ;

    // Boundary conditions and resolution
    // -----------------------------------

    Valeur lim_un (hole1.mp.get_mg()-> get_angu()) ;
    Valeur lim_deux (hole1.mp.get_mg()-> get_angu()) ;

    Scalar n_un_temp (hole1.n_auto) ;
    Scalar n_deux_temp (hole2.n_auto) ;

    switch (bound_nn) {

	case 0 : {

	  lim_un = hole1.boundary_nn_Dir(lim_nn) ;
	  lim_deux = hole2.boundary_nn_Dir(lim_nn) ;
	   
	    n_un_temp = n_un_temp - 1./2. ;
	    n_deux_temp = n_deux_temp - 1./2. ;
 
	    dirichlet_binaire (source_un, source_deux, lim_un, lim_deux, 
			       n_un_temp, n_deux_temp, 0, precision) ;
 	    break ;
	}
	    
	case 1 : {

	    lim_un = hole1.boundary_nn_Neu(lim_nn) ;
	    lim_deux = hole2.boundary_nn_Neu(lim_nn) ;

	    neumann_binaire (source_un, source_deux, lim_un, lim_deux, 
			       n_un_temp, n_deux_temp, 0, precision) ;
	    break ;
	}
	    
	default : {
	    cout << "Unexpected type of boundary conditions for the lapse!" 
		 << endl 
		 << "  bound_nn = " << bound_nn << endl ; 
	    abort() ;
	    break ; 
	}
	    
    } // End of switch  

    n_un_temp = n_un_temp + 1./2. ;
    n_deux_temp = n_deux_temp + 1./2. ;
   
    n_un_temp.raccord(3) ;
    n_deux_temp.raccord(3) ;
    
    // Check: has the Poisson equation been correctly solved ?
    // -----------------------------------------------------
    
    int nz = hole1.mp.get_mg()->get_nzone() ;
    cout << "lapse auto" << endl << norme (n_un_temp) << endl ;    
    Tbl tdiff_nn = diffrel(n_un_temp.laplacian(), source_un) ;
    
    cout << 
	"Relative error in the resolution of the equation for the lapse : "
	 << endl ; 
    for (int l=0; l<nz; l++) {
	cout << tdiff_nn(l) << "  " ; 
    }
    cout << endl ;

    // Relaxation :
    // -------------

    n_un_temp = relax*n_un_temp + (1-relax)*lapse_un_old ;
    n_deux_temp = relax*n_deux_temp + (1-relax)*lapse_deux_old ;
 
    hole1.n_auto = n_un_temp ;
    hole2.n_auto = n_deux_temp ;

    hole1.n_comp_import(hole2) ;
    hole2.n_comp_import(hole1) ;

}


// Resolution for Psi
// -------------------

void Bin_hor::solve_psi (double precision, double relax, int bound_psi) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "-----------------------------------------------" << endl ;
    cout << "Resolution PSI" << endl ;
    
    Scalar psi_un_old (hole1.psi_auto) ;
    Scalar psi_deux_old (hole2.psi_auto) ;
    
    Sym_tensor taa_un = hole1.aa.up_down(hole1.tgam) ;         
    Scalar aa_quad_un = contract(taa_un, 0, 1, hole1.aa_auto, 0, 1) ; 

    Sym_tensor taa_deux = hole2.aa.up_down(hole2.tgam) ;         
    Scalar aa_quad_deux = contract(taa_deux, 0, 1, hole2.aa_auto, 0, 1) ; 
    
    Tensor hdirac1 (contract((hole1.hh).derive_cov(hole1.ff),0,2)) ;
    Tensor hdirac2 (contract((hole2.hh).derive_cov(hole2.ff),0,2)) ;

    // Source 1
    // ---------
    
    Scalar source_un (hole1.mp) ;
    /*
    // Conformally flat source
    source_un.annule_hard() ;
    source_un.set_dzpuis(4) ;
    source_un += - hole1.psi*hole1.get_psi4()* 0.125* aa_quad_un ;
    source_un.std_spectral_base() ;
    */
    
    Scalar tmp_un (hole1.mp) ;

    tmp_un = 0.125* hole1.psi_auto * (hole1.tgam).ricci_scal() 
      - contract(hole1.hh, 0, 1, hole1.psi_auto.derive_cov(hole1.ff)
		 .derive_cov(hole1.ff), 0, 1 ) ;
    tmp_un.inc_dzpuis() ; // dzpuis : 3 -> 4
    
    tmp_un -= contract(hdirac1, 0, hole1.psi_auto
		    .derive_cov(hole1.ff), 0) ;  

    source_un = tmp_un - hole1.psi*hole1.get_psi4()* ( 0.125* aa_quad_un 
       	   - 8.33333333333333e-2* hole1.trK*hole1.trK*hole1.decouple ) ;
    source_un.std_spectral_base() ;
    


    // Source 2
    // ---------
   
    Scalar source_deux (hole2.mp) ;
    /*
    // Conformally flat source
    source_deux.annule_hard() ;
    source_deux.set_dzpuis(4) ;
    source_deux += - hole2.psi*hole2.get_psi4()* 0.125* aa_quad_deux ;
    source_deux.std_spectral_base() ;
    */
    

    Scalar tmp_deux (hole2.mp) ;

    tmp_deux = 0.125* hole2.psi_auto * (hole2.tgam).ricci_scal() 
      - contract(hole2.hh, 0, 1, hole2.psi_auto.derive_cov(hole2.ff)
		 .derive_cov(hole2.ff), 0, 1 ) ;
    tmp_deux.inc_dzpuis() ; // dzpuis : 3 -> 4
    
    tmp_deux -= contract(hdirac2, 0, hole2.psi_auto
		    .derive_cov(hole2.ff), 0) ;  

    source_deux = tmp_deux - hole2.psi*hole2.get_psi4()* ( 0.125* aa_quad_deux 
       	   - 8.33333333333333e-2* hole2.trK*hole2.trK*hole2.decouple ) ;
    source_deux.std_spectral_base() ;
    

    cout << "source psi" << endl << norme(source_un) << endl ;

    // Boundary conditions and resolution :
    // ------------------------------------

    Valeur lim_un (hole1.mp.get_mg()-> get_angu()) ;
    Valeur lim_deux (hole1.mp.get_mg()-> get_angu()) ;

    Scalar psi_un_temp (hole1.psi_auto) ;
    Scalar psi_deux_temp (hole2.psi_auto) ;

    switch (bound_psi) {

	case 0 : {

	    lim_un = hole1.boundary_psi_app_hor() ;
	    lim_deux = hole2.boundary_psi_app_hor() ;

	    neumann_binaire (source_un, source_deux, lim_un, lim_deux, 
			     psi_un_temp, psi_deux_temp, 0, precision) ;
	    break ;
	}

	default : {
	    cout << "Unexpected type of boundary conditions for psi!" 
		 << endl 
		 << "  bound_psi = " << bound_psi << endl ; 
	    abort() ;
	    break ; 
	}
	    
    } // End of switch  

    psi_un_temp = psi_un_temp + 1./2. ;
    psi_deux_temp = psi_deux_temp + 1./2. ;
     
    psi_un_temp.raccord(3) ;
    psi_deux_temp.raccord(3) ;
    
    // Check: has the Poisson equation been correctly solved ?
    // -----------------------------------------------------
    
    int nz = hole1.mp.get_mg()->get_nzone() ;
    cout << "psi auto" << endl << norme (psi_un_temp) << endl ;    
    Tbl tdiff_psi = diffrel(psi_un_temp.laplacian(), source_un) ;
    
    cout << 
	"Relative error in the resolution of the equation for psi : "
	 << endl ; 
    for (int l=0; l<nz; l++) {
	cout << tdiff_psi(l) << "  " ; 
    }
    cout << endl ;

    // Relaxation :
    // -------------

    psi_un_temp = relax*psi_un_temp + (1-relax)*psi_un_old ;
    psi_deux_temp = relax*psi_deux_temp + (1-relax)*psi_deux_old ;
    
    hole1.psi_auto = psi_un_temp ;
    hole2.psi_auto = psi_deux_temp ;

    hole1.psi_comp_import(hole2) ;
    hole2.psi_comp_import(hole1) ;

    hole1.set_der_0x0() ;
    hole2.set_der_0x0() ;

    set_hh_Samaya() ;

}


// Resolution for shift with omega fixed.
// --------------------------------------

void Bin_hor::solve_shift (double precision, double relax, int bound_beta,
			   double omega_eff) {
    
    cout << "------------------------------------------------" << endl ;
    cout << "Resolution shift : Omega = " << omega << endl ;
    
    Sym_tensor taa_un = hole1.aa.up_down(hole1.tgam) ;         
    Scalar aa_quad_un = contract(taa_un, 0, 1, hole1.aa_auto, 0, 1) ; 

    Sym_tensor taa_deux = hole2.aa.up_down(hole2.tgam) ;         
    Scalar aa_quad_deux = contract(taa_deux, 0, 1, hole2.aa_auto, 0, 1) ; 
       
    Tensor hdirac1 (contract((hole1.hh).derive_cov(hole1.ff),0,2)) ;
    Tensor hdirac2 (contract((hole2.hh).derive_cov(hole2.ff),0,2)) ;

    // Source 1
    // ---------

    Vector source_un (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    /*
    // Conformally flat source
    source_un = 2.* contract(hole1.aa_auto, 1, hole1.dn, 0)
      - 12.*hole1.nn*contract(hole1.aa_auto, 1, hole1.dpsi, 0)
      /hole1.psi;
    */

    
    Vector tmp_vect_un (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    
    source_un = 2.* contract(hole1.aa_auto, 1, hole1.dn, 0)
      - 12.*hole1.nn*contract(hole1.aa_auto, 1, hole1.dpsi, 0)
      /hole1.psi;

     
    tmp_vect_un = 2./3.* hole1.trK.derive_con(hole1.tgam)
	* hole1.decouple ;
    tmp_vect_un.inc_dzpuis() ;
   
    source_un += 2.* hole1.nn * ( tmp_vect_un
			- contract(hole1.tgam.connect().get_delta(), 1, 2, 
				     hole1.aa_auto, 0, 1) ) ;

    Vector vtmp_un = contract(hole1.hh, 0, 1, 
                           hole1.beta_auto.derive_cov(hole1.ff)
			   .derive_cov(hole1.ff), 1, 2)
      + 1./3.*contract(hole1.hh, 1, hole1.beta_auto
			       .divergence(hole1.ff).derive_cov(hole1.ff), 0) 
      - hdirac1.derive_lie(hole1.beta_auto) 
      + hole1.gamt_point.divergence(hole1.ff)*hole1.decouple ;      
    vtmp_un.inc_dzpuis() ; // dzpuis: 3 -> 4

    source_un -= vtmp_un ; 
        
    source_un += 2./3.* hole1.beta_auto.divergence(hole1.ff) 
	* hdirac1 ;

    source_un.std_spectral_base() ;
    

    // Source 2
    // ---------

    Vector source_deux (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
    /* 
    // Conformally flat source
    source_deux = 2.* contract(hole2.aa_auto, 1, hole2.dn, 0)
      - 12.*hole2.nn*contract(hole2.aa_auto, 1, hole2.dpsi, 0)
      /hole2.psi;
    */

    
    Vector tmp_vect_deux (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
    
    source_deux = 2.* contract(hole2.aa_auto, 1, hole2.dn, 0)
      - 12.*hole2.nn*contract(hole2.aa_auto, 1, hole2.dpsi, 0)
      /hole2.psi;
     
    tmp_vect_deux = 2./3.* hole2.trK.derive_con(hole2.tgam)
	* hole2.decouple ;
    tmp_vect_deux.inc_dzpuis() ;
   
    source_deux += 2.* hole2.nn * ( tmp_vect_deux
		       - contract(hole2.tgam.connect().get_delta(), 1, 2, 
				     hole2.aa*hole2.decouple, 0, 1) ) ;

    Vector vtmp_deux = contract(hole2.hh, 0, 1, 
                           hole2.beta_auto.derive_cov(hole2.ff)
			   .derive_cov(hole2.ff), 1, 2)
      + 1./3.*contract(hole2.hh, 1, hole2.beta_auto
			       .divergence(hole2.ff).derive_cov(hole2.ff), 0) 
      - hdirac2.derive_lie(hole2.beta_auto) 
      + hole2.gamt_point.divergence(hole2.ff)*hole2.decouple ;      
    vtmp_deux.inc_dzpuis() ; // dzpuis: 3 -> 4

    source_deux -= vtmp_deux ; 
        
    source_deux += 2./3.* hole2.beta_auto.divergence(hole2.ff) 
	* hdirac2 ;

    source_deux.std_spectral_base() ;
    

    Vector source_1 (source_un) ;
    Vector source_2 (source_deux) ;
    source_1.change_triad(hole1.mp.get_bvect_cart()) ;
    source_2.change_triad(hole2.mp.get_bvect_cart()) ;

    cout << "source shift_x" << endl << norme(source_1(1)) << endl ;
    cout << "source shift_y" << endl << norme(source_1(2)) << endl ;
    cout << "source shift_z" << endl << norme(source_1(3)) << endl ;
    
    // Filter for high frequencies.
    for (int i=1 ; i<=3 ; i++) {
	source_un.set(i).filtre(4) ;
	source_deux.set(i).filtre(4) ;
    }
    
    // Boundary conditions 
    // --------------------

    Valeur lim_x_un (hole1.mp.get_mg()-> get_angu()) ;
    Valeur lim_y_un (hole1.mp.get_mg()-> get_angu()) ;
    Valeur lim_z_un (hole1.mp.get_mg()-> get_angu()) ;

    Valeur lim_x_deux (hole2.mp.get_mg()-> get_angu()) ;
    Valeur lim_y_deux (hole2.mp.get_mg()-> get_angu()) ;
    Valeur lim_z_deux (hole2.mp.get_mg()-> get_angu()) ;

    switch (bound_beta) {

	case 0 : {
	    
	    lim_x_un = hole1.boundary_beta_x(omega, omega_eff) ;
	    lim_y_un = hole1.boundary_beta_y(omega, omega_eff) ;
	    lim_z_un = hole1.boundary_beta_z() ;
	    
	    lim_x_deux = hole2.boundary_beta_x(omega, omega_eff) ;
	    lim_y_deux = hole2.boundary_beta_y(omega, omega_eff) ;
	    lim_z_deux = hole2.boundary_beta_z() ;
	    break ;
	}

	default : {
	    cout << "Unexpected type of boundary conditions for beta!" 
		 << endl 
		 << "  bound_beta = " << bound_beta << endl ; 
	    abort() ;
	    break ; 
	}
	    
    } // End of switch  


    // We solve :
    // -----------

    Vector beta_un_old (hole1.beta_auto) ;
    Vector beta_deux_old (hole2.beta_auto) ;
    Vector beta1 (hole1.beta_auto) ;
    Vector beta2 (hole2.beta_auto) ;
    
    poisson_vect_binaire (1./3., source_un, source_deux, 
	lim_x_un, lim_y_un, lim_z_un, 
	lim_x_deux, lim_y_deux, lim_z_deux, 
	beta1, beta2, 0, precision) ;
 
 
    beta1.change_triad(hole1.mp.get_bvect_cart()) ;
    beta2.change_triad(hole2.mp.get_bvect_cart()) ;

    for (int i=1 ; i<=3 ; i++) {
	beta1.set(i).raccord(3) ;
	beta2.set(i).raccord(3) ;
    }

    cout << "shift_auto x" << endl << norme(beta1(1)) << endl ;
    cout << "shift_auto y" << endl << norme(beta1(2)) << endl ;
    cout << "shift_auto z" << endl << norme(beta1(3)) << endl ;

    beta1.change_triad(hole1.mp.get_bvect_spher()) ;
    beta2.change_triad(hole2.mp.get_bvect_spher()) ;

    // Check: has the Poisson equation been correctly solved ?
    // -----------------------------------------------------
    
    int nz = hole1.mp.get_mg()->get_nzone() ;
    Vector lap_beta = (beta1.derive_con(hole1.ff)).divergence(hole1.ff) 
	+ 1./3.* beta1.divergence(hole1.ff).derive_con(hole1.ff) ;
    source_un.dec_dzpuis() ;

    Tbl tdiff_beta_r = diffrel(lap_beta(1), source_un(1)) ; 
    Tbl tdiff_beta_t = diffrel(lap_beta(2), source_un(2)) ; 
    Tbl tdiff_beta_p = diffrel(lap_beta(3), source_un(3)) ; 
    
    cout << 
	"Relative error in the resolution of the equation for beta : "
	 << endl ; 
    cout << "r component : " ;
    for (int l=0; l<nz; l++) {
	cout << tdiff_beta_r(l) << "  " ; 
    }
    cout << endl ;
    cout << "t component : " ;
    for (int l=0; l<nz; l++) {
	cout << tdiff_beta_t(l) << "  " ; 
    }
    cout << endl ;
    cout << "p component : " ;
    for (int l=0; l<nz; l++) {
	cout << tdiff_beta_p(l) << "  " ; 
    }
    cout << endl ;
    
    
    // Relaxation
    // -----------

    Vector beta1_new (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    Vector beta2_new (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;

    // Construction of Omega d/dphi
    // ----------------------------
    
    Vector omdsdp1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
    Scalar yya1 (hole1.mp) ;
    yya1 = hole1.mp.ya ;
    Scalar xxa1 (hole1.mp) ;
    xxa1 = hole1.mp.xa ;
    
    if (fabs(hole1.mp.get_rot_phi()) < 1e-10){ 
      omdsdp1.set(1) = - omega * yya1 ;
      omdsdp1.set(2) = omega * xxa1 ;
      omdsdp1.set(3).annule_hard() ;
    }
    else{
      omdsdp1.set(1) = omega * yya1 ;
      omdsdp1.set(2) = - omega * xxa1 ;
      omdsdp1.set(3).annule_hard() ;
    }
    
    omdsdp1.set(1).set_spectral_va()
      .set_base(*(hole1.mp.get_mg()->std_base_vect_cart()[0])) ;
    omdsdp1.set(2).set_spectral_va()
      .set_base(*(hole1.mp.get_mg()->std_base_vect_cart()[1])) ;
    omdsdp1.set(3).set_spectral_va()
      .set_base(*(hole1.mp.get_mg()->std_base_vect_cart()[2])) ;
    
    omdsdp1.annule_domain(nz-1) ;
    


    Vector omdsdp2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
    Scalar yya2 (hole2.mp) ;
    yya2 = hole2.mp.ya ;
    Scalar xxa2 (hole2.mp) ;
    xxa2 = hole2.mp.xa ;
    
    if (fabs(hole2.mp.get_rot_phi()) < 1e-10){ 
      omdsdp2.set(1) = - omega * yya2 ;
      omdsdp2.set(2) = omega * xxa2 ;
      omdsdp2.set(3).annule_hard() ;
    }
    else{
      omdsdp2.set(1) = omega * yya2 ;
      omdsdp2.set(2) = - omega * xxa2 ;
      omdsdp2.set(3).annule_hard() ;
    }
    
    omdsdp2.set(1).set_spectral_va()
      .set_base(*(hole2.mp.get_mg()->std_base_vect_cart()[0])) ;
    omdsdp2.set(2).set_spectral_va()
      .set_base(*(hole2.mp.get_mg()->std_base_vect_cart()[1])) ;
    omdsdp2.set(3).set_spectral_va()
      .set_base(*(hole2.mp.get_mg()->std_base_vect_cart()[2])) ;
    
    omdsdp2.annule_domain(nz-1) ;
    


    // New shift
    beta1_new = relax*(beta1+hole1.decouple*omdsdp1) + (1-relax)*beta_un_old ;
    beta2_new = relax*(beta2+hole2.decouple*omdsdp2) + (1-relax)*beta_deux_old ;

    hole1.beta_auto = beta1_new ;
    hole2.beta_auto = beta2_new ;

    hole1.beta_comp_import(hole2) ;
    hole2.beta_comp_import(hole1) ;

    // Regularisation of the shifts if necessary
    // -----------------------------------------

    int nnt = hole1.mp.get_mg()->get_nt(1) ;
    int nnp = hole1.mp.get_mg()->get_np(1) ;
    
    int check ;
    check = 0 ;
    for (int k=0; k<nnp; k++)
	for (int j=0; j<nnt; j++){
	    if (fabs((hole1.n_auto+hole1.n_comp).val_grid_point(1, k, j , 0)) < 1e-8){
		check = 1 ;
		break ;
	    }
	}

    if (check == 1){
	double diff_un = hole1.regularisation (hole1.beta_auto, 
					 hole2.beta_auto, omega) ;
	double diff_deux = hole2.regularisation (hole2.beta_auto, 
					   hole1.beta_auto, omega) ;
	hole1.regul = diff_un ;
	hole2.regul = diff_deux ;
    }
    
    else {
	hole1.regul = 0. ;
	hole2.regul = 0. ;
    }
    
}
