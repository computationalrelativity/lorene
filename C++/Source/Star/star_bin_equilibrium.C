
/*
 *   Method of class Star_bin to compute an equilibrium configuration
 *
 *  (see file star.h for documentation).
 */
/*
 *   Copyright (c) 2004 Francois Limousin
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

char star_bin_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.18  2005/02/24 16:04:13  f_limousin
 * Change the name of some variables (for instance dcov_logn --> dlogn).
 * Improve the resolution of the tensorial poisson equation for hh.
 *
 * Revision 1.17  2005/02/18 13:14:18  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.16  2005/02/17 17:32:53  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.15  2005/02/11 18:13:47  f_limousin
 * Important modification : all the poisson equations for the metric
 * quantities are now solved on an affine mapping.
 *
 * Revision 1.14  2004/12/17 16:23:19  f_limousin
 * Modif. comments.
 *
 * Revision 1.13  2004/06/22 12:49:12  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.12  2004/05/27 12:41:00  p_grandclement
 * correction of some shadowed variables
 *
 * Revision 1.11  2004/05/25 14:18:00  f_limousin
 * Include filters
 *
 * Revision 1.10  2004/05/10 10:26:22  f_limousin
 * Minor changes to avoid warnings in the compilation of Lorene
 *
 * Revision 1.9  2004/04/08 16:32:48  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.8  2004/03/25 10:29:26  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.7  2004/03/23 09:56:09  f_limousin
 * Many minor changes
 *
 * Revision 1.6  2004/02/27 21:16:32  e_gourgoulhon
 * Function contract_desal replaced by function contract with
 * argument desaliasing set to true.
 *
 * Revision 1.5  2004/02/27 09:51:51  f_limousin
 * Many minor changes.
 *
 * Revision 1.4  2004/02/21 17:05:13  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.3  2004/01/20 15:17:48  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "cmp.h"
#include "tenseur.h"
#include "metrique.h"
#include "star.h"
#include "param.h"
#include "graphique.h"
#include "utilitaires.h"
#include "tensor.h"
#include "nbr_spx.h"
#include "unites.h"


void Star_bin::equilibrium(double ent_c, int mermax, int mermax_potvit, 
			   int , double , 
			   double relax_potvit, double thres_adapt,
			   const Tbl& fact_resize, Tbl& diff, int ,
			   double omega) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;
    
    // Initializations
    // ---------------
    
    const Mg3d* mg = mp.get_mg() ; 
    int nz = mg->get_nzone() ;	    // total number of domains
    int nr = mp.get_mg()->get_nr(0) ;
    int nt = mp.get_mg()->get_nt(0) ;
    int np = mp.get_mg()->get_np(0) ;
    
    // The following is required to initialize mp_prev as a Map_et:
    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ; 
    
    // Domain and radial indices of points at the surface of the star:
    int l_b = nzet - 1 ; 
    int i_b = mg->get_nr(l_b) - 1 ; 
    int k_b ;
    int j_b ; 
 
    // Value of the enthalpy defining the surface of the star
    double ent_b = 0 ; 
    
    // Error indicators
    // ----------------
    
    double& diff_ent = diff.set(0) ; 
    double& diff_vel_pot = diff.set(1) ; 
    double& diff_logn = diff.set(2) ; 
    double& diff_lnq = diff.set(3) ; 
    double& diff_beta_r = diff.set(4) ; 
    double& diff_beta_t = diff.set(5) ; 
    double& diff_beta_p = diff.set(6) ; 
    double& diff_h11 = diff.set(7) ; 
    double& diff_h21 = diff.set(8) ; 
    double& diff_h31 = diff.set(9) ; 
    double& diff_h22 = diff.set(10) ; 
    double& diff_h32 = diff.set(11) ; 
    double& diff_h33 = diff.set(12) ; 



    // Parameters for the function Map_et::adapt
    // -----------------------------------------

    Param par_adapt ;
    int nitermax = 100 ;
    int niter ; 
    int adapt_flag = 1 ;    //  1 = performs the full computation, 
    //  0 = performs only the rescaling by 
    //      the factor alpha_r
    //##    int nz_search = nzet + 1 ;  // Number of domains for searching the 
    // enthalpy
    int nz_search = nzet ;	// Number of domains for searching the enthalpy
				//  isosurfaces

    double precis_secant = 1.e-14 ; 
    double alpha_r ; 
    double reg_map = 1. ; // 1 = regular mapping, 0 = contracting mapping

    Tbl ent_limit(nz) ; 


    par_adapt.add_int(nitermax, 0) ; // maximum number of iterations to 
    // locate zeros by the secant method
    par_adapt.add_int(nzet, 1) ;    // number of domains where the adjustment 
    // to the isosurfaces of ent is to be 
    // performed
    par_adapt.add_int(nz_search, 2) ;	// number of domains to search 	
                                        // the enthalpy isosurface
    par_adapt.add_int(adapt_flag, 3) ; //  1 = performs the full computation, 
    //  0 = performs only the rescaling by 
    //      the factor alpha_r
    par_adapt.add_int(j_b, 4) ; //  theta index of the collocation point 
    //  (theta_*, phi_*)
    par_adapt.add_int(k_b, 5) ; //  theta index of the collocation point 
    //  (theta_*, phi_*)

    par_adapt.add_int_mod(niter, 0) ;  // number of iterations actually used in
    //  the secant method
    
    par_adapt.add_double(precis_secant, 0) ; // required absolute precision in 
    // the determination of zeros by 
    // the secant method
    par_adapt.add_double(reg_map, 1)	;  // 1. = regular mapping, 
                                           // 0 = contracting mapping
    
    par_adapt.add_double(alpha_r, 2) ;	    // factor by which all the radial 
					    // distances will be multiplied 
    	   
    par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces. 
 
        
    Cmp ssjm1logn (ssjm1_logn) ;
    Cmp ssjm1lnq (ssjm1_lnq) ;
    Cmp ssjm1khi (ssjm1_khi) ;
    Tenseur ssjm1wshift(mp, 1, CON, mp.get_bvect_cart()) ;
    ssjm1wshift.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      ssjm1wshift.set(i) = Cmp(ssjm1_wshift(i+1)) ;
    }
    
    double precis_poisson = 1.e-16 ;     

    // Parameters for the function Scalar::poisson for logn_auto
    // ---------------------------------------------------------------
 
    Param par_logn ;    

    par_logn.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_logn.add_double(relax_poisson,  0) ; // relaxation parameter
    par_logn.add_double(precis_poisson, 1) ; // required precision
    par_logn.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_logn.add_cmp_mod( ssjm1logn ) ; 
    
    // Parameters for the function Scalar::poisson for lnq_auto
    // ---------------------------------------------------------------
    
    Param par_lnq ; 
    
    par_lnq.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_lnq.add_double(relax_poisson,  0) ; // relaxation parameter
    par_lnq.add_double(precis_poisson, 1) ; // required precision
    par_lnq.add_int_mod(niter, 0) ; // number of iterations actually used -
    par_lnq.add_cmp_mod( ssjm1lnq ) ; 
 
    // Parameters for the function Vector::poisson for shift method 2 
    // ---------------------------------------------------------------
    
    Param par_beta ; 
    
    par_beta.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_beta.add_double(relax_poisson,  0) ; // relaxation parameter
    par_beta.add_double(precis_poisson, 1) ; // required precision
    par_beta.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_beta.add_cmp_mod(ssjm1khi) ; 
    par_beta.add_tenseur_mod(ssjm1wshift) ; 
  

    // Computation of all derived quantities. It must be performed before
    // the adaptation of the mapping to avoid discontinuities. 
    // Derivatives of N and logN
    //--------------------------
    
    Vector dlogn_auto = logn_auto.derive_cov(flat) ;
    Vector dlogn_auto_u = logn_auto.derive_con(flat) ;
    Vector dlogn_u = logn.derive_con(flat) ;
    
    const Vector& dnn = nn.derive_cov(flat) ;
    const Tensor& dnn_dd = dnn.derive_cov(flat) ;
    
    // Derivatives of lnq, Q  and beta 
    //----------------------------------
    
    Vector dlnpsi_auto = 0.5*(lnq_auto - logn_auto).derive_cov(flat) ;
    Vector dlnpsi_auto_u = 0.5*(lnq_auto - logn_auto).derive_con(flat) ;
    Vector dlnpsi_u = 0.5*(lnq - logn).derive_con(flat) ;
    
    const Vector& dlnq_u = lnq.derive_con(flat) ;    
    const Vector& dlnq_auto_u = lnq_auto.derive_con(flat) ;
    
    const Tensor& dbeta = beta.derive_cov(flat) ;
    Tensor dbeta_dd = dbeta.derive_cov(flat) ;
    dbeta_dd.inc_dzpuis() ;

    Scalar psi2 (pow(psi4, 0.5)) ;
    psi2.std_spectral_base() ;
    
    Scalar qq = exp(lnq) ;
    qq.std_spectral_base() ;
    
    const Vector& dqq = qq.derive_cov(flat) ;
    Tensor dqq_dd = dqq.derive_cov(flat) ;
    dqq_dd.inc_dzpuis() ;
    Tensor dqq_ud = dqq.derive_con(flat) ;
    dqq_ud.inc_dzpuis() ;
    
    
    // Derivatives of hh, gtilde... 
    //------------------------------
    
    const Tensor& dhh = hh.derive_cov(flat) ;
    const Tensor& dhh_auto = hh_auto.derive_cov(flat) ;
    const Tensor& dhh_auto_u = hh_auto.derive_con(flat) ;
    
    Tensor dhh_dd = dhh.derive_cov(flat) ;
    dhh_dd.inc_dzpuis() ;
    
    Sym_tensor gtilde_cov = gtilde.cov() ;
    Sym_tensor gtilde_con = gtilde.con() ;
    const Tensor& dgtilde = gtilde_cov.derive_cov(flat) ;
    
    Connection gamijk (gtilde, flat) ;
    const Tensor& deltaijk = gamijk.get_delta() ;
    
    // Declaration of all sources 
    //---------------------------
    
    Scalar source_Hij(mp) ;
    Scalar source_Qij(mp) ;
    Scalar source_Sij(mp) ;
    Sym_tensor source_hh(mp, CON, mp.get_bvect_spher()) ;
    
    
    Tensor source_1 (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source_2 (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source_3 (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source_4 (mp, 2, CON, mp.get_bvect_spher()) ;
    
    
    Tensor source1_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source2_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source3_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source4_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source5_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source6_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source7_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
    
    
    Tensor source1_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source2_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source3_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source4_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source5_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source6_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source7_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source8_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source9_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source10_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source11_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
    
    
    Tensor source1_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source2_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source3_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source4_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source5_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source6_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
    Tensor source7_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
	    
	    
    // Source for hh
    //--------------
    
    source_1 = contract(dhh_auto_u, 1, dlogn + 2*dlnpsi, 0) ;
    
    source_2 = - contract(dhh_auto_u, 2, dlogn + 2*dlnpsi, 0) ;
    
    source_3 =  aa_auto.derive_lie(beta) ;
    source_3.inc_dzpuis() ;
    
    source_4 = 0.666666666666667 * beta.divergence(flat) * aa_auto ;
    
    
    // Source terms for Hij
    //---------------------
    
    source1_Hij = 8.*nn*contract(hh_auto,1,dlnpsi,0)*dlnpsi_u ;
    
    source2_Hij = - contract(hh_auto, 1 , dqq_ud, 0) / psi2 ;
    
    source3_Hij = 4.*nn*contract(hh_auto,1,dlogn,0)*dlnpsi_u ; 
    
    source4_Hij = 4.*nn*contract(hh_auto,1,dlnpsi,0)*dlogn_u ;
    
    source5_Hij = ( contract(hh_auto, 0, 1, dqq_dd, 0, 1)/psi2
		    - 8. * nn * contract(contract(hh_auto,0
						  , dlnpsi, 0), 0, dlogn, 0)
		    - 8. * nn * contract(contract(hh_auto,0
						  , dlnpsi, 0), 0, dlnpsi, 0))
	* flat.con() * 0.3333333333333333 ; 
    
    source6_Hij = ( qq.laplacian() / psi2 
		    - 8. * nn * contract(dlnpsi, 0, dlogn_u, 0)
		    - 8. * nn * contract(dlnpsi, 0, dlnpsi_u, 0))
	* hh_auto * 0.3333333333333333 ;
    
    source7_Hij = 0.6666666666666667 * qpig * nn * s_euler * hh ;
    
    
    // Source_terms for Qij
    //---------------------
    
    
    source1_Qij = 0.5*nn*contract(hh_auto, 0, 1,dhh_dd, 2, 3) ;
    
    
    source2_Qij = 0.5 * nn * ( - contract(contract(dhh_auto, 1,  
						   dhh, 2), 1, 3)
		       - contract(contract(contract(gtilde_con, 0, 
		      	    dhh_auto, 2), 0, dhh, 2), 1, 3, gtilde_cov, 0, 1)
		       + contract(contract(contract(contract(gtilde_con, 1, 
		       	     dhh_auto, 2), 2, dhh, 2), 1, gtilde_cov, 0), 2, 3)
		       + contract(contract(dgtilde, 0, 1, dhh_auto, 0, 1), 
		       	  0, 1, gtilde_con * gtilde_con, 1, 3)*0.5 ) ;
	      
    source3_Qij = 0.5 * nn * contract(contract(contract(
			      contract(gtilde_con, 1, dhh_auto, 2), 1,
				   dhh, 2), 1, gtilde_cov, 1), 2, 3) ;
    
    source4_Qij = nn * 8. * contract(contract(hh_auto, 1, 
				    dlnpsi, 0)*hh, 2, dlnpsi, 0) ;
	    
    source5_Qij = - contract(contract(hh_auto, 1, 
				      dqq_dd, 1), 1, hh, 1) / psi2 ;
    
    source6_Qij = - 0.5 * contract(contract(hh_auto, 1, dhh, 2)
				   , 2, dqq, 0) / psi2 ;
    
    
    source7_Qij = 0.5 * contract(contract(hh_auto, 1, dhh, 2), 
				 0, dqq, 0) /psi2 ;
    
    
    source8_Qij = 4. * nn * contract(contract(hh_auto, 1, 
					      dlnpsi, 0)*hh, 2, dlogn, 0)
	
	+ 4. * nn * contract(contract(hh_auto, 1, 
				      dlogn, 0)*hh, 2, dlnpsi, 0) ;
    
    source9_Qij = nn * (contract(contract(dhh_auto, 0, 1, 
			      dgtilde, 0, 2), 0, 1, gtilde_con, 0, 1)
		     -  contract(contract(dhh_auto, 0, 1, 
			 dgtilde, 0, 1), 0, 1, gtilde_con, 0, 1) * 0.5 )
	                 * 0.1666666666666667 * gtilde_con ;

	    
    source10_Qij = - 2.666666666666667*nn*( contract(contract(hh, 0, 
			         dlnpsi, 0), 0, dlnpsi, 0) * hh_auto 
	        	    + contract(contract(hh, 0, 
	      	       	dlnpsi, 0), 0, dlogn, 0) * hh_auto ) ;

    
    source11_Qij = contract(hh, 0, 1, dqq_dd, 0, 1) / psi2
	* hh_auto * 0.3333333333333333 ;
    
    
    // Source terms for Sij
    //---------------------
    
    
    source1_Sij = 8. * nn * dlnpsi_auto_u * dlnpsi_u  ;
    
    source2_Sij = - qq.derive_con(flat).derive_con(flat) / psi2 * 0.5 ;
    source2_Sij.inc_dzpuis(1) ;
    
    source3_Sij = 4. * nn * dlnpsi_u * dlogn_auto_u ;
    
    source4_Sij = 0.1666666666666667 * qq.laplacian()/psi2*flat.con() ;
    
    source5_Sij = - ( nn * contract(dlnpsi, 0, 
				    dlogn_auto_u , 0)
		      + nn * contract(dlnpsi_auto, 0, dlnpsi_u, 0))
	* 2.666666666666667 * flat.con() ;
    
    source6_Sij = 2. * nn * contract(contract(gtilde_cov, 0, 
			         aa_auto + aa_comp, 1), 0, aa_auto, 1) ;
    
    
    source7_Sij = - 2. * qpig * nn * ( psi4 * stress_euler 
			  - 0.33333333333333333 * s_euler * flat.con() ) ; 
    
    
	    Coord& sin_theta = mp.sint ;
	    Scalar sint (mp) ;
	    sint = sin_theta ;
	    sint.std_spectral_base() ;

	    Sym_tensor source_dirac (mp, CON, mp.get_bvect_spher()) ;
	    // killing vector / rr
	    Vector killing (mp, CON, mp.get_bvect_spher()) ;
	    killing.set(1) = 0. ;
	    killing.set(2) = 0. ;
	    killing.set(3) = omega*sint ;
	    killing.std_spectral_base() ;
	    source_dirac =  - aa_auto.derive_lie(killing) ;
	    for(int i=1; i<=3; i++) 
		for(int j=1; j<=i; j++) 
		    source_dirac.set(i,j).mult_r() ;
	    source_dirac.inc_dzpuis() ;

	    for(int i=1; i<=3; i++) 
		for(int j=1; j<=i; j++) {

		    source_Hij = ( source1_Hij(i,j) + source1_Hij(j,i) +
			           source2_Hij(j,i) + source2_Hij(i,j) +
				   source3_Hij(j,i) + source3_Hij(i,j) +
				   source4_Hij(j,i) + source4_Hij(i,j) +
				   source5_Hij(i,j) + source6_Hij(i,j)) / psi4
			+ source7_Hij(i,j) ;


		    source_Qij = ( source1_Qij(i,j) + source2_Qij(i,j) +
				   source3_Qij(j,i) + source4_Qij(i,j) +
				   source5_Qij(i,j) + source6_Qij(i,j) +
				   source6_Qij(j,i) + source7_Qij(i,j) +
				   source8_Qij(i,j) + source9_Qij(i,j) +
				   source10_Qij(i,j) + source11_Qij(i,j)) 
			/ psi4 ;
		   

		    source_Sij = ( 0*source1_Sij(i,j) + 0*source2_Sij(i,j) +
				   0*source3_Sij(i,j) + 0*source3_Sij(j,i) +
				   0*source4_Sij(i,j) + 0*source5_Sij(i,j)) / psi4
			+ 0*source6_Sij(i,j) + source7_Sij(i,j) ;
			    

		    source_hh.set(i,j) = 0*source_1(i,j) + 0*source_1(j,i) 
			+ 0*source_2(i,j) - 2 * psi4 / nn * (
			    0*source_3(i,j) + 0*source_4(i,j) +
			    0*source_Hij + 0*source_Qij + source_Sij ) 
			+0 * source_dirac(i,j) ;
		    

		    cout << "i = " << i << ", j = " << j << endl ;
		    cout << "----------------------------" << endl ;
		    
		    cout << "source1" << endl << norme(source_1(i,j)/(nr*nt*np)) << endl ;
		    cout << "source2" << endl << norme(source_2(i,j)/(nr*nt*np)) << endl ;
		    cout << "source3" << endl << norme(source_3(i,j)/(nr*nt*np)) << endl ;
		    cout << "source4" << endl << norme(source_4(i,j)/(nr*nt*np)) << endl ;
		    cout << "source H" << endl << norme(source_Hij/(nr*nt*np)) << endl ;
		    cout << "source Q" << endl << norme(source_Qij/(nr*nt*np)) << endl ;
		    cout << "source S" << endl << norme(source_Sij/(nr*nt*np)) << endl ;
		    cout << "source dirac" << endl << norme(source_dirac(i,j)/(nr*nt*np)) << endl ;
		    cout << "source tot" << endl << norme(source_hh(i,j)/(nr*nt*np)) << endl ;
		    
		} 
  
  	    Sym_tensor source_ana (mp, CON, mp.get_bvect_spher()) ;
	    source_ana = dlogn_auto_u * dlnq_u - 0.3333333333333333 *
		contract(dlogn_auto, 0, dlnq_u, 0) * flat.con() ;

//		    des_meridian(source_hh, 0., 4., "source_hh", 20) ;

    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // ------------------
    
    Scalar pot_ext = logn_comp + pot_centri + loggam ;
    
    Scalar ent_jm1 = ent ;	// Enthalpy at previous step
    
    Scalar source(mp) ; // source term in the equation for hh_auto, 
                            // logn_auto and lnq_auto
			    
    Vector source_beta(mp, CON, mp.get_bvect_spher()) ;  // source term 
    // in the equation for beta_auto



    //=========================================================================
    // 			Start of iteration
    //=========================================================================

    for(int mer=0 ; mer<mermax ; mer++ ) {

	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_ent = " << diff_ent << endl ;    

	//-----------------------------------------------------
	// Resolution of the elliptic equation for the velocity
	// scalar potential
	//-----------------------------------------------------
	
	double precis_poisson = 1e-16 ;
	if (irrotational) {
	    diff_vel_pot = velocity_potential(mermax_potvit, precis_poisson, 
					      relax_potvit) ; 
	    
	}

	diff_vel_pot = 0. ; // to avoid the warning 

	//-----------------------------------------------------
	// Computation of the new radial scale
	//--------------------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface
	
	// Values at the center of the star:
	double logn_auto_c  = logn_auto.val_grid_point(0, 0, 0, 0) ; 
	double pot_ext_c  = pot_ext.val_grid_point(0, 0, 0, 0) ; 

	// Search for the reference point (theta_*, phi_*) [notation of
	//  Bonazzola, Gourgoulhon & Marck PRD 58, 104020 (1998)]
	//  at the surface of the star
	// ------------------------------------------------------------
	double alpha_r2 = 0 ; 
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {
		
		double pot_ext_b  = pot_ext.val_grid_point(l_b, k, j, i_b) ; 
		double logn_auto_b  = logn_auto.val_grid_point(l_b, k, j, i_b) ;
		// See Eq (100) from Gourgoulhon et al. (2001)
		double alpha_r2_jk = ( ent_c - ent_b + pot_ext_c - pot_ext_b) /
 
		    ( logn_auto_b - logn_auto_c ) ;
		  
		if (alpha_r2_jk > alpha_r2) {
		    alpha_r2 = alpha_r2_jk ; 
		    k_b = k ; 
		    j_b = j ; 
		}

	    }
	}
      	
	alpha_r = sqrt(alpha_r2) ;
		
	cout << "k_b, j_b, alpha_r: " << k_b << "  " << j_b << "  " 
	     <<  alpha_r << endl ;

	// New value of logn_auto 
	// ----------------------

	logn_auto = alpha_r2 * logn_auto ;
	logn_auto_c  = logn_auto.val_grid_point(0, 0, 0, 0) ;

	//------------------------------------------------------------
	// Change the values of the inner points of the second domain
	// by those of the outer points of the first domain
	//------------------------------------------------------------

	logn_auto.set_spectral_va().smooth(nzet, logn_auto.set_spectral_va()) ;

	//------------------------------------------
	// First integral	-->  enthalpy in all space
	// See Eq (98) from Gourgoulhon et al. (2001)
	//-------------------------------------------

	ent = (ent_c + logn_auto_c + pot_ext_c) - logn_auto - pot_ext ;

	(ent.set_spectral_va()).smooth(nzet, ent.set_spectral_va()) ;

	//----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//----------------------------------------------------
	
	// Shall the adaptation be performed (cusp) ?
	// ------------------------------------------
	
	double dent_eq = ent.dsdr().val_point(ray_eq(),M_PI/2.,0.) ;
	double dent_pole = ent.dsdr().val_point(ray_pole(),0.,0.) ;
	double rap_dent = fabs( dent_eq / dent_pole ) ; 
	cout << "| dH/dr_eq / dH/dr_pole | = " << rap_dent << endl ; 
	
	if ( rap_dent < thres_adapt ) {
	    adapt_flag = 0 ;	// No adaptation of the mapping 
	    cout << "******* FROZEN MAPPING  *********" << endl ; 
	}
	else{
	    adapt_flag = 1 ;	// The adaptation of the mapping is to be
	    //  performed
	}

	ent_limit.set_etat_qcq() ; 
	for (int l=0; l<nzet; l++) {	// loop on domains inside the star
	    ent_limit.set(l) = ent.val_grid_point(l, k_b, j_b, i_b) ; 
	}
	ent_limit.set(nzet-1) = ent_b  ; 

	Map_et mp_prev = mp_et ; 

	Cmp ent_cmp(ent) ;
	mp.adapt(ent_cmp, par_adapt) ; 
	ent = ent_cmp ;

	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	Cmp ent_cmp2 (ent) ;
	mp.reevaluate_symy(&mp_prev, nzet+1, ent_cmp2) ; 
	ent = ent_cmp2 ;

	double ent_s_max = -1 ; 
	int k_s_max = -1 ; 
	int j_s_max = -1 ; 
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {
		double xx = fabs( ent.val_grid_point(l_b, k, j, i_b) ) ;
		if (xx > ent_s_max) {
		    ent_s_max = xx ; 
		    k_s_max = k ; 
		    j_s_max = j ; 
		}
	    }
	}
	cout << "Max. abs(enthalpy) at the boundary between domains nzet-1"
	     << " and nzet : " << endl ; 
	cout << "   " << ent_s_max << " reached for k = " << k_s_max <<
	    " and j = " << j_s_max << endl ; 
	

	// Readjustment of the external boundary of domain l=nzet
	// to keep a fixed ratio with respect to star's surface

	
    	if (nz == 4 && nzet == 1) {
	    double rr_in_1 = mp.val_r(1,-1., M_PI/2, 0.) ; 
	    double rr_out_1 = mp.val_r(1, 1., M_PI/2, 0.) ; 
	    double rr_out_2 = mp.val_r(2, 1., M_PI/2, 0.) ; 

	    mp.resize(1, rr_in_1/rr_out_1 * fact_resize(0)) ; 
	    mp.resize(2, rr_in_1/rr_out_2 * fact_resize(1)) ; 
	}
	else{

	    if (nz == 5 && nzet == 1) {
		double rr_in_1 = mp.val_r(1,-1., M_PI/2, 0.) ; 
		double rr_out_1 = mp.val_r(1, 1., M_PI/2, 0.) ; 
		double rr_out_2 = mp.val_r(2, 1., M_PI/2, 0.) ; 
		double rr_out_3 = mp.val_r(3, 1., M_PI/2, 0.) ; 
		
		double fact_resize_0 ;
		if (fact_resize(0) > 2.4) fact_resize_0 = fact_resize(0)/2. ;
		else fact_resize_0 = fact_resize(0)/2. + 0.5 ;


		mp.resize(1, rr_in_1/rr_out_1 * fact_resize_0) ; 
		mp.resize(2, rr_in_1/rr_out_2 * fact_resize(0)) ; 
		mp.resize(3, rr_in_1/rr_out_3 * fact_resize(1)) ; 
	    }
	    else{		
		int n_resize ;
		//      	if (nz > 4) {
		//       	  n_resize = nz - 4 ;
		if (nz > 4) {
		    n_resize = nz - 3 ;
		}
		else {
		    n_resize = nzet ;
		}
		
		double rr_in = mp.val_r(nzet,-1., M_PI/2, 0.) ; 
		double rr_out = mp.val_r(n_resize,1., M_PI/2, 0.) ; 
		
		mp.resize(n_resize, rr_in/rr_out * fact_resize(0)) ; 
	    }
	}
	

	//----------------------------------------------------
	// Equation of state  
	//----------------------------------------------------
	
	equation_of_state() ; 	// computes new values for nbar (n), ener (e) 
				// and press (p) from the new ent (H)
	
	//---------------------------------------------------------
	// Matter source terms in the gravitational field equations	
	//---------------------------------------------------------

	hydro_euler() ;		// computes new values for ener_euler (E), 
				// s_euler (S) and u_euler (U^i)

	
	//--------------------------------------------------------
	// Poisson equation for logn_auto (nu_auto)
	//--------------------------------------------------------
  
	// Source 
	//--------
    
	source = qpig * psi4 % (ener_euler + s_euler) 
	    + psi4 % (aa_quad_auto + aa_quad_comp) ;

	source -= contract(dlogn_auto_u, 0, dlogn, 0, true) 
	    + 2. * contract(dlogn_auto_u, 0, dlnpsi, 0, true) ;
		    
	Scalar temp = - contract(hh_auto, 0, 1, dnn_dd/nn, 0, 1) ;
	temp.inc_dzpuis(1) ;
	temp.annule_domain(nz-1) ;

	source += temp - 2. * contract(hh_auto, 0, 1, dlnpsi 
				  * dlogn, 0, 1) ;

	cout << "moyenne de la source pour logn_auto" << endl ;
	cout <<  norme(source/(nr*nt*np)) << endl ;

//	source_tot.filtre(4) ;

	// Resolution of the Poisson equation 
	// ----------------------------------

	source.poisson(par_logn, logn_auto) ; 
	ssjm1_logn = ssjm1logn ;

	cout << "logn_auto" << endl << norme(logn_auto/(nr*nt*np)) << endl ;
  
	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	Tbl tdiff_logn = diffrel(logn_auto.laplacian(), source) ;

	cout << 
	    "Relative error in the resolution of the equation for logn_auto : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_logn(l) << "  " ; 
	}
	cout << endl ;
	diff_logn = max(abs(tdiff_logn)) ; 


	//--------------------------------------------------------
	// Poisson equation for lnq_auto
	//--------------------------------------------------------

	// Source
	//--------

	source = qpig * psi4 % s_euler 
	    + 0.75 * psi4 % (aa_quad_auto + aa_quad_comp) ;
	
	source -= 0.5 * contract(dlogn_auto_u, 0, dlogn, 0, true) 
	    + 0.5 * contract(dlnq_auto_u, 0, dlnq, 0, true) ;

	source += 0.0625 * contract(gtilde_con, 0, 1, contract(
					dhh_auto, 0, 1, dgtilde, 0, 1), 0, 1) 
			   
	    - 0.125 * contract(gtilde_con, 0, 1, contract(dhh_auto, 
					0, 1, dgtilde, 0, 2), 0, 1) 

	    + 2. * contract(hh_auto, 0, 1, dlnpsi * dlnpsi, 
			    0, 1) 

	    + 2 * contract(hh_auto, 0, 1, dlnpsi * dlogn, 0, 1) 
	
	    - contract(hh_auto, 0, 1, dqq_dd, 0, 1) / qq ;


	cout << "moyenne de la source pour lnq_auto" << endl ;
	cout <<  norme(source/(nr*nt*np)) << endl ;
	
//	source_tot.filtre(4) ;

	// Resolution of the Poisson equation 
	// ----------------------------------

	source.poisson(par_lnq, lnq_auto) ; 
	ssjm1_lnq = ssjm1lnq ;

	cout << "lnq_auto" << endl << norme(lnq_auto/(nr*nt*np)) << endl ;

	// Check: has the Poisson equation been correctly solved 
	// -----------------------------------------------------
    
	Tbl tdiff_lnq = diffrel(lnq_auto.laplacian(), source) ;

	cout << 
	    "Relative error in the resolution of the equation for lnq : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_lnq(l) << "  " ; 
	}
	cout << endl ;
	diff_lnq = max(abs(tdiff_lnq)) ; 

	//--------------------------------------------------------
	// Vector Poisson equation for beta_auto
	//--------------------------------------------------------

	// Source
	//-------

	u_euler.change_triad(mp.get_bvect_spher()) ;

	source_beta = (4.*qpig) * nn % psi4 % (ener_euler + press) * u_euler 
	
	    + 2. * nn * contract(aa_auto, 1, dlogn - 6 * dlnpsi, 0)  ;


	source_beta -= 2. * nn * contract(aa_auto, 0, 1, deltaijk, 1, 2) 

	    + contract(hh_auto, 0, 1, dbeta_dd, 1, 2) 

	    + 0.3333333333333333 * contract(contract(hh_auto, 1, 
						     dbeta_dd, 2), 1, 2) ;
	

	// Resolution of the Poisson equation 
	// ----------------------------------

	for (int i=1; i<=3; i++) {
	    if(source_beta(i).dz_nonzero()) {
		assert( source_beta(i).get_dzpuis() == 4 ) ; 
	    }
	    else{
		(source_beta.set(i)).set_dzpuis(4) ; 
	    }
	}

	double lambda = double(1) / double(3) ; 

	beta_auto = source_beta.poisson(lambda, par_beta, 2) ; 
	ssjm1_khi = ssjm1khi ;
	for (int i=0; i<3; i++){
	    ssjm1_wshift.set(i+1) = ssjm1wshift(i) ;
	}


	cout << "beta_auto(r)" << endl << norme(beta_auto(1)/(nr*nt*np)) 
	     << endl ;
	cout << "beta_auto(t)" << endl << norme(beta_auto(2)/(nr*nt*np)) 
	     << endl ;
	cout << "beta_auto(p)" << endl << norme(beta_auto(3)/(nr*nt*np)) 
	     << endl ;
  

	// Check: has the equation for beta_auto been correctly solved ?
	// --------------------------------------------------------------
	
	Vector lap_beta = (beta_auto.derive_con(flat)).divergence(flat) 
	    + lambda* beta_auto.divergence(flat).derive_con(flat) ;
	
	source_beta.dec_dzpuis() ;
	Tbl tdiff_beta_r = diffrel(lap_beta(1), source_beta(1)) ; 
	Tbl tdiff_beta_t = diffrel(lap_beta(2), source_beta(2)) ; 
	Tbl tdiff_beta_p = diffrel(lap_beta(3), source_beta(3)) ; 

	cout << 
	    "Relative error in the resolution of the equation for beta_auto : "
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
	
	diff_beta_r = max(abs(tdiff_beta_r)) ; 
	diff_beta_t = max(abs(tdiff_beta_t)) ; 
	diff_beta_p = max(abs(tdiff_beta_p)) ; 


	if (!conf_flat) {
	   
	    //--------------------------------------------------------
	    // Poisson equation for hh
	    //--------------------------------------------------------
	 
		    
	    // Construction of the affine mapping
	    // ----------------------------------
	    
	    double* bornes = new double[nz+1];
	    bornes[0] = 0. ; 
	    for (int l=1; l<nz; l++) 
		bornes[l] = mp.val_r(l, -1., M_PI/2., 0.) ;
	    bornes[nz] = __infinity ; 
	    Map_af mpaff (*mp.get_mg(), bornes) ;
	    mpaff.set_ori(mp.get_ori_x(), mp.get_ori_y(), mp.get_ori_z()) ; 
	    mpaff.set_rot_phi(mp.get_rot_phi()) ;

	    Metric_flat flat_aff (mpaff.flat_met_spher()) ;


	    // Construction of the source on the previous mapping
	    // --------------------------------------------------

	    Sym_tensor source_old_map(mp_et,CON, mp_et.get_bvect_spher()) ;
	    
	    for(int i=1; i<=3; i++)
		for(int j=i; j<=3; j++){
		    source_old_map.set(i,j).set_etat_qcq() ;
		    source_old_map.set(i,j).set_spectral_va() = source_hh(i,j)
			.get_spectral_va() ;
		}
	    
//	    des_meridian(source_old_map, 0., 4., "source_old_map", 0) ;

	    // Construction of the source on the affine mapping
	    // ------------------------------------------------

	    Sym_tensor copie_souhh (source_old_map) ;
	    for(int i=1; i<=3; i++)
		for(int j=i; j<=3; j++)
		    copie_souhh.set(i,j).set_dzpuis(0) ;
	    copie_souhh.std_spectral_base() ;
	    copie_souhh.change_triad(mp.get_bvect_cart()) ;

	    Sym_tensor source_hh_aff (mpaff, CON, mp.get_bvect_cart()) ;
	    for(int i=1; i<=3; i++)
		for(int j=i; j<=3; j++)
		    source_hh_aff.set(i,j).import(copie_souhh(i,j)) ;
	    source_hh_aff.std_spectral_base() ;
	    source_hh_aff.change_triad(mpaff.get_bvect_spher()) ;
	    for(int i=1; i<=3; i++)
		for(int j=i; j<=3; j++)
		    source_hh_aff.set(i,j).set_dzpuis(4) ;
	    
//	    des_meridian(source_hh_aff, 0., 4., "source_hh_aff", 6) ;
	    

	    cout << "Divergence de source_hh_aff" << endl ;
	    for (int i=1; i<=3; i++){
		cout << "  Comp. " << i << " :  " ;
		for (int l=0; l<nz; l++){
		    cout << norme(source_hh_aff.divergence(flat_aff)(i)/(nr*nt*np))(l) << " " ;
		}
		cout << endl ;
	    }
	    cout << endl ;
	    
	    cout << "Trace de la source_hh_aff" << endl ;
	    for (int l=0; l<nz; l++){
		cout << norme(source_hh_aff.trace(flat_aff)/(nr*nt*np))(l) 
		     << " " ;
	    }
	    cout << endl << endl ;


	    // Resolution of the tensorial poisson equation
	    // --------------------------------------------

	    Sym_tensor_trans source_hht(mpaff, mpaff.get_bvect_spher(), flat_aff) ;
	    source_hht = source_hh_aff ;
	    const Sym_tensor_tt& source_htt = source_hht.tt_part() ;
	    
//	    cout << mpaff << endl ;

	    Sym_tensor hh_aff (mpaff, CON, mpaff.get_bvect_spher()) ;

	    des_meridian(source_htt, 1.2, 1.7, "source_htt", 12) ;
	    	    
	    hh_aff = source_htt.poisson(0) ;
	
	    des_meridian(hh_aff, 1.2, 1.7, "hh_aff", 24) ;

	    hh_auto.change_triad(mp.get_bvect_cart()) ;
	    hh_aff.change_triad(mpaff.get_bvect_cart()) ;
	    for(int i=1; i<=3; i++)
		for(int j=i; j<=3; j++)
		    hh_auto.set(i,j).import(hh_aff(i,j)) ;
	    hh_auto.std_spectral_base() ;
	    hh_auto.change_triad(mp.get_bvect_spher()) ;
	    hh_aff.change_triad(mpaff.get_bvect_spher()) ;


//	des_meridian(hh_auto, 0., 5., "hh_auto", 5) ;
	

	    Sym_tensor_tt lap_hh(mpaff, mpaff.get_bvect_spher(), flat_aff) ; 
	    lap_hh = hh_aff.derive_con(flat_aff).divergence(flat_aff) ;
	    lap_hh.inc_dzpuis() ;

//	des_meridian(lap_hh, 0., 5., "lap_hh", 6) ;
	

	    cout << "hh_TT_auto :" << endl ;
	    for (int i=1; i<=3; i++)
		for (int j=1; j<=i; j++) {
		    cout << "  Comp. " << i << " " << j << " :  " ;
		    for (int l=0; l<nz; l++){
			cout << norme(hh_auto(i,j)/(nr*nt*np))(l) << " " ;
		    }
		    cout << endl ;
		}
	    cout << endl ;

	    cout << "lap_hh :" << endl ;
	    for (int i=1; i<=3; i++)
		for (int j=1; j<=i; j++) {
		    cout << "  Comp. " << i << " " << j << " :  " ;
		    for (int l=0; l<nz; l++){
			cout << norme(lap_hh(i,j)/(nr*nt*np))(l) << " " ;
		    }
		    cout << endl ;
		}
	    cout << endl ;
	
	    cout << "source_tt :" << endl ;
	    for (int i=1; i<=3; i++)
		for (int j=1; j<=i; j++) {
		    cout << "  Comp. " << i << " " << j << " :  " ;
		    for (int l=0; l<nz; l++){
			cout << norme(source_htt(i,j)/(nr*nt*np))(l) << " " ;
		    }
		    cout << endl ;
		}
	    cout << endl ;
	
	    // Check: has the Poisson equation been correctly solved ?
	    // -----------------------------------------------------
	
	    cout << "Relative error in the resolution of the equation for hh"
		 << endl ;	
	
	    Tbl tdiff_h11 = diffrel(lap_hh(1,1), source_htt(1,1)) ;
	    cout << " Comp 1 1 :  " ; 
	    for (int l=0; l<nz; l++) {
		cout << tdiff_h11(l) << "  " ; 
	    }
	    cout << endl ;
	    diff_h11 = max(abs(tdiff_h11)) ; 
	
	
	    Tbl tdiff_h21 = diffrel(lap_hh(2,1), source_htt(2,1)) ;
	    cout << " Comp 2 1 :  "  ; 
	    for (int l=0; l<nz; l++) {
		cout << tdiff_h21(l) << "  " ; 
	    }
	    cout << endl ;
	    diff_h21 = max(abs(tdiff_h21)) ; 
	    
	    
	    Tbl tdiff_h31 = diffrel(lap_hh(3,1), source_htt(3,1)) ;
	    cout << " Comp 3 1 :  " ; 
	    for (int l=0; l<nz; l++) {
		cout << tdiff_h31(l) << "  " ; 
	    }
	    cout << endl ;
	    diff_h31 = max(abs(tdiff_h31)) ; 
	
	    Tbl tdiff_h22 = diffrel(lap_hh(2,2), source_htt(2,2)) ;
	    cout << " Comp 2 2 :  " ; 
	    for (int l=0; l<nz; l++) {
		cout << tdiff_h22(l) << "  " ; 
	    }
	    cout << endl ;
	    diff_h22 = max(abs(tdiff_h22)) ; 
	
	
	    Tbl tdiff_h32 = diffrel(lap_hh(3,2) ,source_htt(3,2)) ;
	    cout << " Comp 3 2 :  " ; 
	    for (int l=0; l<nz; l++) {
		cout << tdiff_h32(l) << "  " ; 
	    }
	    cout << endl ;
	    diff_h32 = max(abs(tdiff_h32)) ; 
	
	
	    Tbl tdiff_h33 = diffrel(lap_hh(3,3),source_htt(3,3)) ;
	    cout << " Comp 3 3 :  " ; 
	    for (int l=0; l<nz; l++) {
		cout << tdiff_h33(l) << "  " ;
	    }
	    cout << endl ;
	    diff_h33 = max(abs(tdiff_h33)) ; 
	    
	}

	// End of relativistic equations	
	   
	   
	//-------------------------------------------------
	//  Relative change in enthalpy
	//-------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent, ent_jm1 ) ; 
	diff_ent = diff_ent_tbl(0) ; 
	for (int l=1; l<nzet; l++) {
	    diff_ent += diff_ent_tbl(l) ; 
	}
	diff_ent /= nzet ; 
	
	
	ent_jm1 = ent ; 


    } // End of main loop
    
    //=========================================================================
    // 			End of iteration
    //=========================================================================

}  





    
    
