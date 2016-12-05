/*
 * Methods of Bin_star::helical
 *
 * (see file star.h for documentation)
 *
 */

/*
 *   Copyright (c) 2006 Francois Limousin
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:17:47  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:52:45  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:12:59  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/08/19 06:41:59  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.3  2006/08/01 14:26:50  f_limousin
 * Small changes
 *
 * Revision 1.2  2006/06/05 17:05:57  f_limousin
 * *** empty log message ***
 *
 * Revision 1.1  2006/04/11 14:25:15  f_limousin
 * New version of the code : improvement of the computation of some
 * critical sources, estimation of the dirac gauge, helical symmetry...
 *

 *
 * $Header$ *
 */


// C headers
#include <cmath>

// Lorene headers
#include "cmp.h"
#include "tenseur.h"
#include "metrique.h"
#include "binary.h"
#include "param.h"
#include "graphique.h"
#include "utilitaires.h"
#include "tensor.h"
#include "nbr_spx.h"
#include "unites.h"

namespace Lorene {
void Binary::helical(){

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

 
    Sym_tensor lie_aij_1 (star1.mp, CON, star1.mp.get_bvect_cart()) ;
    Sym_tensor lie_aij_2 (star2.mp, CON, star2.mp.get_bvect_cart()) ;

    Scalar lie_K_1 (star1.mp) ;
    Scalar lie_K_2 (star2.mp) ;

    for (int ll=1; ll<=2; ll++) {
	
	Star_bin star_i (*et[ll-1]) ;
	
	Map& mp = star_i.mp ;
	const Mg3d* mg = mp.get_mg() ; 
	int nz = mg->get_nzone() ;	    // total number of domains
	
	Metric_flat flat = star_i.flat ;
	Metric gtilde = star_i.gtilde ;
	Scalar nn = star_i.nn ;
	Scalar psi4 = star_i.psi4 ;

	// -------------------------------
	// AUXILIARY QUANTITIES
	// -------------------------------

	// Derivatives of N and logN
	//--------------------------

	const Vector dcov_logn_auto = star_i.logn_auto.derive_cov(flat) ;
	
	Tensor dcovdcov_logn_auto = (star_i.logn_auto.derive_cov(flat))
	                                  .derive_cov(flat) ;
	dcovdcov_logn_auto.inc_dzpuis() ;

	// Derivatives of lnq, phi and Q
	//-------------------------------

	const Scalar phi (0.5 * (star_i.lnq - star_i.logn)) ;
	const Scalar phi_auto (0.5 * (star_i.lnq_auto - star_i.logn_auto)) ;

	const Vector dcov_phi_auto = phi_auto.derive_cov(flat) ;
	
	const Vector dcov_lnq = 2*star_i.dcov_phi + star_i.dcov_logn ;
	const Vector dcon_lnq = 2*star_i.dcon_phi + star_i.dcon_logn ;
	const Vector dcov_lnq_auto = star_i.lnq_auto.derive_cov(flat) ;
     	Tensor dcovdcov_lnq_auto = dcov_lnq_auto.derive_cov(flat) ;
	dcovdcov_lnq_auto.inc_dzpuis() ;

	Scalar qq = exp(star_i.lnq) ;
	qq.std_spectral_base() ;
	const Vector& dcov_qq = qq.derive_cov(flat) ;

	Tensor dcovdcov_beta_auto = star_i.beta_auto.derive_cov(flat)
	  .derive_cov(flat) ;
	dcovdcov_beta_auto.inc_dzpuis() ;


	// Derivatives of hij, gtilde... 
	//------------------------------

	Scalar psi2 (pow(star_i.psi4, 0.5)) ;
	psi2.std_spectral_base() ;

	const Tensor& dcov_hij = star_i.hij.derive_cov(flat) ;
	const Tensor& dcon_hij = star_i.hij.derive_con(flat) ;
	const Tensor& dcov_hij_auto = star_i.hij_auto.derive_cov(flat) ;

	const Sym_tensor& gtilde_cov = star_i.gtilde.cov() ;
	const Sym_tensor& gtilde_con = star_i.gtilde.con() ;
	const Tensor& dcov_gtilde = gtilde_cov.derive_cov(flat) ;

	// H^i and its derivatives ( = O in Dirac gauge)
	// ---------------------------------------------

	double lambda_dirac = 0. ;

	const Vector hdirac = lambda_dirac * star_i.hij.divergence(flat) ;
	const Vector hdirac_auto = lambda_dirac * 
	    star_i.hij_auto.divergence(flat) ;

	Tensor dcov_hdirac = lambda_dirac * hdirac.derive_cov(flat) ;
	dcov_hdirac.inc_dzpuis() ;
	Tensor dcov_hdirac_auto = lambda_dirac * hdirac_auto.derive_cov(flat) ;
	dcov_hdirac_auto.inc_dzpuis() ;
	Tensor dcon_hdirac_auto = lambda_dirac * hdirac_auto.derive_con(flat) ;
	dcon_hdirac_auto.inc_dzpuis() ;


	// Function exp(-(r-r_0)^2/sigma^2)
	// --------------------------------
	
	double r0 = mp.val_r(nz-2, 1, 0, 0) ;
	double sigma = 1.*r0 ;
	double om = omega ;
	
	Scalar rr (mp) ;
	rr = mp.r ;
	
	Scalar ff (mp) ;
	ff = exp( -(rr - r0)*(rr - r0)/sigma/sigma ) ;
	for (int ii=0; ii<nz-1; ii++)
	    ff.set_domain(ii) = 1. ;
	ff.set_outer_boundary(nz-1, 0) ;
	ff.std_spectral_base() ;
	

	// omdsdp
	// ---------------------------------

	Vector omdsdp (mp, CON, mp.get_bvect_cart()) ;
	Scalar yya (mp) ;
	yya = mp.ya ;
	Scalar xxa (mp) ;
	xxa = mp.xa ;
	Scalar zza (mp) ;
	zza = mp.za ;

	if (fabs(mp.get_rot_phi()) < 1e-10){ 
	    omdsdp.set(1) = - om * yya * ff ;
	    omdsdp.set(2) = om * xxa * ff ;
	    omdsdp.set(3).annule_hard() ;
	}
	else{
	    omdsdp.set(1) = om * yya * ff ;
	    omdsdp.set(2) = - om * xxa * ff ;
	    omdsdp.set(3).annule_hard() ;
	}
 
	omdsdp.set(1).set_outer_boundary(nz-1, 0) ;
	omdsdp.set(2).set_outer_boundary(nz-1, 0) ;	    
	omdsdp.std_spectral_base() ;


	// Computation of helical A^{ij}
	// ------------------------------

	const Tensor& dbeta = star_i.beta_auto.derive_con(gtilde) ;
 	Scalar div_beta = star_i.beta_auto.divergence(gtilde) ;
	
	Sym_tensor tkij_a (star_i.tkij_auto) ;
	for (int i=1; i<=3; i++) 
	    for (int j=1; j<=i; j++) {
		tkij_a.set(i, j) = dbeta(i, j) + dbeta(j, i) - 
		    double(2) /double(3) * div_beta * (gtilde.con())(i,j) ; 
	    }
	
	tkij_a = tkij_a - star_i.hij_auto.derive_lie(omdsdp) ;	
	tkij_a = 0.5 * tkij_a / nn ;   

	Sym_tensor tkij_auto_cov = tkij_a.up_down(gtilde) ;
	
	Scalar a2_auto = contract(tkij_auto_cov, 0, 1, tkij_a, 0, 1,true) ; 

	// COMP 
	const Tensor& dbeta_comp = star_i.beta_comp.derive_con(gtilde) ;
       	Scalar divbeta_comp = star_i.beta_comp.divergence(gtilde) ;
	
	Sym_tensor tkij_c (star_i.tkij_comp) ;
	for (int i=1; i<=3; i++) 
	    for (int j=i; j<=3; j++) {		
		tkij_c.set(i, j) = dbeta_comp(i, j) + dbeta_comp(j, i) - 
		    double(2) /double(3) * divbeta_comp * (gtilde.con())(i,j) ; 
	    }
	
	tkij_c = tkij_c - star_i.hij_comp.derive_lie(omdsdp) ;	
	tkij_c = 0.5 * tkij_c / nn ;
	
  	Scalar a2_comp = contract(tkij_auto_cov, 0, 1, tkij_c, 0, 1,true) ; 


//	tkij_a = star_i.tkij_auto ;
//	tkij_c = star_i.tkij_comp ;


	// Sources for N
	// ---------------

	Scalar source1(mp) ;
	Scalar source2(mp) ;
	Scalar source3(mp) ;
	Scalar source4(mp) ;
	Scalar source_tot(mp) ;

	source1 = qpig * star_i.psi4 % (star_i.ener_euler + star_i.s_euler) ; 

	source2 = star_i.psi4 % (a2_auto + a2_comp) ;

	source3 = - contract(dcov_logn_auto, 0, star_i.dcon_logn, 0, true) 
	    - 2. * contract(contract(gtilde_con, 0, star_i.dcov_phi, 0), 
			    0, dcov_logn_auto, 0, true) ;
		    
	source4 = - contract(star_i.hij, 0, 1, dcovdcov_logn_auto + 
			     dcov_logn_auto*star_i.dcov_logn, 0, 1) ;

	source_tot = source1 + source2 + source3 + source4 ;
	
	if (ll == 1) {
	lie_K_1 = nn / star_i.psi4 * (source_tot - star_i.logn_auto
					   .laplacian()) ;
	lie_K_1.dec_dzpuis(4) ;
	}
	if (ll == 2) {
	lie_K_2 = nn / star_i.psi4 * (source_tot - star_i.logn_auto
					   .laplacian()) ;
	lie_K_2.dec_dzpuis(4) ;
	}


	// Sources for hij
	// --------------

	Scalar source_tot_hij(mp) ;
	Tensor source_Sij(mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor source_Rij(mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor tens_temp(mp, 2, CON, mp.get_bvect_cart()) ;
	
	Tensor source_1 (mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor source_2 (mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor source_3a (mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor source_3b (mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor source_4 (mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor source_5 (mp, 2, CON, mp.get_bvect_cart()) ;
	Tensor source_6 (mp, 2, CON, mp.get_bvect_cart()) ;
	

	source_1 = contract(dcon_hij, 1, dcov_lnq_auto, 0) ;
	
	source_2 = - contract(dcon_hij, 2, dcov_lnq_auto, 0) 
	    - 2./3. * contract(hdirac, 0, dcov_lnq_auto, 0) * flat.con() ;
	    
	// Lie derivative of A^{ij}
	// --------------------------

	Scalar decouple_logn = (star_i.logn_auto - 1.e-8)/
	    (star_i.logn - 2.e-8) ;

	// Construction of Omega d/dphi
	// ----------------------------
	    
	// Construction of D_k \Phi^i
	Itbl type (2) ;
	type.set(0) = CON ;
	type.set(1) = COV ;

	Tensor dcov_omdsdphi (mp, 2, type, mp.get_bvect_cart()) ;
	dcov_omdsdphi.set(1,1) = 0. ;
	dcov_omdsdphi.set(2,1) = om * ff ;
	dcov_omdsdphi.set(3,1) = 0. ;
	dcov_omdsdphi.set(2,2) = 0. ;
	dcov_omdsdphi.set(3,2) = 0. ;
	dcov_omdsdphi.set(3,3) = 0. ;
	dcov_omdsdphi.set(1,2) = -om *ff ;
	dcov_omdsdphi.set(1,3) = 0. ;
	dcov_omdsdphi.set(2,3) = 0. ;
	dcov_omdsdphi.std_spectral_base() ;

	source_3a = contract(tkij_a, 0, dcov_omdsdphi, 1) ;
	source_3a.inc_dzpuis() ;

	// Source 3b
	// ------------

	source_3b = - contract(omdsdp, 0, tkij_a.derive_cov(flat), 2) ; 

 
	// Source 4
	// ---------

	source_4 = - tkij_a.derive_lie(star_i.beta) ;
	source_4.inc_dzpuis() ;
	source_4 += - 2./3. * star_i.beta.divergence(flat) * tkij_a ;

	source_5 = dcon_hdirac_auto ;
 	    
	source_6 = - 2./3. * hdirac_auto.divergence(flat) * flat.con() ;
	source_6.inc_dzpuis() ;

	// Source terms for Sij
	//---------------------
	    
	source_Sij = 8. * nn / psi4 * phi_auto.derive_con(gtilde) * 
	    contract(gtilde_con, 0, star_i.dcov_phi, 0) ;

	source_Sij += 4. / psi4 * phi_auto.derive_con(gtilde) * 
	    nn * contract(gtilde_con, 0, star_i.dcov_logn, 0) +
	    4. / psi4 * nn * contract(gtilde_con, 0, star_i.dcov_logn, 0) *
	    phi_auto.derive_con(gtilde) ;

	source_Sij += - nn / (3.*psi4) * gtilde_con * 
	    ( 0.25 * contract(gtilde_con, 0, 1, contract(dcov_hij_auto, 0, 1,
					       dcov_gtilde, 0, 1), 0, 1)
	      - 0.5 * contract(gtilde_con, 0, 1, contract(dcov_hij_auto, 0, 1,
					       dcov_gtilde, 0, 2), 0, 1)) ;

	source_Sij += - 8.*nn / (3.*psi4) * gtilde_con * 
	 contract(dcov_phi_auto, 0, contract(gtilde_con, 0, star_i.dcov_phi, 0), 0) ;

	tens_temp = nn / (3.*psi4) * hdirac.divergence(flat)*star_i.hij_auto ;
	tens_temp.inc_dzpuis() ;

	source_Sij += tens_temp ;
	   
	source_Sij += - 8./(3.*psi4) * contract(dcov_phi_auto, 0,
	   nn*contract(gtilde_con, 0, star_i.dcov_logn, 0), 0) * gtilde_con ;

	source_Sij += 2.*nn* contract(gtilde_cov, 0, 1, tkij_a *
				(tkij_a+tkij_c), 1, 3) ;

	source_Sij += - 2. * qpig * nn * ( psi4 * star_i.stress_euler 
		      - 0.33333333333333333 * star_i.s_euler * gtilde_con ) ; 
	    
	source_Sij += - 1./(psi4*psi2) * contract(gtilde_con, 1, 
			contract(gtilde_con, 1, qq*dcovdcov_lnq_auto + 
				 qq*dcov_lnq_auto*dcov_lnq, 0), 1) ;

	source_Sij += - 0.5/(psi4*psi2) * contract(contract(star_i.hij, 1,
				    dcov_hij_auto, 2), 1, dcov_qq, 0) -
	    0.5/(psi4*psi2) * contract(contract(dcov_hij_auto, 2,
				    star_i.hij, 1), 1, dcov_qq, 0) ;
					
	source_Sij += 0.5/(psi4*psi2) * contract(contract(star_i.hij, 0,
					 dcov_hij_auto, 2), 0, dcov_qq, 0) ;

	source_Sij += 1./(3.*psi4*psi2)*contract(gtilde_con, 0, 1,
       qq*dcovdcov_lnq_auto + qq*dcov_lnq_auto*dcov_lnq, 0, 1)
	    *gtilde_con ;

	source_Sij += 1./(3.*psi4*psi2) * contract(hdirac, 0, 
					    dcov_qq, 0)*star_i.hij_auto ;

	// Source terms for Rij
	//---------------------

	source_Rij = contract(star_i.hij, 0, 1, dcov_hij_auto.derive_cov(flat),
 			      2, 3) ;
	source_Rij.inc_dzpuis() ;


	source_Rij += - contract(star_i.hij_auto, 1, dcov_hdirac, 1) -
	    contract(dcov_hdirac, 1, star_i.hij_auto, 1) ;
	    
	source_Rij += contract(hdirac, 0, dcov_hij_auto, 2) ;

	source_Rij += - contract(contract(dcov_hij_auto, 1, dcov_hij, 2),
				 1, 3) ;

	source_Rij += - contract(gtilde_cov, 0, 1, contract(contract(
		gtilde_con, 0, dcov_hij_auto, 2), 0, dcov_hij, 2), 1, 3) ;

	source_Rij += contract(contract(contract(contract(gtilde_cov, 0, 
	      dcov_hij_auto, 1), 2, gtilde_con, 1), 0, dcov_hij, 1), 0, 3) +
	    contract(contract(contract(contract(gtilde_cov, 0, 
		dcov_hij_auto, 1), 0, dcov_hij, 1), 0, 3), 0, gtilde_con, 1) ;

	source_Rij += 0.5 * contract(gtilde_con*gtilde_con, 1, 3, 
		 contract(dcov_hij_auto, 0, 1, dcov_gtilde, 0, 1), 0, 1) ;

	source_Rij = source_Rij * 0.5 ;

	for(int i=1; i<=3; i++) 
	    for(int j=1; j<=i; j++) {

		source_tot_hij = source_1(i,j) + source_1(j,i) 
		    + source_2(i,j) + 2.*psi4/nn * (
			source_4(i,j) - source_Sij(i,j)) 
		    - 2.* source_Rij(i,j) +
		    source_5(i,j) + source_5(j,i) + source_6(i,j) ;
		source_tot_hij.dec_dzpuis() ;
		    
		source3 = 2.*psi4/nn * (source_3a(i,j) + source_3a(j,i)  
				 + source_3b(i,j)) ; 

		source_tot_hij = source_tot_hij + source3 ;


		if (ll == 1){
		    if(i==1 && j==1) {
			Scalar lapl (star_i.hij_auto(1,1).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_1.set(1,1) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==2 && j==1) {
			Scalar lapl (star_i.hij_auto(2,1).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_1.set(2,1) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==3 && j==1) {
			Scalar lapl (star_i.hij_auto(3,1).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_1.set(3,1) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==2 && j==2) {
			Scalar lapl (star_i.hij_auto(2,2).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_1.set(2,2) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==3 && j==2) {
			Scalar lapl (star_i.hij_auto(3,2).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_1.set(3,2) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==3 && j==3) {
			Scalar lapl (star_i.hij_auto(3,3).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_1.set(3,3) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		}

		if (ll == 2){
		    if(i==1 && j==1) {
			Scalar lapl (star_i.hij_auto(1,1).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_2.set(1,1) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==2 && j==1) {
			Scalar lapl (star_i.hij_auto(2,1).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_2.set(2,1) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==3 && j==1) {
			Scalar lapl (star_i.hij_auto(3,1).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_2.set(3,1) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==2 && j==2) {
			Scalar lapl (star_i.hij_auto(2,2).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_2.set(2,2) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==3 && j==2) {
			Scalar lapl (star_i.hij_auto(3,2).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_2.set(3,2) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		    if(i==3 && j==3) {
			Scalar lapl (star_i.hij_auto(3,3).laplacian()) ;
			lapl.dec_dzpuis() ;
			lie_aij_2.set(3,3) = nn / (2.*psi4) * 
			    (lapl-source_tot_hij) ;
		    }
		}
	    }
    }

    lie_aij_1.dec_dzpuis(3) ;
    lie_aij_2.dec_dzpuis(3) ;
		
    int nz = star1.mp.get_mg()->get_nzone() ;
  
    // Construction of an auxiliar grid and mapping (Last domain is at lambda)
    double* bornes = new double [6] ;
    bornes[nz] = __infinity ;
    bornes[4] = M_PI / omega ;
    bornes[3] = M_PI / omega * 0.5 ;
    bornes[2] = M_PI / omega  * 0.2 ;
    bornes[1] = M_PI / omega  * 0.1 ;
    bornes[0] = 0 ;
    
    Map_af mapping (*(star1.mp.get_mg()), bornes) ;
    
    delete [] bornes ; 
 
    
    Sym_tensor lie_aij2_1 (star1.mp, CON, star1.mp.get_bvect_cart()) ;
    Sym_tensor lie_aij_tot_1 (star1.mp, CON, star1.mp.get_bvect_cart()) ;
    Sym_tensor lie_aij_tot (mapping, CON, mapping.get_bvect_cart()) ;

    Scalar lie_K2_1 (star1.mp) ;
    lie_K2_1.set_etat_qcq() ;
    Scalar lie_K_tot_1 (star1.mp) ;


    // Importation on the mapping 1
    // -------------------------------

    lie_K2_1.import(lie_K_2) ;
    lie_K_tot_1 = lie_K_1 + lie_K2_1 ;
    lie_K_tot_1.inc_dzpuis(2) ;

    lie_aij_2.change_triad(star1.mp.get_bvect_cart()) ;    
    for(int i=1; i<=3; i++) 
	for(int j=1; j<=i; j++) {
	    lie_aij2_1.set(i,j).import(lie_aij_2(i,j)) ;
	    lie_aij2_1.set(i,j).set_spectral_va().set_base(lie_aij_2(i,j).
					    get_spectral_va().get_base()) ;
	}

    lie_aij_tot_1 = lie_aij_1 + lie_aij2_1 ;
    lie_aij_tot_1.inc_dzpuis(2) ;

    
    Sym_tensor lie_kij_tot (star1.mp, CON, star1.mp.get_bvect_cart()) ;
    lie_kij_tot = lie_aij_tot_1/star1.psi4 + 1./3.*star1.gamma.con()*
	lie_K_tot_1 ;


    cout << " IN THE CENTER OF STAR 1 " << endl 
	 << " ----------------------- " << endl ;
    /*
    cout << " components xx, xy, yy, xz, yz, zz" << endl ;
    for(int i=1; i<=3; i++) 
	for(int j=1; j<=i; j++) {
	    Scalar resu(lie_kij_tot(i,j)*lie_kij_tot(i,j)) ;
	    cout << "i = " << i << ", j = " << j << endl ; 
	    cout << "norme de la diff " << endl 
		 << norme(resu)/(nr*nt*np) << endl ;

	    // Computation of the integral
	    // -----------------------------

	    Tbl integral (nz) ;
	    integral.annule_hard()  ;
	    Tbl integ (resu.integrale_domains()) ;
	    for (int mm=0; mm<nz; mm++) 
		for (int pp=0; pp<=mm; pp++) 
		    integral.set(mm) += integ(pp) ; 
	    cout << sqrt(integral) / sqrt(omega) << endl ;  // To get 
	    // dimensionless quantity
	}
    */

    cout << "L2 norm of L_k K^{ab} " << endl ;
    Scalar determinant (pow(star1.get_gamma().determinant(), 0.5)) ;
    determinant.std_spectral_base() ;
    Scalar resu(2.*contract(lie_kij_tot, 0, 1, 
			 lie_kij_tot.up_down(star1.gamma), 0, 1)
		*determinant) ;
    Tbl integral (nz) ;
    integral.annule_hard() ;
    Tbl integ (resu.integrale_domains()) ;
    for (int mm=0; mm<nz; mm++) 
	for (int pp=0; pp<=mm; pp++) 
	    integral.set(mm) += integ(pp) ; 
    cout << sqrt(integral) * sqrt(mass_adm()*ggrav) << endl ;  
    cout << sqrt(integral) / sqrt(omega) << endl ;  // To get 
    // dimensionless quantity

    
    cout << "omega = " << omega << endl ;
    cout << "mass_adm = " << mass_adm() << endl ;
    

    lie_kij_tot.dec_dzpuis(2) ;
    
    cout << "Position du centre de l'etoile x/lambda = "
	 << -star1.get_mp().get_ori_x() * omega / M_PI << " in Lorene units"
	 << endl ;
   


    // Importation on the mapping defined in the center of mass
    // -------------------------------------------------------------
/*
    lie_aij_tot_1.change_triad(mapping.get_bvect_cart()) ;    
    for(int i=1; i<=3; i++) 
	for(int j=1; j<=i; j++) {
	    lie_aij_tot.set(i,j).import(lie_aij_tot_1(i,j)) ;
	    lie_aij_tot.set(i,j).set_spectral_va().set_base(lie_aij_tot_1(i,j).
					    get_spectral_va().get_base()) ;
	}

    lie_aij_tot.inc_dzpuis(2) ;

    cout << " IN THE CENTER OF MASS : " << endl 
	 << " ----------------------- " << endl ;

    cout << " components xx, xy, yy, xz, yz, zz" << endl ;
    for(int i=1; i<=3; i++) 
	for(int j=1; j<=i; j++) {
	    Scalar resu(lie_aij_tot(i,j)*lie_aij_tot(i,j)) ;
	    cout << "i = " << i << ", j = " << j << endl ; 
	    cout << "norme de la diff " << endl 
		 << norme(resu)/(nr*nt*np) << endl ;

	    Tbl integral (nz) ;
	    integral.annule_hard()  ;
	    Tbl integ (resu.integrale_domains()) ;
	    for (int mm=0; mm<nz; mm++) 
		for (int pp=0; pp<=mm; pp++) 
		    integral.set(mm) += integ(pp) ; 
	    cout << sqrt(integral) / sqrt(omega) << endl ;  // To get 
	    // dimensionless quantity
	}

    cout << "L2 norm of L_k K^{ab} " << endl ;
    resu = contract(lie_aij_tot, 0, 1, 
		    lie_aij_tot.up_down(star1.gtilde), 0, 1) ;
    integral.annule_hard()  ;
    integ = resu.integrale_domains() ;
    for (int mm=0; mm<nz; mm++) 
	for (int pp=0; pp<=mm; pp++) 
	    integral.set(mm) += integ(pp) ; 
    cout << sqrt(integral) / sqrt(omega) << endl ;  // To get 
    // dimensionless quantity
    */

    cout << "Omega M = " << omega * mass_adm()*ggrav << endl ;

}
}
