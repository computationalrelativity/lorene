/*
 *  Method for vector Poisson equation inverting eqs. for V^r and eta as a block.
 *  Space considered is R3 minus a ball.
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2007  Jerome Novak
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

char vector_poisson_boundary2_C[] = "$Header$" ;

/* 
 *
 * * $Header$
 *
 */

// C headers
#include <assert.h>
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "metric.h"
#include "diff.h"
#include "param_elliptic.h"
#include "proto.h"
#include "graphique.h"
#include "map.h"
#include "utilitaires.h"



void Vector::poisson_boundary2(double lam, Vector& resu, Scalar boundvr, Scalar boundeta, Scalar boundmu, double dir_vr, double neum_vr, double dir_eta, double neum_eta, double dir_mu, double neum_mu) const {

    const Map_af* mpaff = dynamic_cast<const Map_af*>(mp) ;
#ifndef NDEBUG 
    for (int i=0; i<3; i++)
	assert(cmp[i]->check_dzpuis(4)) ;
    // All this has a meaning only for spherical components:
    const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
    assert(bvs != 0x0) ; 
    //## ... and affine mapping, for the moment!
    assert( mpaff != 0x0) ;
#endif

     if (fabs(lam + 1.) < 1.e-8) {

       cout <<" lambda = -1 is not supported by this code!" << endl;
       abort(); 
     }



     // Extraction of Mtbl_cf objects from boundary informations; 
     Scalar bound_vr2 = boundvr; 
     bound_vr2.set_spectral_va().ylm(); 
     Scalar bound_eta2 = boundeta; 
     bound_eta2.set_spectral_va().ylm();
     Scalar bound_mu2 = boundmu; bound_mu2.set_spectral_va().ylm();
     // Necessaire pour poisoon_dirichlet:
     Valeur vbound_mu2 = bound_mu2.get_spectral_va();


     
     const  Mtbl_cf *bound_vr = bound_vr2.get_spectral_va().c_cf; 
     const  Mtbl_cf *bound_mu = bound_mu2.get_spectral_va().c_cf;
     const Mtbl_cf *bound_eta = bound_eta2.get_spectral_va().c_cf;




    // Some working objects
    //---------------------
  const Mg3d& mg = *mpaff->get_mg() ;
  int nz = mg.get_nzone() ; int nzm1 = nz - 1;
  assert(mg.get_type_r(0) == RARE) ;
  assert(mg.get_type_r(nzm1) == UNSURR) ;
  Scalar S_r = *cmp[0] ;
  Scalar S_eta = eta() ;
  Scalar het(*mpaff) ; 
  Scalar vr(*mpaff) ; 
  bool all_zero = false ; // Artefact from R3 case; to remove eventually
  if (S_r.get_etat() == ETATZERO) {
      if (S_eta.get_etat() == ETATZERO) {

	  S_r.annule_hard() ;
	  S_r.set_spectral_base(S_eta.get_spectral_base()) ;
	  S_eta.annule_hard() ;
	  S_eta.set_spectral_base(S_r.get_spectral_base()) ;
	  // A little bit stupid for the moment; 
	  // But calculation HAS to be made even with zero-source!!

// 	  vr.set_etat_zero() ;
// 	  het.set_etat_zero() ;
// 	  all_zero = true ;
      }
      else {
	  S_r.annule_hard() ;
	  S_r.set_spectral_base(S_eta.get_spectral_base()) ;
      }
  }
  if ((S_eta.get_etat() == ETATZERO)&&(!all_zero)) {
      S_eta.annule_hard() ;
      S_eta.set_spectral_base(S_r.get_spectral_base()) ;
  }
  if (!all_zero) {
  S_r.set_spectral_va().ylm() ;
  S_eta.set_spectral_va().ylm() ;
  const Base_val& base = S_eta.get_spectral_va().base ;
  Mtbl_cf sol_part_eta(mg, base) ; sol_part_eta.annule_hard() ;
  Mtbl_cf sol_part_vr(mg, base) ; sol_part_vr.annule_hard() ;
  Mtbl_cf sol_hom_un_eta(mg, base) ; sol_hom_un_eta.annule_hard() ;
  Mtbl_cf sol_hom_deux_eta(mg, base) ; sol_hom_deux_eta.annule_hard() ;
  Mtbl_cf sol_hom_trois_eta(mg, base) ; sol_hom_trois_eta.annule_hard() ;
  Mtbl_cf sol_hom_quatre_eta(mg, base) ; sol_hom_quatre_eta.annule_hard() ;
  Mtbl_cf sol_hom_un_vr(mg, base) ; sol_hom_un_vr.annule_hard() ;
  Mtbl_cf sol_hom_deux_vr(mg, base) ; sol_hom_deux_vr.annule_hard() ;
  Mtbl_cf sol_hom_trois_vr(mg, base) ; sol_hom_trois_vr.annule_hard() ;
  Mtbl_cf sol_hom_quatre_vr(mg, base) ; sol_hom_quatre_vr.annule_hard() ;
 
  // Solution of the l=0 part for Vr
  //---------------------------------
  Scalar sou_l0 = (*cmp[0]) / (1. + lam) ;
  sou_l0.set_spectral_va().ylm() ;

 Param_elliptic param_l0(sou_l0) ;
  for (int l=0; l<nz; l++){
      param_l0.set_poisson_vect_r(l, true) ;
  }


   Scalar vrl0 = sou_l0.sol_elliptic_boundary(param_l0, *bound_vr, 1., 0.) ;



  // Build-up & inversion of the system for (eta, V^r) in each domain (except for the nucleus!)
  //-----------------------------------------------------------------

  // Nucleus
  //--------
  int nr = mg.get_nr(1) ;
  int nt = mg.get_nt(1) ;
  int np = mg.get_np(1) ;
  double alpha = mpaff->get_alpha()[1] ; double alp2 = alpha*alpha ;
  double beta = mpaff->get_beta()[1] ;

 
  int l_q = 0 ; int m_q = 0; int base_r = 0 ;

//   // Loop on l and m
//   //----------------
//   for (int k=0 ; k<np+1 ; k++) {
//       for (int j=0 ; j<nt ; j++) { 
// 	  base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
// 	  if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 0) ) {
// 	      int aa = 0 ; int bb = 0 ;  int nr0 = 0 ;
// 	      if (base_r == R_CHEBP) {
// 		  nr0 = nr - 1 ;
// 		  aa = 0 ; bb = 1 ;
// 	      }
// 	      else {
// 		  assert (base_r == R_CHEBI) ;
// 		  nr0 = nr - 2 ;
// 		  aa = 2 ; bb = 1 ;
// 	      }
// 	      int d0 = nr - nr0 ;
// 	      int nrtot = 2*nr0 ;
// 	      Matrice oper(2*nr, 2*nr) ; oper.set_etat_qcq() ;
// 	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
// 	      Diff_dsdx2 d2(base_r, nr) ; const Matrice& md2 = d2.get_matrice() ;
// 	      Diff_sxdsdx xd(base_r, nr) ; const Matrice& mxd = xd.get_matrice() ;
// 	      Diff_sx2 s2(base_r, nr) ; const Matrice& ms2 = s2.get_matrice() ;

// 	      // Building the operator
// 	      //----------------------
// 	      for (int lin=0; lin<nr; lin++) { //eq.1 
// 		  for (int col=0; col<nr; col++) 
// 		      oper.set(lin,col) = 
// 			  (md2(lin,col) + 2*mxd(lin,col) 
// 			   -(lam+1)*l_q*(l_q+1)*ms2(lin,col)) / alp2 ;
// 		  for (int col=0; col<nr; col++) 
// 		      oper.set(lin,col+nr) = 
// 			  (lam*mxd(lin,col) + 2*(1+lam)*ms2(lin,col)) / alp2 ;
// 	      }
// 	      for (int lin=0; lin<nr; lin++) { //eq.2
// 		  for (int col=0; col<nr; col++)
// 		      oper.set(lin+nr,col) = 
// 			  (-lam*l_q*(l_q+1)*mxd(lin,col) 
// 			   +(lam+2)*l_q*(l_q+1)*ms2(lin,col)) / alp2 ;
// 		  for (int col=0; col<nr; col++)
// 		      oper.set(lin+nr, col+nr) = 
// 			  ((lam+1)*(md2(lin,col) + 2*mxd(lin,col)) 
// 			   - (2*(lam+1) + l_q*(l_q+1))*ms2(lin,col)) / alp2 ;
// 	      }
// 	      bool pb_eta = ( ( fabs( lam*double(l_q+3) + 2 ) < 0.01) && (l_q <=2) ) ;
// 	      if (!pb_eta) {
// 		  for (int col=0; col<nr; col++) 
// 		      oper.set(nr0-1, col) = 1 ;
// 		  for (int col=0; col<nr; col++) 
// 		      oper.set(nr0-1,nr+col) = 0 ;		  
// 	      }
// 	      if ((l_q > 2)||pb_eta) {
// 		  for (int col=0; col<nr; col++) 
// 		      oper.set(nr+nr0-1, col) = 0 ;
// 		  for (int col=0; col<nr; col++) 
// 		      oper.set(nr+nr0-1,nr+col) = 1 ;		  
// 	      }

// 	      Matrice op2(nrtot, nrtot) ;
// 	      op2.set_etat_qcq() ;
// 	      for (int i=0; i<nr0; i++) {
// 		  for (int col=0; col<nr0;col++) 
// 		      op2.set(i,col) = (aa*col+bb)*oper(i,col+1) 
// 			  + (aa*(col+1)+bb)*oper(i,col) ;
// 		  for (int col=0; col<nr0;col++) 
// 		      op2.set(i,col+nr0) = (aa*col+bb)*oper(i,col+nr+1)
// 			  + (aa*(col+1)+bb)*oper(i,col+nr) ;
// 	      }
	      
// 	      for (int i=nr0; i<nrtot; i++) {
// 		  for (int col=0; col<nr0;col++) 
// 		      op2.set(i,col) = (aa*col+bb)*oper(i+d0,col+1) 
// 			  + (aa*(col+1)+bb)*oper(i+d0,col) ;
// 		  for (int col=0; col<nr0;col++) 
// 		      op2.set(i,col+nr0) = (aa*col+bb)*oper(i+d0,col+nr+1)
// 			  + (aa*(col+1)+bb)*oper(i+d0,col+nr) ;
// 	      }
// 	      op2.set_lu() ;

//  	      // Filling the r.h.s
//  	      //------------------
//  	      for (int i=0; i<nr0; i++)  //eq.1
//  		  sec_membre.set(i) = (*S_eta.get_spectral_va().c_cf)(0, k, j, i) ;
//  	      if (!pb_eta) sec_membre.set(nr0-1) = 0 ;
//  	      for (int i=0; i<nr0; i++) //eq.2
//  		  sec_membre.set(i+nr0) 
//  		      = (*S_r.get_spectral_va().c_cf)(0, k, j, i) ;
//  	      if ((l_q > 2)||pb_eta) sec_membre.set(nrtot-1) = 0 ;

//  	      // Inversion of the "big" operator
//  	      //--------------------------------
//  	      Tbl big_res = op2.inverse(sec_membre) ;
	      
//  	      // Putting coefficients of eta and Vr to individual arrays
//  	      //--------------------------------------------------------
//  	      sol_part_eta.set(0, k, j, 0) = (aa+bb)*big_res(0) ;
// 	      for (int i=1; i<nr0; i++)
// 		  sol_part_eta.set(0, k, j, i) = (aa*(i-1)+bb)*big_res(i-1) 
// 		      + (aa*(i+1)+bb)*big_res(i);
// 	      sol_part_eta.set(0, k, j, nr0) = big_res(nr0-1)*(aa*(nr0-1) + bb) ;
// 	      sol_part_vr.set(0, k, j, 0) = (aa+bb)*big_res(nr0) ;
// 	      for (int i=1; i<nr0; i++)
// 		  sol_part_vr.set(0, k, j, i) = (aa*(i-1)+bb)*big_res(nr0+i-1) 
// 		      + (aa*(i+1)+bb)*big_res(nr0+i);
// 	      sol_part_vr.set(0, k, j, nr0) = big_res(nrtot-1)*(aa*(nr0-1)+bb) ;
// 	      if (base_r == R_CHEBI) {
// 		  sol_part_eta.set(0, k, j, nr-1) = 0 ;
// 		  sol_part_vr.set(0, k, j, nr-1) = 0 ;
// 	      }

//  	      // Homogeneous solutions (only r^(l-1) and r^(l+1) in the nucleus)
// 	      if (l_q<=2) {
// 		  double fac_eta = lam*double(l_q+3) + 2. ;
// 		  double fac_vr = double(l_q+1)*(lam*l_q - 2.) ;
// 		  Tbl sol_hom1 = solh(nr, l_q-1, 0., base_r) ;
// 		  Tbl sol_hom2 = solh(nr, l_q+1, 0., base_r) ;
// 		  for (int i=0 ; i<nr ; i++) {
// 		      sol_hom_un_eta.set(0, k, j, i) = sol_hom1(i) ;
// 		      sol_hom_un_vr.set(0, k, j, i) = l_q*sol_hom1(i) ;
// 		      sol_hom_trois_eta.set(0, k, j, i) = fac_eta*sol_hom2(i) ; 
// 		      sol_hom_trois_vr.set(0, k, j, i) = fac_vr*sol_hom2(i) ; 
// 		  }
// 	      }
// 	      else {
// 		  for (int i=0; i<nrtot; i++) sec_membre.set(i) = 0 ;
// 		  // First homogeneous sol.
// 		  sec_membre.set(nr0-1) = 1 ;
// 		  big_res = op2.inverse(sec_membre) ;

// 		  sol_hom_un_eta.set(0, k, j, 0) = (aa+bb)*big_res(0) ;
// 		  for (int i=1; i<nr0; i++)
// 		      sol_hom_un_eta.set(0, k, j, i) = (aa*(i-1)+bb)*big_res(i-1) 
// 			  + (aa*(i+1)+bb)*big_res(i);
// 		  sol_hom_un_eta.set(0, k, j, nr0) = big_res(nr0-1)*(aa*(nr0-1) + bb) ;
// 		  sol_hom_un_vr.set(0, k, j, 0) = (aa+bb)*big_res(nr0) ;
// 		  for (int i=1; i<nr0; i++)
// 		      sol_hom_un_vr.set(0, k, j, i) = (aa*(i-1)+bb)*big_res(nr0+i-1) 
// 			  + (aa*(i+1)+bb)*big_res(nr0+i);
// 		  sol_hom_un_vr.set(0, k, j, nr0) = big_res(nrtot-1)*(aa*(nr0-1)+bb) ;
// 		  if (base_r == R_CHEBI) {
// 		      sol_hom_un_eta.set(0, k, j, nr-1) = 0 ;
// 		      sol_hom_un_vr.set(0, k, j, nr-1) = 0 ;
// 		  }

// 		  sec_membre.set(nr0-1) = 0 ;
// 		  sec_membre.set(nrtot-1) = 1 ;
// 		  big_res = op2.inverse(sec_membre) ;

// 		  sol_hom_trois_eta.set(0, k, j, 0) = (aa+bb)*big_res(0) ;
// 		  for (int i=1; i<nr0; i++)
// 		      sol_hom_trois_eta.set(0, k, j, i) = (aa*(i-1)+bb)*big_res(i-1) 
// 			  + (aa*(i+1)+bb)*big_res(i);
// 		  sol_hom_trois_eta.set(0, k, j, nr0) = big_res(nr0-1)*(aa*(nr0-1)+bb) ;
// 		  sol_hom_trois_vr.set(0, k, j, 0) = (aa+bb)*big_res(nr0) ;
// 		  for (int i=1; i<nr0; i++)
// 		      sol_hom_trois_vr.set(0, k, j, i) = (aa*(i-1)+bb)*big_res(nr0+i-1) 
// 			  + (aa*(i+1)+bb)*big_res(nr0+i);
// 		  sol_hom_trois_vr.set(0, k, j, nr0) = big_res(nrtot-1)*(aa*(nr0-1)+bb) ;
// 		  if (base_r == R_CHEBI) {
// 		      sol_hom_trois_eta.set(0, k, j, nr-1) = 0 ;
// 		      sol_hom_trois_vr.set(0, k, j, nr-1) = 0 ;
// 		  }
		  
// 	      }
// 	  }
//       }
//  }    


  // Shells
  //-------
  for (int zone=1 ; zone<nzm1 ; zone++) {
      nr = mg.get_nr(zone) ; 
      assert (nr > 5) ;
      assert(nt == mg.get_nt(zone)) ;
      assert(np == mg.get_np(zone)) ;
      alpha = mpaff->get_alpha()[zone] ;
      beta = mpaff->get_beta()[zone] ;
      double ech = beta / alpha ;

  // Loop on l and m
  //----------------
      for (int k=0 ; k<np+1 ; k++) {
      for (int j=0 ; j<nt ; j++) {
	  base.give_quant_numbers(zone, k, j, m_q, l_q, base_r) ;
	  if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q != 0) ) {
	      int dege1 = 2 ; //degeneracy of eq.1
	      int dege2 = 2 ;  //degeneracy of eq.2
	      int nr_eq1 = nr - dege1 ; //Eq.1 is for eta
	      int nr_eq2 = nr - dege2 ; //Eq.2 is for V^r
	      int nrtot = 2*nr ;
	      Matrice oper(nrtot, nrtot) ; oper.set_etat_qcq() ;
	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
	      Diff_x2dsdx2 x2d2(base_r, nr); const Matrice& m2d2 = x2d2.get_matrice() ;
	      Diff_xdsdx2 xd2(base_r, nr) ; const Matrice& mxd2 = xd2.get_matrice() ;
	      Diff_dsdx2 d2(base_r, nr) ; const Matrice& md2 = d2.get_matrice() ;
	      Diff_xdsdx xd(base_r, nr) ; const Matrice& mxd = xd.get_matrice() ;
	      Diff_dsdx d1(base_r, nr) ; const Matrice& md = d1.get_matrice() ;
	      Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;

	      // Building the operator
	      //----------------------
	      for (int lin=0; lin<nr_eq1; lin++) { 
		  for (int col=0; col<nr; col++) 
		      oper.set(lin,col) 
			  = m2d2(lin,col) + 2*ech*mxd2(lin,col) + ech*ech*md2(lin,col) 
			  + 2*(mxd(lin,col) + ech*md(lin,col)) 
			  - (lam+1)*l_q*(l_q+1)*mid(lin,col) ;
		  for (int col=0; col<nr; col++) 
		      oper.set(lin,col+nr) 
			  = lam*(mxd(lin,col) + ech*md(lin,col)) 
			  + 2*(1+lam)*mid(lin,col) ; 
	      }
	      oper.set(nr_eq1, 0) = 1 ;
	      for (int col=1; col<2*nr; col++) 
		  oper.set(nr_eq1, col) = 0 ;
	      oper.set(nr_eq1+1, 0) = 0 ;
	      oper.set(nr_eq1+1, 1) = 1 ;
	      for (int col=2; col<2*nr; col++) 
		  oper.set(nr_eq1+1, col) = 0 ;
	      for (int lin=0; lin<nr_eq2; lin++) {
		  for (int col=0; col<nr; col++)
		      oper.set(lin+nr,col) 
			  = -l_q*(l_q+1)*(lam*(mxd(lin,col) + ech*md(lin,col))
					  - (lam+2)*mid(lin,col)) ;
		  for (int col=0; col<nr; col++)
		      oper.set(lin+nr, col+nr) 
			  = (lam+1)*(m2d2(lin,col) + 2*ech*mxd2(lin,col) 
				     + ech*ech*md2(lin,col) 
				     + 2*(mxd(lin,col) + ech*md(lin,col)))
			  -(2*(lam+1)+l_q*(l_q+1))*mid(lin,col) ;		  
	      }
	      for (int col=0; col<nr; col++) 
		  oper.set(nr+nr_eq2, col) = 0 ;
	      oper.set(nr+nr_eq2, nr) = 1 ;
	      for (int col=nr+1; col<2*nr; col++) 
		  oper.set(nr+nr_eq2, col) = 0 ;
	      for (int col=0; col<nr+1; col++) 
		  oper.set(nr+nr_eq2+1, col) = 0 ;
	      oper.set(nr+nr_eq2+1, nr+1) = 1 ;
	      for (int col=nr+2; col<2*nr; col++) 
		  oper.set(nr+nr_eq2+1, col) = 0 ;

	      oper.set_lu() ;
	      
	      // Filling the r.h.s
	      //------------------
	      Tbl sr(nr) ; sr.set_etat_qcq() ; 
	      Tbl seta(nr) ; seta.set_etat_qcq() ;
	      for (int i=0; i<nr; i++) {
		  sr.set(i) = (*S_r.get_spectral_va().c_cf)(zone, k, j, i);
		  seta.set(i) = (*S_eta.get_spectral_va().c_cf)(zone, k, j, i) ;
	      }
	      Tbl xsr= sr ;  Tbl x2sr= sr ;
	      Tbl xseta= seta ; Tbl x2seta = seta ;
	      multx2_1d(nr, &x2sr.t, base_r) ; multx_1d(nr, &xsr.t, base_r) ;
	      multx2_1d(nr, &x2seta.t, base_r) ; multx_1d(nr, &xseta.t, base_r) ;
	      
	      for (int i=0; i<nr_eq1; i++) 
		  sec_membre.set(i) = alpha*(alpha*x2seta(i) + 2*beta*xseta(i))
		      + beta*beta*seta(i);
	      sec_membre.set(nr_eq1) = 0 ;
	      sec_membre.set(nr_eq1+1) = 0 ;
	      for (int i=0; i<nr_eq2; i++)
		  sec_membre.set(i+nr) = beta*beta*sr(i) 
		      + alpha*(alpha*x2sr(i) + 2*beta*xsr(i)) ;
	      sec_membre.set(nr+nr_eq2) = 0 ;
	      sec_membre.set(nr+nr_eq2+1) = 0 ;

	      // Inversion of the "big" operator
	      //--------------------------------
	      Tbl big_res = oper.inverse(sec_membre) ;
	      for (int i=0; i<nr; i++) {
		  sol_part_eta.set(zone, k, j, i) = big_res(i) ;
		  sol_part_vr.set(zone, k, j, i) = big_res(nr+i) ;		  
	      }
		  
	      // Getting the homogeneous solutions
	      //----------------------------------
	      sec_membre.annule_hard() ;
	      sec_membre.set(nr_eq1) = 1 ;
	      big_res = oper.inverse(sec_membre) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_hom_un_eta.set(zone, k, j, i) = big_res(i) ;
		  sol_hom_un_vr.set(zone, k, j, i) = big_res(nr+i) ;
	      }
	      sec_membre.set(nr_eq1) = 0 ;
	      sec_membre.set(nr_eq1+1) = 1 ;
	      big_res = oper.inverse(sec_membre) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_hom_deux_eta.set(zone, k, j, i) = big_res(i) ;
		  sol_hom_deux_vr.set(zone, k, j, i) = big_res(nr+i) ;
	      }
	      sec_membre.set(nr_eq1+1) = 0 ;
	      sec_membre.set(nr+nr_eq2) = 1 ;
	      big_res = oper.inverse(sec_membre) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_hom_trois_eta.set(zone, k, j, i) = big_res(i) ;
		  sol_hom_trois_vr.set(zone, k, j, i) = big_res(nr+i) ;
	      }
	      sec_membre.set(nr+nr_eq2) = 0 ;
	      sec_membre.set(nr+nr_eq2+1) = 1 ;
	      big_res = oper.inverse(sec_membre) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_hom_quatre_eta.set(zone, k, j, i) = big_res(i) ;
		  sol_hom_quatre_vr.set(zone, k, j, i) = big_res(nr+i) ;
	      }

	  } 
      }
      }
  }

  // Compactified external domain
  //-----------------------------
  nr = mg.get_nr(nzm1) ;
  assert(nt == mg.get_nt(nzm1)) ;
  assert(np == mg.get_np(nzm1)) ;
  alpha = mpaff->get_alpha()[nzm1] ; alp2 = alpha*alpha ;
  assert (nr > 4) ;
  
  // Loop on l and m
  //----------------
  for (int k=0 ; k<np+1 ; k++) {
      for (int j=0 ; j<nt ; j++) { 
	  base.give_quant_numbers(nzm1, k, j, m_q, l_q, base_r) ;
	  if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q != 0) ) {
	      int dege1 = 3; //degeneracy of eq.1
	      int dege2 = (l_q == 1 ? 2 : 3); //degeneracy of eq.2
	      int nr_eq1 = nr - dege1 ; //Eq.1 is for eta
	      int nr_eq2 = nr - dege2 ; //Eq.2 is the div-free condition
	      int nrtot = 2*nr ;
	      Matrice oper(nrtot, nrtot) ; oper.set_etat_qcq() ;
	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
	      Diff_dsdx2 d2(base_r, nr) ; const Matrice& md2 = d2.get_matrice() ;
	      Diff_sxdsdx sxd(base_r, nr) ; const Matrice& mxd = sxd.get_matrice() ;
	      Diff_sx2 sx2(base_r, nr) ; const Matrice& ms2 = sx2.get_matrice() ;

	      // Building the operator
	      //----------------------
	      for (int lin=0; lin<nr_eq1; lin++) { 
		  for (int col=0; col<nr; col++) 
		      oper.set(lin,col) 
			  = (md2(lin,col) - (lam+1)*l_q*(l_q+1)*ms2(lin,col))/alp2 ;
			  for (int col=0; col<nr; col++) 
			      oper.set(lin,col+nr) = 
				  (-lam*mxd(lin,col) + 2*(1+lam)*ms2(lin,col)) / alp2 ;
	      }
	      for (int col=0; col<nr; col++)
		  oper.set(nr_eq1, col) = 1 ;
	      for (int col=nr; col<nrtot; col++)
		  oper.set(nr_eq1, col) = 0 ;
	      int pari = -1 ;
	      for (int col=0; col<nr; col++) {
		  oper.set(nr_eq1+1, col) = pari*col*col ;
		  pari = pari ;
	      }
	      for (int col=nr; col<nrtot; col++)
		  oper.set(nr_eq1+1, col) = 0 ;
	      oper.set(nr_eq1+2, 0) = 1 ;
	      for (int col=1; col<nrtot; col++)
		  oper.set(nr_eq1+2, col) = 0 ;
	      for (int lin=0; lin<nr_eq2; lin++) {
		  for (int col=0; col<nr; col++)
		      oper.set(lin+nr,col) = (l_q*(l_q+1)*(lam*mxd(lin,col) 
			       + (lam+2)*ms2(lin,col))) / alp2 ;
		  for (int col=0; col<nr; col++)
		      oper.set(lin+nr, col+nr) = ((lam+1)*md2(lin,col) 
			     - (2*(lam+1) + l_q*(l_q+1))*ms2(lin,col)) / alp2 ;
	      }
	      for (int col=0; col<nr; col++)
		  oper.set(nr+nr_eq2, col) = 0 ;
	      for (int col=nr; col<nrtot; col++)
		  oper.set(nr+nr_eq2, col) = 1 ;
	      for (int col=0; col<nr; col++)
		  oper.set(nr+nr_eq2+1, col) = 0 ;
	      pari = -1 ;
	      for (int col=0; col<nr; col++) {
		  oper.set(nr+nr_eq2+1, nr+col) = pari*col*col ;
		  pari = pari ;
	      }
	      if (l_q > 1) {
		  for (int col=0; col<nr; col++)
		      oper.set(nr+nr_eq2+2, col) = 0 ;
		  oper.set(nr+nr_eq2+2, nr) = 1 ;
		  for (int col=nr+1; col<nrtot; col++)
		      oper.set(nr+nr_eq2+2, col) = 0 ;
	      }
	      oper.set_lu() ;

	      // Filling the r.h.s
	      //------------------
	      for (int i=0; i<nr_eq1; i++) 
		  sec_membre.set(i) = (*S_eta.get_spectral_va().c_cf)(nzm1, k, j, i) ;
	      for (int i=nr_eq1; i<nr; i++)
		  sec_membre.set(i) = 0 ;
	      for (int i=0; i<nr_eq2; i++)
		  sec_membre.set(i+nr) =(*S_r.get_spectral_va().c_cf)(nzm1, k, j, i);
	      for (int i=nr_eq2; i<nr; i++)
		  sec_membre.set(nr+i) = 0 ;
	      Tbl big_res = oper.inverse(sec_membre) ;
		  
	      // Putting coefficients of h and v to individual arrays
	      //-----------------------------------------------------
	      for (int i=0; i<nr; i++) {
		  sol_part_eta.set(nzm1, k, j, i) = big_res(i) ;
		  sol_part_vr.set(nzm1, k, j, i) = big_res(i+nr) ;
	      }

	      // Homogeneous solutions (only 1/r^(l+2) and 1/r^l in the CED)
	      //------------------------------------------------------------
	      if (l_q == 1) {		  
		  Tbl sol_hom1 = solh(nr, 0, 0., base_r) ;
		  Tbl sol_hom2 = solh(nr, 2, 0., base_r) ;
		  double fac_eta = lam + 2. ;
		  double fac_vr = 2*lam + 2. ;
		  for (int i=0 ; i<nr ; i++) {
		      sol_hom_deux_eta.set(nzm1, k, j, i) = sol_hom2(i) ;
		      sol_hom_quatre_eta.set(nzm1, k, j, i) = fac_eta*sol_hom1(i) ;
		      sol_hom_deux_vr.set(nzm1, k, j, i) = -2*sol_hom2(i) ;
		      sol_hom_quatre_vr.set(nzm1, k, j, i) = fac_vr*sol_hom1(i) ;
		  }
	      }
	      else {
		  sec_membre.annule_hard() ;
		  sec_membre.set(nr-1) = 1 ;
		  big_res = oper.inverse(sec_membre) ;

		  for (int i=0; i<nr; i++) {
		      sol_hom_deux_eta.set(nzm1, k, j, i) = big_res(i) ;
		      sol_hom_deux_vr.set(nzm1, k, j, i) = big_res(nr+i) ;
		  }
		  sec_membre.set(nr-1) = 0 ;
		  sec_membre.set(2*nr-1) = 1 ;
		  big_res = oper.inverse(sec_membre) ;
		  for (int i=0; i<nr; i++) {
		      sol_hom_quatre_eta.set(nzm1, k, j, i) = big_res(i) ;
		      sol_hom_quatre_vr.set(nzm1, k, j, i) = big_res(nr+i) ;
		  }
	      }
	  } 
      }
  }

  // Now let's match everything ...
  //-------------------------------

  // Resulting V^r & eta
  vr.set_etat_qcq() ;
  vr.set_spectral_base(base) ;
  vr.set_spectral_va().set_etat_cf_qcq() ;
  Mtbl_cf& cf_vr = *vr.set_spectral_va().c_cf ;
  cf_vr.annule_hard() ;
  het.set_etat_qcq() ;
  het.set_spectral_base(base) ;
  het.set_spectral_va().set_etat_cf_qcq() ;
  Mtbl_cf& cf_eta = *het.set_spectral_va().c_cf ;
  cf_eta.annule_hard() ;
  int taille = 4*nzm1 -2 ; // Two less than R3 case 
  Tbl sec_membre(taille) ; 
  Matrice systeme(taille, taille) ; 
  systeme.set_etat_qcq() ;
  int ligne ;  int colonne ;
  
  // Loop on l and m
  //----------------
  for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) {
	  base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	  if ((nullite_plm(j, nt, k, np, base) == 1)&&(l_q != 0)) {
		
	      ligne = 0 ;
	      colonne = 0 ;
	      systeme.annule_hard() ;
	      sec_membre.annule_hard() ;

// 	      //Nucleus 
// 	      nr = mg.get_nr(0) ;
// 	      alpha = mpaff->get_alpha()[0] ;
// 	      // value of at x=1 of eta ...
// 	      systeme.set(ligne, colonne) = sol_hom_un_eta.val_out_bound_jk(0, j, k) ;
// 	      systeme.set(ligne, colonne+1) = sol_hom_trois_eta.val_out_bound_jk(0, j, k) ;
// 	      sec_membre.set(ligne) = -sol_part_eta.val_out_bound_jk(0, j, k) ;
// 	      ligne++ ;
// 	      // ... and of its couterpart for V^r
// 	      systeme.set(ligne, colonne) = sol_hom_un_vr.val_out_bound_jk(0, j, k) ;
// 	      systeme.set(ligne, colonne+1) = sol_hom_trois_vr.val_out_bound_jk(0, j, k) ;
// 	      sec_membre.set(ligne) = -sol_part_vr.val_out_bound_jk(0,j,k) ; 
// 	      ligne++ ; 

//               //derivatives
// 	      int pari = (base_r == R_CHEBP ? 0 : 1) ;
// 	      for (int i=0; i<nr; i++) {
// 		  systeme.set(ligne, colonne) 
// 		      += (2*i+pari)*(2*i+pari)*sol_hom_un_eta(0, k, j, i)/alpha ;
// 		  systeme.set(ligne, colonne+1) 
// 		      += (2*i+pari)*(2*i+pari)*sol_hom_trois_eta(0, k, j, i)/alpha ;
// 		  sec_membre.set(ligne) 
// 		      -= (2*i+pari)*(2*i+pari)* sol_part_eta(0, k, j, i)/alpha ;
// 	      }
// 	      ligne++ ;
// 	      // ... and of its couterpart for V^r
// 	      for (int i=0; i<nr; i++) {
// 		  systeme.set(ligne, colonne) 
// 		      += (2*i+pari)*(2*i+pari)*sol_hom_un_vr(0, k, j, i)/alpha ;
// 		  systeme.set(ligne, colonne+1) 
// 		      += (2*i+pari)*(2*i+pari)*sol_hom_trois_vr(0, k, j, i)/alpha ;
// 		  sec_membre.set(ligne) 
// 		      -= (2*i+pari)*(2*i+pari)* sol_part_vr(0, k, j, i)/alpha ;
// 	      }
// 	      colonne += 2 ; 

      	      //shell 1 (boundary condition)
	      int zone = 1;
		  nr = mg.get_nr(zone) ;
		  alpha = mpaff->get_alpha()[zone] ;
	// 	  ligne -= 3 ;
		  // Here we prepare for a Robyn-type boundary condition on eta. 
		  // Parameters are dir_eta and neum_eta
		  systeme.set(ligne, colonne) 
		      = -dir_eta*sol_hom_un_eta.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      = -dir_eta * sol_hom_deux_eta.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      = -dir_eta*sol_hom_trois_eta.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      = -dir_eta*sol_hom_quatre_eta.val_in_bound_jk(zone, j, k)  ;
 
 


		  sec_membre.set(ligne) += dir_eta* sol_part_eta.val_in_bound_jk(zone, j, k) ;
		  //  ligne++ ;


		  int pari = -1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  -= neum_eta*pari*i*i*sol_hom_un_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  -= neum_eta*pari*i*i*sol_hom_deux_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  -= neum_eta*pari*i*i*sol_hom_trois_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  -= neum_eta*pari*i*i*sol_hom_quatre_eta(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  += neum_eta*pari*i*i* sol_part_eta(zone, k, j, i)/alpha ;
		      pari = -pari ;
		  }

		  sec_membre.set(ligne) -= (*bound_eta).val_in_bound_jk(1,j,k);
		  ligne++ ;
 
                 

		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) 
		      = - dir_vr*sol_hom_un_vr.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      = - dir_vr*sol_hom_deux_vr.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      = -dir_vr*sol_hom_trois_vr.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      = -dir_vr*sol_hom_quatre_vr.val_in_bound_jk(zone, j, k)  ;
		  sec_membre.set(ligne) += dir_vr*sol_part_vr.val_in_bound_jk(zone, j, k) ;
		  // ligne++ ;

	
	
		  pari = -1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  -= neum_vr*pari*i*i*sol_hom_un_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  -= neum_vr*pari*i*i*sol_hom_deux_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  -= neum_vr*pari*i*i*sol_hom_trois_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  -= neum_vr*pari*i*i*sol_hom_quatre_vr(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  += neum_vr*pari*i*i* sol_part_vr(zone, k, j, i)/alpha ;
		      pari = -pari ;
		  }

                  sec_membre.set(ligne) -= (*bound_mu).val_in_bound_jk(1,j,k);
		  ligne++ ;



	
		  // Values in 1: beginning with solution in eta...
			
		  systeme.set(ligne, colonne) 
		      += sol_hom_un_eta.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      += sol_hom_deux_eta.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      += sol_hom_trois_eta.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      += sol_hom_quatre_eta.val_out_bound_jk(zone, j, k)  ;
		  sec_membre.set(ligne) -= sol_part_eta.val_out_bound_jk(zone, j, k) ;
		  ligne++ ;
		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) 
		      += sol_hom_un_vr.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      += sol_hom_deux_vr.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      += sol_hom_trois_vr.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      += sol_hom_quatre_vr.val_out_bound_jk(zone, j, k)  ;
		  sec_membre.set(ligne) -= sol_part_vr.val_out_bound_jk(zone, j, k) ;
		  ligne++ ;

		  //derivative at 1 
		  pari = 1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  += pari*i*i*sol_hom_un_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  += pari*i*i*sol_hom_deux_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  += pari*i*i*sol_hom_trois_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  += pari*i*i*sol_hom_quatre_eta(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  -= pari*i*i* sol_part_eta(zone, k, j, i)/alpha ;
		      pari = pari ;
		  }
		  ligne++ ;
		  // ... and their couterparts for V^r
		  pari = 1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  += pari*i*i*sol_hom_un_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  += pari*i*i*sol_hom_deux_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  += pari*i*i*sol_hom_trois_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  += pari*i*i*sol_hom_quatre_vr(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  -= pari*i*i* sol_part_vr(zone, k, j, i)/alpha ;
		      pari = pari ;
		  }
		  colonne += 4 ;

		  // Other shells (WARNING! mapping must contain at least two shells!!)
		  if ( nz <= 2){
		    cout <<"WARNING!! Mapping must contain at least 2 shells!!" << endl;}
		  for (int zone=2 ; zone<nzm1 ; zone++) {

	  nr = mg.get_nr(zone) ;
		  alpha = mpaff->get_alpha()[zone] ;
		  ligne -= 3 ;
		  systeme.set(ligne, colonne) 
		      = -sol_hom_un_eta.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      = -sol_hom_deux_eta.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      = -sol_hom_trois_eta.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      = -sol_hom_quatre_eta.val_in_bound_jk(zone, j, k)  ;
		  sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(zone, j, k) ;
		  ligne++ ;
		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) 
		      = -sol_hom_un_vr.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      = -sol_hom_deux_vr.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      = -sol_hom_trois_vr.val_in_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      = -sol_hom_quatre_vr.val_in_bound_jk(zone, j, k)  ;
		  sec_membre.set(ligne) += sol_part_vr.val_in_bound_jk(zone, j, k) ;
		  ligne++ ;

		  //derivative of (x+echelle)^(l-1) at -1 
		  pari = -1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  -= pari*i*i*sol_hom_un_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  -= pari*i*i*sol_hom_deux_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  -= pari*i*i*sol_hom_trois_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  -= pari*i*i*sol_hom_quatre_eta(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  += pari*i*i* sol_part_eta(zone, k, j, i)/alpha ;
		      pari = -pari ;
		  }
		  ligne++ ;
		  // ... and their couterparts for V^r
		  pari = -1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  -= pari*i*i*sol_hom_un_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  -= pari*i*i*sol_hom_deux_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  -= pari*i*i*sol_hom_trois_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  -= pari*i*i*sol_hom_quatre_vr(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  += pari*i*i* sol_part_vr(zone, k, j, i)/alpha ;
		      pari = -pari ;
		  }
		  ligne++ ;
			
		  systeme.set(ligne, colonne) 
		      += sol_hom_un_eta.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      += sol_hom_deux_eta.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      += sol_hom_trois_eta.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      += sol_hom_quatre_eta.val_out_bound_jk(zone, j, k)  ;
		  sec_membre.set(ligne) -= sol_part_eta.val_out_bound_jk(zone, j, k) ;
		  ligne++ ;
		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) 
		      += sol_hom_un_vr.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+1) 
		      += sol_hom_deux_vr.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+2) 
		      += sol_hom_trois_vr.val_out_bound_jk(zone, j, k)  ;
		  systeme.set(ligne, colonne+3) 
		      += sol_hom_quatre_vr.val_out_bound_jk(zone, j, k)  ;
		  sec_membre.set(ligne) -= sol_part_vr.val_out_bound_jk(zone, j, k) ;
		  ligne++ ;

		  //derivative at 1 
		  pari = 1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  += pari*i*i*sol_hom_un_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  += pari*i*i*sol_hom_deux_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  += pari*i*i*sol_hom_trois_eta(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  += pari*i*i*sol_hom_quatre_eta(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  -= pari*i*i* sol_part_eta(zone, k, j, i)/alpha ;
		      pari = pari ;
		  }
		  ligne++ ;
		  // ... and their couterparts for V^r
		  pari = 1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  += pari*i*i*sol_hom_un_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+1) 
			  += pari*i*i*sol_hom_deux_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+2) 
			  += pari*i*i*sol_hom_trois_vr(zone, k, j, i)/alpha ;
		      systeme.set(ligne, colonne+3) 
			  += pari*i*i*sol_hom_quatre_vr(zone, k, j, i)/alpha ;
		      sec_membre.set(ligne) 
			  -= pari*i*i* sol_part_vr(zone, k, j, i)/alpha ;
		      pari = pari ;
		  }
		  colonne += 4 ;
		  }		    
	      		    
	      //Compactified external domain
	      nr = mg.get_nr(nzm1) ;

	      alpha = mpaff->get_alpha()[nzm1] ;
	      ligne -= 3 ;
	      //value of (x-1)^(l+2) at -1 :
	      systeme.set(ligne, colonne) 
		  -= sol_hom_deux_eta.val_in_bound_jk(nzm1, j, k) ;
	      systeme.set(ligne, colonne+1) 
		  -= sol_hom_quatre_eta.val_in_bound_jk(nzm1, j, k) ;
	      sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(nzm1, j, k) ;
	      //... and of its couterpart for V^r
	      systeme.set(ligne+1, colonne) 
		  -= sol_hom_deux_vr.val_in_bound_jk(nzm1, j, k) ;
	      systeme.set(ligne+1, colonne+1) 
		  -= sol_hom_quatre_vr.val_in_bound_jk(nzm1, j, k) ;
	      sec_membre.set(ligne+1) += sol_part_vr.val_in_bound_jk(nzm1, j, k) ;
			
	      ligne += 2 ;
	      //derivative of (x-1)^(l+2) at -1 :
		  pari = 1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  -= 4*alpha*pari*i*i*sol_hom_deux_eta(nzm1, k, j, i) ;
		      systeme.set(ligne, colonne+1) 
			  -= 4*alpha*pari*i*i*sol_hom_quatre_eta(nzm1, k, j, i) ;
		      sec_membre.set(ligne) 
			  += 4*alpha*pari*i*i* sol_part_eta(nzm1, k, j, i) ;
		      pari = -pari ;
		  }
		  ligne++ ;
		  // ... and their couterparts for V^r
		  pari = 1 ;
		  for (int i=0; i<nr; i++) {
		      systeme.set(ligne, colonne) 
			  -= 4*alpha*pari*i*i*sol_hom_deux_vr(nzm1, k, j, i) ;
		      systeme.set(ligne, colonne+1) 
			  -= 4*alpha*pari*i*i*sol_hom_quatre_vr(nzm1, k, j, i) ;
		      sec_membre.set(ligne) 
			  += 4*alpha*pari*i*i* sol_part_vr(nzm1, k, j, i) ;
		      pari = -pari ;
		  }
			
	      // Solution of the system giving the coefficients for the homogeneous 
	      // solutions
	      //-------------------------------------------------------------------
	      systeme.set_lu() ;
	      Tbl facteurs(systeme.inverse(sec_membre)) ;
	      int conte = 0 ;

	      // everything is put to the right place, the same combination of hom.
	      // solutions (with some l or -(l+1) factors) must be used for V^r
	      //-------------------------------------------------------------------
	      nr = mg.get_nr(0) ; //nucleus
	      for (int i=0 ; i<nr ; i++) {
		cf_eta.set(0, k, j, i) = 0.;
		cf_vr.set(0, k, j, i) = 0.;
	      }
	 //      conte += 2 ;
	      for (int zone=1 ; zone<nzm1 ; zone++) { //shells
		  nr = mg.get_nr(zone) ;
		  for (int i=0 ; i<nr ; i++) {
		      cf_eta.set(zone, k, j, i) = 
			  sol_part_eta(zone, k, j, i)
			  +facteurs(conte)*sol_hom_un_eta(zone, k, j, i) 
			  +facteurs(conte+1)*sol_hom_deux_eta(zone, k, j, i) 
			  +facteurs(conte+2)*sol_hom_trois_eta(zone, k, j, i) 
			  +facteurs(conte+3)*sol_hom_quatre_eta(zone, k, j, i) ;
		      cf_vr.set(zone, k, j, i) = sol_part_vr(zone, k, j, i)
			  +facteurs(conte)*sol_hom_un_vr(zone, k, j, i) 
			  +facteurs(conte+1)*sol_hom_deux_vr(zone, k, j, i) 
			  +facteurs(conte+2)*sol_hom_trois_vr(zone, k, j, i) 
			  +facteurs(conte+3)*sol_hom_quatre_vr(zone, k, j, i) ;
		  }
		  conte+=4 ;
	      }
	      nr = mg.get_nr(nz-1) ; //compactified external domain
	      for (int i=0 ; i<nr ; i++) {
		  cf_eta.set(nzm1, k, j, i) = sol_part_eta(nzm1, k, j, i)
		      +facteurs(conte)*sol_hom_deux_eta(nzm1, k, j, i) 
		      +facteurs(conte+1)*sol_hom_quatre_eta(nzm1, k, j, i) ;
		  cf_vr.set(nzm1, k, j, i) = sol_part_vr(nzm1, k, j, i)
		      +facteurs(conte)*sol_hom_deux_vr(nzm1, k, j, i) 
		      +facteurs(conte+1)*sol_hom_quatre_vr(nzm1, k, j, i) ;

	      }
	  } // End of nullite_plm  
      } //End of loop on theta
  vr.set_spectral_va().ylm_i() ;
  vr += vrl0 ;
  het.set_spectral_va().ylm_i() ;
  }

  
  Scalar pipo(*mpaff); 
  pipo.annule_hard(); 
  pipo.std_spectral_base(); 
  pipo.set_dzpuis(3);
  Param_elliptic param_mu(mu());


  

    resu.set_vr_eta_mu(vr, het, mu().sol_elliptic_boundary(param_mu, (*bound_mu), dir_mu , neum_mu)) ;
    //  resu.set_vr_eta_mu(vr, het, mu().poisson_dirichlet(vbound_mu2, 0)) ;


 
    Scalar mu2 =  mu().sol_elliptic_boundary(param_mu, *bound_mu, dir_mu , neum_mu);
    Scalar mu3 = mu2.laplacian();


  Scalar mu2b = mu().poisson_dirichlet(vbound_mu2, 0);
  Scalar mu3b = mu2b.laplacian();
 
 
  cout << "tests de résolution de l'équation en mu seule" << endl;
   
  cout << "difference absolue pour resolution" << endl;
 
  //  maxabs (mu() - mu3); 
 
  maxabs(mu() - mu3);
 

  cout << "decomposition spectrale de la condition aux bords" << endl;
 
  bound_mu->display();  
 
  cout << "solution obtenue" << endl;
  //  mu2.set_spectral_va().ylm();
  //  mu2.spectral_display();
 
  mu2.set_spectral_va().ylm();
  mu2.spectral_display();





  // resu.set_vr_eta_mu(vr, het, mu().poisson()) ;

  return ;
  
}
