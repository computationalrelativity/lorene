/*
 *  Method for vector Poisson equation inverting eqs. for V^r and eta as a block.
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005  Jerome Novak
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

char vector_poisson_block_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/12/30 13:39:38  j_novak
 * Changed the Galerkin base in the nucleus to (hopefully) stabilise the solver
 * when used in an iteration. Similar changes in the CED too.
 *
 * Revision 1.2  2005/02/15 15:43:18  j_novak
 * First version of the block inversion for the vector Poisson equation (method 6).
 *
 *
 * $Header$
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

void Vector::poisson_block(double lam, Vector& resu) const {

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
 	const Metric_flat& mets = mp->flat_met_spher() ;
 	Vector_divfree sou(*mp, *triad, mets) ;
 	for (int i=1; i<=3; i++) sou.set(i) = *cmp[i-1] ;
 	resu = sou.poisson() ;
 	return ;
     }

    // Some working objects
    //---------------------
  const Mg3d& mg = *mpaff->get_mg() ;
  int nz = mg.get_nzone() ; int nzm1 = nz - 1;
  assert(mg.get_type_r(0) == RARE) ;
  assert(mg.get_type_r(nzm1) == UNSURR) ;
  Scalar S_r = *cmp[0] ;
  if (S_r.get_etat() == ETATZERO) S_r.annule_hard() ;
  Scalar S_eta = eta() ;
  if (S_eta.get_etat() == ETATZERO) S_eta.annule_hard() ;
  S_r.set_spectral_va().ylm() ;
  S_eta.set_spectral_va().ylm() ;
  const Base_val& base = S_eta.get_spectral_va().base ;
  Mtbl_cf sol_part_eta(mg, base) ; sol_part_eta.annule_hard() ;
  Mtbl_cf sol_part_vr(mg, base) ; sol_part_vr.annule_hard() ;
  Mtbl_cf solution_hom_un(mg, base) ; solution_hom_un.annule_hard() ;
  Mtbl_cf solution_hom_deux(mg, base) ; solution_hom_deux.annule_hard() ;
  Mtbl_cf solution_hom_trois(mg, base) ; solution_hom_trois.annule_hard() ;
  Mtbl_cf solution_hom_quatre(mg, base) ; solution_hom_quatre.annule_hard() ;

  Scalar sou_l0 = (*cmp[0]) / (1. + lam) ;
  Param_elliptic param_l0(sou_l0) ;
  for (int l=0; l<nz; l++)
      param_l0.set_poisson_vect_r(l, true) ;
  Scalar vrl0 = sou_l0.sol_elliptic(param_l0) ;

  // Build-up & inversion of the system for (eta, V^r) in each domain
  //-----------------------------------------------------------------

  // Nucleus
  //--------
  int nr = mg.get_nr(0) ;
  int nt = mg.get_nt(0) ;
  int np = mg.get_np(0) ;
  double alpha = mpaff->get_alpha()[0] ; double alp2 = alpha*alpha ;
  double beta = mpaff->get_beta()[0] ;
  int l_q = 0 ; int m_q = 0; int base_r = 0 ;

  // Loop on l and m
  //----------------
  for (int k=0 ; k<np+1 ; k++) {
      for (int j=0 ; j<nt ; j++) { 
	  base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	  if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q != 0) ) {
	      int aa = 0 ; int bb = 0 ;  int nr0 = 0 ;
	      if (base_r == R_CHEBP) {
		  nr0 = nr-1 ;
		  aa = 0 ; bb = 1 ;
	      }
	      else {
		  assert (base_r == R_CHEBI) ;
		  nr0 = nr - 2 ;
		  aa = 2 ; bb = 1 ;
	      }
	      int d0 = nr - nr0 ;
	      int nrtot = 2*nr0 ;
	      Matrice oper(2*nr, 2*nr) ; oper.set_etat_qcq() ;
	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
	      Diff_dsdx2 d2(base_r, nr) ; const Matrice& md2 = d2.get_matrice() ;
	      Diff_sxdsdx xd(base_r, nr) ; const Matrice& mxd = xd.get_matrice() ;
	      Diff_sx2 s2(base_r, nr) ; const Matrice& ms2 = s2.get_matrice() ;

	      // Building the operator
	      //----------------------
	      for (int lin=0; lin<nr; lin++) { //eq.1 
		  for (int col=0; col<nr; col++) 
		      oper.set(lin,col) = 
			  (md2(lin,col) + 2*mxd(lin,col) 
			   -(lam+1)*l_q*(l_q+1)*ms2(lin,col)) / alp2 ;
		  for (int col=0; col<nr; col++) 
		      oper.set(lin,col+nr) = 
			  (lam*mxd(lin,col) + 2*(1+lam)*ms2(lin,col)) / alp2 ;
	      }
	      for (int lin=0; lin<nr; lin++) { //eq.2
		  for (int col=0; col<nr; col++)
		      oper.set(lin+nr,col) = 
			  (-lam*l_q*(l_q+1)*mxd(lin,col) 
			   +(lam+2)*l_q*(l_q+1)*ms2(lin,col)) / alp2 ;
		  for (int col=0; col<nr; col++)
		      oper.set(lin+nr, col+nr) = 
			  ((lam+1)*(md2(lin,col) + 2*mxd(lin,col)) 
			   - (2*(lam+1) + l_q*(l_q+1))*ms2(lin,col)) / alp2 ;
	      }
	      bool pb_eta = ( ( fabs( lam*double(l_q+3) + 2 ) < 0.01) && (l_q <=2) ) ;
	      if (!pb_eta) {
		  for (int col=0; col<nr; col++) 
		      oper.set(nr0-1, col) = 1 ;
		  for (int col=0; col<nr; col++) 
		      oper.set(nr0-1,nr+col) = 0 ;		  
	      }
	      if ((l_q > 2)||pb_eta) {
		  for (int col=0; col<nr; col++) 
		      oper.set(nr+nr0-1, col) = 0 ;
		  for (int col=0; col<nr; col++) 
		      oper.set(nr+nr0-1,nr+col) = 1 ;		  
	      }

	      Matrice op2(nrtot, nrtot) ;
	      op2.set_etat_qcq() ;
	      for (int i=0; i<nr0; i++) {
		  for (int col=0; col<nr0;col++) 
		      op2.set(i,col) = (aa*col+bb)*oper(i,col+1) 
			  + (aa*(col+1)+bb)*oper(i,col) ;
		  for (int col=0; col<nr0;col++) 
		      op2.set(i,col+nr0) = (aa*col+bb)*oper(i,col+nr+1)
			  + (aa*(col+1)+bb)*oper(i,col+nr) ;
	      }

	      for (int i=nr0; i<nrtot; i++) {
		  for (int col=0; col<nr0;col++) 
		      op2.set(i,col) = (aa*col+bb)*oper(i+d0,col+1) 
			  + (aa*(col+1)+bb)*oper(i+d0,col) ;
		  for (int col=0; col<nr0;col++) 
		      op2.set(i,col+nr0) = (aa*col+bb)*oper(i+d0,col+nr+1)
			  + (aa*(col+1)+bb)*oper(i+d0,col+nr) ;
	      }
 	      op2.set_lu() ;

	      // Filling the r.h.s
	      //------------------
	      for (int i=0; i<nr0; i++)  //eq.1
		  sec_membre.set(i) = (*S_eta.get_spectral_va().c_cf)(0, k, j, i) ;
	      if (!pb_eta) sec_membre.set(nr0-1) = 0 ;
	      for (int i=0; i<nr0; i++) //eq.2
		  sec_membre.set(i+nr0) 
		      = (*S_r.get_spectral_va().c_cf)(0, k, j, i) ;
	      if ((l_q > 2)||pb_eta) sec_membre.set(nrtot-1) = 0 ;

	      // Inversion of the "big" operator
	      //--------------------------------
	      Tbl big_res = op2.inverse(sec_membre) ;
	      Tbl res_eta(nr) ;	 res_eta.set_etat_qcq() ;
	      Tbl res_vr(nr) ;  res_vr.set_etat_qcq() ;
	      
	      // Putting coefficients of eta and Vr to individual arrays
	      //--------------------------------------------------------
	      res_eta.set(0) = (aa+bb)*big_res(0) ;
	      for (int i=1; i<nr0; i++)
		  res_eta.set(i) = (aa*(i-1)+bb)*big_res(i-1) 
		      + (aa*(i+1)+bb)*big_res(i);
	      res_eta.set(nr0) = big_res(nr0-1)*(aa*(nr0-1) + bb) ;
	      res_vr.set(0) = (aa+bb)*big_res(nr0) ;
	      for (int i=1; i<nr0; i++)
		  res_vr.set(i) = (aa*(i-1)+bb)*big_res(nr0+i-1) 
		      + (aa*(i+1)+bb)*big_res(nr0+i);
	      res_vr.set(nr0) = big_res(nrtot-1)*(aa*(nr0-1)+bb) ;
	      if (base_r == R_CHEBI) {
		  res_eta.set(nr-1) = 0 ;
		  res_vr.set(nr-1) = 0 ;
	      }

	      // Homogeneous solution (only r^(l-1) and r^(l+1) in the nucleus)
	      Tbl sol_hom1 = solh(nr, l_q-1, 0., base_r) ;
	      Tbl sol_hom2 = solh(nr, l_q+1, 0., base_r) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_part_eta.set(0, k, j, i) = res_eta(i) ;
		  sol_part_vr.set(0, k, j, i) = res_vr(i) ;
		  solution_hom_un.set(0, k, j, i) = sol_hom1(i) ;
		  solution_hom_deux.set(0, k, j, i) = 0. ; 
		  solution_hom_trois.set(0, k, j, i) = sol_hom2(i) ; 
		  solution_hom_quatre.set(0, k, j, i) = 0. ; 
	      }
	  }
      }
  }	    


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
	      int nrtot = nr_eq1 + nr_eq2 ;
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
		  for (int col=dege1; col<nr; col++) 
		      oper.set(lin,col-dege1) 
			  = m2d2(lin,col) + 2*ech*mxd2(lin,col) + ech*ech*md2(lin,col) 
			  + 2*(mxd(lin,col) + ech*md(lin,col)) 
			  - (lam+1)*l_q*(l_q+1)*mid(lin,col) ;
		  for (int col=dege2; col<nr; col++) 
		      oper.set(lin,col-dege2+nr_eq1) 
			  = lam*(mxd(lin,col) + ech*md(lin,col)) + 2*(1+lam)*mid(lin,col) ; 
	      }
	      for (int lin=0; lin<nr_eq2; lin++) {
		  for (int col=dege1; col<nr; col++)
		      oper.set(lin+nr_eq1,col-dege1) 
			  = -l_q*(l_q+1)*(lam*(mxd(lin,col) + ech*md(lin,col))
					  - (lam+2)*mid(lin,col)) ;
		  for (int col=dege2; col<nr; col++)
		      oper.set(lin+nr_eq1, col-dege2+nr_eq1) 
			  = (lam+1)*(m2d2(lin,col) + 2*ech*mxd2(lin,col) 
				     + ech*ech*md2(lin,col) 
				     + 2*(mxd(lin,col) + ech*md(lin,col)))
			  -(2*(lam+1)+l_q*(l_q+1))*mid(lin,col) ;
	      }
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
	      for (int i=0; i<nr_eq2; i++)
		  sec_membre.set(i+nr_eq1) = beta*beta*sr(i) 
		      + alpha*(alpha*x2sr(i) + 2*beta*xsr(i)) ;

	      // Inversion of the "big" operator
	      //--------------------------------
	      Tbl big_res = oper.inverse(sec_membre) ;
	      Tbl res_eta(nr) ;	 res_eta.set_etat_qcq() ;
	      Tbl res_vr(nr) ;  res_vr.set_etat_qcq() ;
		  
	      // Putting coefficients of h and v to individual arrays
	      //-----------------------------------------------------
	      for (int i=0; i<dege1; i++)
		  res_eta.set(i) = 0 ;
	      for (int i=dege1; i<nr; i++)
		  res_eta.set(i) = big_res(i-dege1) ;
	      for (int i=0; i<dege2; i++)
		  res_vr.set(i) = 0 ;
	      for (int i=dege2; i<nr; i++)
		  res_vr.set(i) = big_res(i-dege2+nr_eq1) ;

	      //homogeneous solutions
	      Tbl sol_hom1 = solh(nr, l_q-1, ech, base_r) ;
	      Tbl sol_hom2 = solh(nr, l_q+1, ech, base_r) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_part_eta.set(zone, k, j, i) = res_eta(i) ;
		  sol_part_vr.set(zone, k, j, i) = res_vr(i) ;
		  solution_hom_un.set(zone, k, j, i) = sol_hom1(0,i) ;
		  solution_hom_deux.set(zone, k, j, i) = sol_hom2(1,i) ;
		  solution_hom_trois.set(zone, k, j, i) = sol_hom2(0,i) ;
		  solution_hom_quatre.set(zone, k, j, i) = sol_hom1(1,i) ;
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
	      if ( fabs( lam*double(l_q+3) + 2 ) < 0.01) {
		  dege1 = dege2 ;
		  dege2 = 3 ;
	      }
	      int nr_eq1 = nr - dege1 ; //Eq.1 is for eta
	      int nr_eq2 = nr - dege2 ; //Eq.2 is the div-free condition
	      int nrtot = nr_eq1 + nr_eq2 ;
	      Matrice oper(nrtot, nrtot) ; oper.set_etat_qcq() ;
	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
	      Diff_dsdx2 d2(base_r, nr) ; const Matrice& md2 = d2.get_matrice() ;
	      Diff_sxdsdx sxd(base_r, nr) ; const Matrice& mxd = sxd.get_matrice() ;
	      Diff_sx2 sx2(base_r, nr) ; const Matrice& ms2 = sx2.get_matrice() ;

	      // Building the operator
	      //----------------------
	      for (int lin=0; lin<nr_eq1; lin++) { 
		  for (int col=dege1; col<nr; col++) 
		      oper.set(lin,col-dege1) 
			  = (md2(lin,col) - (lam+1)*l_q*(l_q+1)*ms2(lin,col))/alp2 ;
			  for (int col=dege2; col<nr; col++) 
			      oper.set(lin,col-dege2+nr_eq1) = 
				  (-lam*mxd(lin,col) + 2*(1+lam)*ms2(lin,col)) / alp2 ;
	      }
	      for (int lin=0; lin<nr_eq2; lin++) {
		  for (int col=dege1; col<nr; col++)
		      oper.set(lin+nr_eq1,col-dege1) 
			  = (l_q*(l_q+1)*(lam*mxd(lin,col) 
					  + (lam+2)*ms2(lin,col))) / alp2 ;
		  for (int col=dege2; col<nr; col++)
		      oper.set(lin+nr_eq1, col-dege2+nr_eq1) 
			  = ((lam+1)*md2(lin,col) 
			     - (2*(lam+1) + l_q*(l_q+1))*ms2(lin,col)) / alp2 ;
	      }
	      oper.set_lu() ;

	      // Filling the r.h.s
	      //------------------
	      for (int i=0; i<nr_eq1; i++) 
		  sec_membre.set(i) = (*S_eta.get_spectral_va().c_cf)(nzm1, k, j, i) ;
	      for (int i=0; i<nr_eq2; i++)
		  sec_membre.set(i+nr_eq1) =(*S_r.get_spectral_va().c_cf)(nzm1, k, j, i);
	      Tbl big_res = oper.inverse(sec_membre) ;
 	      Tbl res_eta(nr) ;	 res_eta.set_etat_qcq() ;
 	      Tbl res_vr(nr) ;  res_vr.set_etat_qcq() ;
		  
	      // Putting coefficients of h and v to individual arrays
	      //-----------------------------------------------------
	      for (int i=0; i<dege1; i++)
		  res_eta.set(i) = 0 ;
	      for (int i=dege1; i<nr; i++)
		  res_eta.set(i) = big_res(i-dege1) ;
	      for (int i=0; i<dege2; i++)
		  res_vr.set(i) = 0 ;
	      for (int i=dege2; i<nr; i++)
		  res_vr.set(i) = big_res(i-dege2+nr_eq1) ;
	      double somme = 0 ;
	      for (int i=0 ; i<nr ; i++)
		  somme += i*i*res_eta(i) ;
	      double somme_deux = somme ;
	      for (int i=0 ; i<nr ; i++)
		  somme_deux -= res_eta(i) ;
	      res_eta.set(1) = -somme ;
	      res_eta.set(0) = somme_deux ;
	      somme = 0 ;
	      for (int i=0 ; i<nr ; i++)
		  somme += i*i*res_vr(i) ;
	      somme_deux = somme ;
	      for (int i=0 ; i<nr ; i++)
		  somme_deux -= res_vr(i) ;
	      res_vr.set(1) = -somme ;
	      res_vr.set(0) = somme_deux ;


	      // Homogeneous solution (only 1/r^(l+2) and 1/r^l in the CED)
	      Tbl sol_hom1 = solh(nr, l_q-1, 0., base_r) ;
	      Tbl sol_hom2 = solh(nr, l_q+1, 0., base_r) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_part_eta.set(nzm1, k, j, i) = res_eta(i) ;
		  sol_part_vr.set(nzm1, k, j, i) = res_vr(i) ;
		  solution_hom_un.set(nzm1, k, j, i) = 0. ;
		  solution_hom_deux.set(nzm1, k, j, i) = sol_hom2(i) ;
		  solution_hom_trois.set(nzm1, k, j, i) = 0. ;
		  solution_hom_quatre.set(nzm1, k, j, i) = sol_hom1(i) ;
	      }
	  }	    
      }
  }

  // Now let's match everything ...
  //-------------------------------

  // Resulting V^r & eta
  Scalar vr(*mpaff) ; vr.set_etat_qcq() ;
  vr.set_spectral_base(base) ;
  vr.set_spectral_va().set_etat_cf_qcq() ;
  Mtbl_cf& cf_vr = *vr.set_spectral_va().c_cf ;
  cf_vr.annule_hard() ;
  Scalar het(*mpaff) ; het.set_etat_qcq() ;
  het.set_spectral_base(base) ;
  het.set_spectral_va().set_etat_cf_qcq() ;
  Mtbl_cf& cf_eta = *het.set_spectral_va().c_cf ;
  cf_eta.annule_hard() ;
  int taille = 4*nzm1 ;
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
		
	      double f3_eta = lam*double(l_q) + 3.*lam + 2. ;
	      double f4_eta = 2. + 2.*lam - lam*double(l_q) ;
	      double f3_vr = double(l_q+1)*(lam*double(l_q) - 2.) ;
	      double f4_vr = double(l_q)*(lam*double(l_q) + lam + 2.) ;
	      ligne = 0 ;
	      colonne = 0 ;
	      sec_membre.annule_hard() ;
	      for (int l=0; l<taille; l++) 
		  for (int c=0; c<taille; c++)
		      systeme.set(l,c) = 0 ;
	      //Nucleus 
	      nr = mg.get_nr(0) ;
	      alpha = mpaff->get_alpha()[0] ;
	      // value of x^(l-1) at 1 ...
	      systeme.set(ligne, colonne) = 1. ;
	      // value of x^(l+1) at 1 ...
	      systeme.set(ligne, colonne+1) = f3_eta ;
	      for (int i=0 ; i<nr ; i++)
		  sec_membre.set(ligne) -= sol_part_eta(0, k, j, i) ;
	      ligne++ ;
	      // ... and of its couterpart for V^r
	      systeme.set(ligne, colonne) = l_q;
	      systeme.set(ligne, colonne+1) = f3_vr ;
	      for (int i=0; i<nr; i++)
		  sec_membre.set(ligne) -= sol_part_vr(0,k,j,i) ; 
	      ligne++ ; //derivatives
	      // derivative of x^(l-1) at 1 ...
	      systeme.set(ligne, colonne) = double(l_q-1)/alpha ;
	      // derivative of x^(l+1) at 1 ...
	      systeme.set(ligne, colonne+1) = f3_eta*double(l_q+1)/alpha ;
	      if (base_r == R_CHEBP)
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= 4*i*i/alpha * sol_part_eta(0, k, j, i) ;
	      else
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= 
			  (2*i+1)*(2*i+1)/alpha * sol_part_eta(0, k, j, i) ;
	      ligne++ ;
	      // ... and of its couterpart for V^r
	      systeme.set(ligne, colonne) = l_q*double(l_q-1)/alpha ;
	      systeme.set(ligne, colonne+1) = f3_vr*double(l_q+1)/alpha ;
	      if (base_r == R_CHEBP)
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= 4*i*i/alpha * sol_part_vr(0, k, j, i) ;
	      else
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= 
			  (2*i+1)*(2*i+1)/alpha * sol_part_vr(0, k, j, i) ;
	      colonne += 2 ; 

      	      //shells
	      for (int zone=1 ; zone<nzm1 ; zone++) {
		  nr = mg.get_nr(zone) ;
		  alpha = mpaff->get_alpha()[zone] ;
		  double echelle = mpaff->get_beta()[zone]/alpha ;
		  ligne -= 3 ;
		  //value of (x+echelle)^(l-1) at -1 
		  systeme.set(ligne, colonne) = -pow(echelle-1., double(l_q-1)) ;
		  // value of 1/(x+echelle) ^(l+2) at -1 
		  systeme.set(ligne, colonne+1) = -1/pow(echelle-1., double(l_q+2)) ;
		  //value of (x+echelle)^(l+1) at -1 
		  systeme.set(ligne, colonne+2) = -f3_eta*pow(echelle-1., double(l_q+1));
		  // value of 1/(x+echelle) ^l at -1 
		  systeme.set(ligne, colonne+3) = -f4_eta/pow(echelle-1., double(l_q)) ;
		  for (int i=0 ; i<nr ; i++)
		      if (i%2 == 0)
			  sec_membre.set(ligne) += sol_part_eta(zone, k, j, i) ;
		      else sec_membre.set(ligne) -= sol_part_eta(zone, k, j, i) ;
		  ligne++ ;
		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) = -l_q*pow(echelle-1., double(l_q-1)) ;
		  systeme.set(ligne, colonne+1) = (l_q+1)/pow(echelle-1., double(l_q+2));
		  systeme.set(ligne, colonne+2) = -f3_vr*pow(echelle-1., double(l_q+1)) ;
		  systeme.set(ligne, colonne+3) = -f4_vr/pow(echelle-1., double(l_q));
 		  for (int i=0 ; i<nr ; i++)
		      if (i%2 == 0)
			  sec_membre.set(ligne) += sol_part_vr(zone, k, j, i) ;
		      else sec_membre.set(ligne) -= sol_part_vr(zone, k, j, i) ;
		  ligne++ ;

		  //derivative of (x+echelle)^(l-1) at -1 
		  systeme.set(ligne, colonne) 
		      = -(l_q-1)*pow(echelle-1., double(l_q-2))/alpha ;
		  // derivative of 1/(x+echelle) ^(l+2) at -1 
		  systeme.set(ligne, colonne+1) 
		      = (l_q+2)/pow(echelle-1., double(l_q+3))/alpha ;
		  // derivative of (x+echelle)^(l+1) at -1 
		  systeme.set(ligne, colonne+2) 
		      = -f3_eta*(l_q+1)*pow(echelle-1., double(l_q))/alpha;
		  // derivative of 1/(x+echelle) ^l at -1 
		  systeme.set(ligne, colonne+3) 
		      = (f4_eta*l_q/pow(echelle-1., double(l_q+1)))/alpha ;
		  for (int i=0 ; i<nr ; i++)
		      if (i%2 == 0) sec_membre.set(ligne) 
					-= i*i/alpha*sol_part_eta(zone, k, j, i) ;
		      else sec_membre.set(ligne) +=
			       i*i/alpha*sol_part_eta(zone, k, j, i) ;
		  ligne++ ;
		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) 
		      = -l_q*(l_q-1)*pow(echelle-1., double(l_q-2))/alpha ;
		  systeme.set(ligne, colonne+1) 
		      = -(l_q+1)*(l_q+2)/pow(echelle-1., double(l_q+3))/alpha ;
		  systeme.set(ligne, colonne+2) 
		      = -f3_vr*(l_q+1)*pow(echelle-1., double(l_q))/alpha ;
		  systeme.set(ligne, colonne+3) 
		      = (f4_vr*l_q/pow(echelle-1., double(l_q+1)))/alpha ;
		  for (int i=0 ; i<nr ; i++)
		      if (i%2 == 0) sec_membre.set(ligne) 
					-= i*i/alpha*sol_part_vr(zone, k, j, i) ;
		      else sec_membre.set(ligne) +=
			       i*i/alpha*sol_part_vr(zone, k, j, i) ;
		  ligne++ ;
			
		  //value of (x+echelle)^(l-1) at 1 
		  systeme.set(ligne, colonne) = pow(echelle+1., double(l_q-1)) ;
		  // value of 1/(x+echelle) ^(l+2) at 1 
		  systeme.set(ligne, colonne+1) = 1./pow(echelle+1., double(l_q+2)) ;
		  //value of (x+echelle)^(l+1) at 1 
		  systeme.set(ligne, colonne+2) = f3_eta*pow(echelle+1., double(l_q+1));
		  // value of 1/(x+echelle) ^l at 1 
		  systeme.set(ligne, colonne+3) = f4_eta/pow(echelle+1., double(l_q)) ;
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= sol_part_eta(zone, k, j, i) ;
		  ligne++ ;
		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) = l_q*pow(echelle+1., double(l_q-1)) ;
		  systeme.set(ligne, colonne+1) 
		      = -double(l_q+1) / pow(echelle+1., double(l_q+2));
		  systeme.set(ligne, colonne+2) = f3_vr*pow(echelle+1., double(l_q+1)) ;
		  systeme.set(ligne, colonne+3) = f4_vr/pow(echelle+1., double(l_q));
 		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= sol_part_vr(zone, k, j, i) ;
		  ligne++ ;

		  //derivative of (x+echelle)^(l-1) at 1 
		  systeme.set(ligne, colonne) 
		      = (l_q-1) * pow(echelle+1., double(l_q-2))/alpha ;
		  // derivative of 1/(x+echelle) ^(l+2) at 1 
		  systeme.set(ligne, colonne+1) 
		      = -(l_q+2) / pow(echelle+1., double(l_q+3))/alpha ;
		  // derivative of (x+echelle)^(l+1) at 1 
		  systeme.set(ligne, colonne+2) 
		      = f3_eta*(l_q+1) * pow(echelle+1., double(l_q))/alpha;
		  // derivative of 1/(x+echelle) ^l at 1 
		  systeme.set(ligne, colonne+3) 
		      = -f4_eta*l_q / pow(echelle+1., double(l_q+1))/alpha ;
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= i*i/alpha*sol_part_eta(zone, k, j, i) ;
		  ligne++ ;
		  // ... and their couterparts for V^r
		  systeme.set(ligne, colonne) 
		      = l_q*(l_q-1) * pow(echelle+1., double(l_q-2))/alpha ;
		  systeme.set(ligne, colonne+1) 
		      = (l_q+1)*(l_q+2) / pow(echelle+1., double(l_q+3))/alpha ;
		  systeme.set(ligne, colonne+2) 
		      = f3_vr*(l_q+1) * pow(echelle+1., double(l_q))/alpha ;
		  systeme.set(ligne, colonne+3) 
		      = -f4_vr*l_q / pow(echelle+1., double(l_q+1))/alpha ;
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= i*i/alpha*sol_part_vr(zone, k, j, i) ;

		  colonne += 4 ;
	      }		    
	      //Compactified external domain
	      nr = mg.get_nr(nzm1) ;

	      alpha = mpaff->get_alpha()[nzm1] ;
	      ligne -= 3 ;
	      //value of (x-1)^(l+2) at -1 :
	      systeme.set(ligne, colonne) = -pow(-2, double(l_q+2)) ;
	      //value of (x-1)^l at -1 :
	      systeme.set(ligne, colonne+1) = -f4_eta*pow(-2, double(l_q)) ;
	      for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0) sec_membre.set(ligne) += sol_part_eta(nzm1, k, j, i) ;
		  else sec_membre.set(ligne) -= sol_part_eta(nzm1, k, j, i) ;
	      //... and of its couterpart for V^r
	      systeme.set(ligne+1, colonne) = double(l_q+1)*pow(-2, double(l_q+2)) ;
	      systeme.set(ligne+1, colonne+1) = -f4_vr*pow(-2, double(l_q)) ;
	      for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0) sec_membre.set(ligne+1) += sol_part_vr(nzm1, k, j, i) ;
		  else sec_membre.set(ligne+1) -= sol_part_vr(nzm1, k, j, i) ;
			
	      ligne += 2 ;
	      //derivative of (x-1)^(l+2) at -1 :
	      systeme.set(ligne, colonne) = alpha*(l_q+2)*pow(-2, double(l_q+3)) ;
	      //derivative of (x-1)^l at -1 :
	      systeme.set(ligne, colonne+1) = alpha*l_q*f4_eta*pow(-2, double(l_q+1)) ;
	      for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0) sec_membre.set(ligne) 
				    -= -4*alpha*i*i*sol_part_eta(nzm1, k, j, i) ;
		  else sec_membre.set(ligne) 
			   += -4*alpha*i*i*sol_part_eta(nzm1, k, j, i) ;
	      //... and of its couterpart for V^r
	      systeme.set(ligne+1, colonne) 
		  = -alpha*double((l_q+1)*(l_q+2))*pow(-2, double(l_q+3)) ;
	      systeme.set(ligne+1, colonne+1) 
		  = alpha*double(l_q)*f4_vr*pow(-2, double(l_q+1)) ;
	      for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0) sec_membre.set(ligne+1) 
				    -= -4*alpha*i*i*sol_part_vr(nzm1, k, j, i) ;
		  else sec_membre.set(ligne+1) 
			   += -4*alpha*i*i*sol_part_vr(nzm1, k, j, i) ;
			
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
		  cf_eta.set(0, k, j, i) = sol_part_eta(0, k, j, i)
		      +facteurs(conte)*solution_hom_un(0, k, j, i) 
		      +facteurs(conte+1)*f3_eta*solution_hom_trois(0, k, j, i) ;
		  cf_vr.set(0, k, j, i) = sol_part_vr(0, k, j, i)
		      +double(l_q)*facteurs(conte)*solution_hom_un(0, k, j, i) 
		      +facteurs(conte+1)*f3_vr*solution_hom_trois(0, k, j, i) ;
	      }
	      conte += 2 ;
	      for (int zone=1 ; zone<nzm1 ; zone++) { //shells
		  nr = mg.get_nr(zone) ;
		  for (int i=0 ; i<nr ; i++) {
		      cf_eta.set(zone, k, j, i) = 
			  sol_part_eta(zone, k, j, i)
			  +facteurs(conte)*solution_hom_un(zone, k, j, i) 
			  +facteurs(conte+1)*solution_hom_deux(zone, k, j, i) 
			  +facteurs(conte+2)*f3_eta*solution_hom_trois(zone, k, j, i) 
			  +facteurs(conte+3)*f4_eta*solution_hom_quatre(zone, k, j, i) ;
		      cf_vr.set(zone, k, j, i) = sol_part_vr(zone, k, j, i)
			  +double(l_q)*facteurs(conte)*solution_hom_un(zone, k, j, i) 
			  -double(l_q+1)*facteurs(conte+1)*solution_hom_deux(zone, k, j, i) 
			  +f3_vr*facteurs(conte+2)*solution_hom_trois(zone, k, j, i) 
			  +f4_vr*facteurs(conte+3)*solution_hom_quatre(zone, k, j, i) ;
		  }
		  conte+=4 ;
	      }
	      nr = mg.get_nr(nz-1) ; //compactified external domain
	      for (int i=0 ; i<nr ; i++) {
		  cf_eta.set(nzm1, k, j, i) = sol_part_eta(nzm1, k, j, i)
		      +facteurs(conte)*solution_hom_deux(nzm1, k, j, i) 
		      +f4_eta*facteurs(conte+1)*solution_hom_quatre(nzm1, k, j, i) ;
		  cf_vr.set(nzm1, k, j, i) = sol_part_vr(nzm1, k, j, i)
		      -double(l_q+1)*facteurs(conte)*solution_hom_deux(nzm1, k, j, i) 
		      +f4_vr*facteurs(conte+1)*solution_hom_quatre(nzm1, k, j, i) ;

	      }
	  } // End of nullite_plm  
      } //End of loop on theta
  vr.set_spectral_va().ylm_i() ;
  vr += vrl0 ;
  het.set_spectral_va().ylm_i() ;

  resu.set_vr_eta_mu(vr, het, mu().poisson()) ;

  return ;

}
