/*
 *  Resolution of the divergence-free vector Poisson equation
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char vector_df_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2005/02/09 13:20:27  j_novak
 * Completely new way of solving the vector Poisson equation in the Div_free
 * case: the system is inverted "as a whole" for V^r and eta. This works only
 * with Map_af...
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>
#include <math.h>

// Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"
#include "param.h"

Vector_divfree Vector_divfree::poisson(Param& par ) const {

   // All this has a meaning only for spherical components:
#ifndef NDEBUG 
    const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
    assert(bvs != 0x0) ; 
#endif
    
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    double relax = par.get_double() ; 
    double precis = par.get_double(1) ;     
    Cmp& ss_khi = par.get_cmp_mod(0) ;
    Cmp& ss_mu = par.get_cmp_mod(1) ;
      
    // Solution for the r-component
    // ----------------------------
    
    Scalar source_r = *(cmp[0]) ; 
    source_r.mult_r() ; 
    
    Param par_khi ;
    par_khi.add_int(nitermax, 0) ;
    par_khi.add_int_mod(niter, 0) ;
    par_khi.add_double(relax, 0) ;
    par_khi.add_double(precis, 1) ;
    par_khi.add_cmp_mod(ss_khi, 0) ;

    Scalar khi (*mp) ;
    khi.set_etat_zero() ;

    source_r.poisson(par_khi, khi) ; 
    khi.div_r() ;   // khi now contains V^r
    
    // Solution for mu
    // ---------------
    
    Param par_mu ;
    par_mu.add_int(nitermax, 0) ;
    par_mu.add_int_mod(niter, 0) ;
    par_mu.add_double(relax, 0) ;
    par_mu.add_double(precis, 1) ;
    par_mu.add_cmp_mod(ss_mu, 0) ;

    Scalar mu_resu (*mp) ;
    mu_resu.set_etat_zero() ;
 
    mu().poisson(par_mu, mu_resu) ;

    // Final result
    // ------------
    
    Vector_divfree resu(*mp, *triad, *met_div) ; 
    
    resu.set_vr_mu(khi, mu_resu) ; 
    
    return resu ;

}

/*
 * In the case without parameters, first is solved the equation for mu and then
 * the system of equations for (eta, V^r) is inverted as a whole:
 * d2 eta / dr2 + 2/r d eta / dr - 1/r d V^r / dr = S(eta)
 * d V^r / dr + 2/r V^r - l(l+1)/r eta = 0 (div free condition)
 * 
 * There is no l=0 contribution (divergence free in all space!).
 * In practise the system is inverted for h(r) and v(r) , such that
 * eta = r^2 h and V^r = r^2 v in the nucleus and the shells,
 * in the compactified domain one has eta = u^2 h and V^r = u^2 v (where u=1/r); 
 * this method is used only to get particular solutions.
 * 
 * Homogeneous solutions are known analitically: r^(l-1) and/or 1/r^(l+2)
 * Once the matching of eta is done, there is no choice left for V^r
 * because of the div-free condition.
 */
Vector_divfree Vector_divfree::poisson() const {

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

  // Final result
  // ------------
    Vector_divfree resu(*mpaff, *triad, *met_div) ; 

  // Solution for mu
  // ---------------
    Scalar mu_resu = mu().poisson() ;
    
    Scalar f_r(*mpaff) ;
    if (cmp[0]->get_etat() == ETATZERO) { // no work needed ...
	f_r.set_etat_zero() ;
	resu.set_vr_mu(f_r, mu_resu) ; 
	return resu ;
    }
  
    // Some working objects
    //---------------------
  const Mg3d& mg = *mpaff->get_mg() ;
  int nz = mg.get_nzone() ; int nzm1 = nz - 1;
  assert(mg.get_type_r(0) == RARE) ;
  assert(mg.get_type_r(nzm1) == UNSURR) ;
  Scalar S_r = *cmp[0] ;
  Scalar S_eta = eta() ;
  S_r.set_spectral_va().ylm() ;
  S_eta.set_spectral_va().ylm() ;
  const Base_val& base = S_eta.get_spectral_va().base ;
  Mtbl_cf sol_part_eta(mg, base) ; sol_part_eta.annule_hard() ;
  Mtbl_cf sol_part_vr(mg, base) ; sol_part_vr.annule_hard() ;
  Mtbl_cf solution_hom_un(mg, base) ; solution_hom_un.annule_hard() ;
  Mtbl_cf solution_hom_deux(mg, base) ; solution_hom_deux.annule_hard() ;

  // Build-up & inversion of the system for (eta, V^r) in each domain
  //-----------------------------------------------------------------

  // Nucleus
  //--------
  int nr = mg.get_nr(0) ;
  int nt = mg.get_nt(0) ;
  int np = mg.get_np(0) ;
  double alpha = mpaff->get_alpha()[0] ;
  double beta = mpaff->get_beta()[0] ;
  int l_q = 0 ; int m_q = 0; int base_r = 0 ;
  int nr0 = nr - 1 ; //one degree of freedom less because of division by r^2

  // Loop on l and m
  //----------------
  for (int k=0 ; k<np+1 ; k++) {
      for (int j=0 ; j<nt ; j++) { 
	  base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	  if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q != 0) ) {
	      int dege1 = (l_q <3 ? 0 : 1) ; //degeneracy of eq.1
	      int dege2 = 0 ;  //degeneracy of eq.2
	      int nr_eq1 = nr0 - dege1 ; //Eq.1 is for h (eta/r^2)
	      int nr_eq2 = nr0 - dege2 ; //Eq.2 is the div-free condition
	      int nrtot = nr_eq1 + nr_eq2 ;
	      Matrice oper(nrtot, nrtot) ; oper.set_etat_qcq() ;
	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
	      Diff_x2dsdx2 d2(base_r, nr) ; const Matrice& md2 = d2.get_matrice() ;
	      Diff_xdsdx xd(base_r, nr) ; const Matrice& mxd = xd.get_matrice() ;
	      Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;

	      // Building the operator
	      //----------------------
	      for (int lin=0; lin<nr_eq1; lin++) { //eq.1 
		  for (int col=dege1; col<nr0; col++) 
		      oper.set(lin,col-dege1) 
			  = md2(lin,col) + 6*mxd(lin,col) + 6*mid(lin,col) ;
		  for (int col=dege2; col<nr0; col++) 
		      oper.set(lin,col-dege2+nr_eq1) = -mxd(lin,col) -2*mid(lin,col) ;
	      }
	      for (int lin=0; lin<nr0; lin++) { //eq.2
		  for (int col=dege1; col<nr0; col++)
		      oper.set(lin+nr_eq1,col-dege1) = -l_q*(l_q+1)*mid(lin,col) ;
		  for (int col=dege2; col<nr0; col++)
		      oper.set(lin+nr_eq1, col-dege2+nr_eq1) = mxd(lin,col) + 4*mid(lin,col) ;
	      }
	      oper.set_lu() ;

	      // Filling the r.h.s
	      //------------------
	      for (int i=0; i<nr_eq1; i++)  //eq.1
		  sec_membre.set(i) = (*S_eta.get_spectral_va().c_cf)(0, k, j, i) ;
	      for (int i=0; i<nr0; i++) //eq.2
		  sec_membre.set(i+nr_eq1) = 0 ;

	      // Inversion of the "big" operator
	      //--------------------------------
	      Tbl big_res = oper.inverse(sec_membre) ;
	      Tbl res_eta(nr) ;	 res_eta.set_etat_qcq() ;
	      Tbl res_vr(nr) ;  res_vr.set_etat_qcq() ;
	      
	      // Putting coefficients of h and v to individual arrays
	      //-----------------------------------------------------
	      for (int i=0; i<dege1; i++)
		  res_eta.set(i) = 0 ;
	      for (int i=dege1; i<nr0; i++)
		  res_eta.set(i) = big_res(i-dege1) ;
	      res_eta.set(nr0) = 0 ;
	      for (int i=0; i<dege2; i++)
		  res_vr.set(i) = 0 ;
	      for (int i=dege2; i<nr0; i++)
		  res_vr.set(i) = big_res(i-dege2+nr_eq1) ;
	      res_vr.set(nr0) = 0 ;

	      // Multiplication by xi^2 (the alpha^2 factor comes soon)
	      //-------------------------------------------------------
	      multx2_1d(nr, &res_eta.t, base_r) ;
	      multx2_1d(nr, &res_vr.t, base_r) ;

	      // Homogeneous solution (only r^(l-1) in the nucleus)
	      Tbl sol_hom = solh(nr, l_q-1, 0., base_r) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_part_eta.set(0, k, j, i) = alpha*alpha*res_eta(i) ;
		  sol_part_vr.set(0, k, j, i) = alpha*alpha*res_vr(i) ;
		  solution_hom_un.set(0, k, j, i) = sol_hom(i) ;
		  solution_hom_deux.set(0, k, j, i) = 0. ; 
	      }
	  }
      }
  }	    


  // Shells
  //-------
  for (int zone=1 ; zone<nzm1 ; zone++) {
      nr = mg.get_nr(zone) ; 
      assert (nr > 5) ;
      nr0 = nr - 2 ; //two degrees of freedom less because of division by r^2
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
	      int dege1 = 0 ; //degeneracy of eq.1
	      int dege2 = 1 ;  //degeneracy of eq.2
	      int nr_eq1 = nr0 - dege1 ; //Eq.1 is for eta
	      int nr_eq2 = nr0 - dege2 ; //Eq.2 is the div-free condition
	      int nrtot = nr_eq1 + nr_eq2 + 1;
	      Matrice oper(nrtot, nrtot) ; oper.set_etat_qcq() ;
	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
	      Diff_x2dsdx2 x2d2(base_r, nr) ; const Matrice& m2d2 = x2d2.get_matrice() ;
	      Diff_xdsdx2 xd2(base_r, nr) ; const Matrice& mxd2 = xd2.get_matrice() ;
	      Diff_dsdx2 d2(base_r, nr) ; const Matrice& md2 = d2.get_matrice() ;
	      Diff_xdsdx xd(base_r, nr) ; const Matrice& mxd = xd.get_matrice() ;
	      Diff_dsdx d1(base_r, nr) ; const Matrice& md = d1.get_matrice() ;
	      Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;

	      // Building the operator
	      //----------------------
	      for (int lin=0; lin<nr_eq1; lin++) { 
		  for (int col=dege1; col<nr0; col++) 
		      oper.set(lin,col-dege1) 
			  = m2d2(lin,col) + ech*(2*mxd2(lin,col) + ech*md2(lin,col))
			  +6*(mxd(lin,col) +ech*md(lin,col) + mid(lin,col)) ;
			  for (int col=dege2; col<nr0+1; col++) 
			      oper.set(lin,col-dege2+nr_eq1) = 
				  -mxd(lin,col) - ech*md(lin,col) - 2*mid(lin,col) ;
	      }
	      for (int lin=0; lin<nr_eq2; lin++) {
		  for (int col=dege1; col<nr0; col++)
		      oper.set(lin+nr_eq1,col-dege1) = -l_q*(l_q+1)*mid(lin,col) ;
		  for (int col=dege2; col<nr0+1; col++)
		      oper.set(lin+nr_eq1, col-dege2+nr_eq1) = 
			  mxd(lin,col) + ech*md(lin,col) + 4*mid(lin,col) ;
	      }
	      //Additional line to avoid spurious homogeneous solutions
	      //this line is the first one of the V^r eq.
	      for (int col=dege1; col<nr0; col++)
		  oper.set(nrtot-1, col-dege1) = 0 ;
	      for (int col=dege2; col<nr0+1; col++)
		  oper.set(nrtot-1, col-dege2+nr_eq1) = 
		      m2d2(0,col) + ech*(2*mxd2(0,col) + ech*md2(0,col))
		      +8*(mxd(0,col) +ech*md(0,col)) 
		      +(12 - l_q*(l_q+1))*mid(0,col) ;
	      oper.set_lu() ;
	      
	      // Filling the r.h.s
	      //------------------
	      for (int i=0; i<nr_eq1; i++) 
		  sec_membre.set(i) = (*S_eta.get_spectral_va().c_cf)(zone, k, j, i) ;
	      for (int i=0; i<nr_eq2; i++)
		  sec_membre.set(i+nr_eq1) = 0 ;
	      sec_membre.set(nr_eq1+nr_eq2) = 
		  (*S_r.get_spectral_va().c_cf)(zone, k, j, 0) ;

	      // Inversion of the "big" operator
	      //--------------------------------
	      Tbl big_res = oper.inverse(sec_membre) ;
	      Tbl res_eta(nr) ;	 res_eta.set_etat_qcq() ;
	      Tbl res_vr(nr) ;  res_vr.set_etat_qcq() ;
		  
	      // Putting coefficients of h and v to individual arrays
	      //-----------------------------------------------------
	      for (int i=0; i<dege1; i++)
		  res_eta.set(i) = 0 ;
	      for (int i=dege1; i<nr0; i++)
		  res_eta.set(i) = big_res(i-dege1) ;
	      res_eta.set(nr0) = 0 ;
	      res_eta.set(nr0+1) = 0 ;
	      for (int i=0; i<dege2; i++)
		  res_vr.set(i) = 0 ;
	      for (int i=dege2; i<nr0+1; i++)
		  res_vr.set(i) = big_res(i-dege2+nr_eq1) ;
	      res_vr.set(nr0+1) = 0 ;

	      // Multiplication by r^2 
	      //-----------------------
	      Tbl res_v1 = res_vr ; Tbl res_v2 = res_vr ;
	      Tbl res_e1 = res_eta ; Tbl res_e2 = res_eta ;
	      multx2_1d(nr, &res_e2.t, base_r) ; multx_1d(nr, &res_e1.t, base_r) ;
	      multx2_1d(nr, &res_v2.t, base_r) ; multx_1d(nr, &res_v1.t, base_r) ;
	      
	      //homogeneous solutions
	      Tbl sol_hom1 = solh(nr, l_q-1, ech, base_r) ;
	      Tbl sol_hom2 = solh(nr, l_q+1, ech, base_r) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_part_eta.set(zone, k, j, i) = alpha*alpha*res_e2(i) +
		      2*alpha*beta*res_e1(i) + beta*beta*res_eta(i) ;
		  sol_part_vr.set(zone, k, j, i) = alpha*alpha*res_v2(i) +
		      2*alpha*beta*res_v1(i) + beta*beta*res_vr(i) ;
		  solution_hom_un.set(zone, k, j, i) = sol_hom1(0,i) ;
		  solution_hom_deux.set(zone, k, j, i) = sol_hom2(1,i) ;
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
  alpha = mpaff->get_alpha()[nzm1] ;
  assert (nr > 4) ;
  nr0 = nr - 2 ;  //two degrees of freedom less because of division by r^2

  // Loop on l and m
  //----------------
  for (int k=0 ; k<np+1 ; k++) {
      for (int j=0 ; j<nt ; j++) { 
	  base.give_quant_numbers(nzm1, k, j, m_q, l_q, base_r) ;
	  if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q != 0) ) {
	      int dege1 = 0; //degeneracy of eq.1
	      int dege2 = 1; //degeneracy of eq.2
	      int nr_eq1 = nr0 - dege1 ; //Eq.1 is for eta
	      int nr_eq2 = nr0 - dege2 ; //Eq.2 is the div-free condition
	      int nrtot = nr_eq1 + nr_eq2 ;
	      Matrice oper(nrtot, nrtot) ; oper.set_etat_qcq() ;
	      Tbl sec_membre(nrtot) ; sec_membre.set_etat_qcq() ;
	      Diff_x2dsdx2 x2d2(base_r, nr) ; const Matrice& m2d2 = x2d2.get_matrice() ;
	      Diff_xdsdx xd(base_r, nr) ; const Matrice& mxd = xd.get_matrice() ;
	      Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;

	      // Building the operator
	      //----------------------
	      for (int lin=0; lin<nr_eq1; lin++) { 
		  for (int col=dege1; col<nr0; col++) 
		      oper.set(lin,col-dege1) 
			  = m2d2(lin,col) + 4*mxd(lin,col) + 2*mid(lin,col) ;
			  for (int col=dege2; col<nr0; col++) 
			      oper.set(lin,col-dege2+nr_eq1) = 
				  mxd(lin,col) + 2*mid(lin,col) ;
	      }
	      for (int lin=0; lin<nr_eq2; lin++) {
		  for (int col=dege1; col<nr0; col++)
		      oper.set(lin+nr_eq1,col-dege1) = l_q*(l_q+1)*mid(lin,col) ;
		  for (int col=dege2; col<nr0; col++)
		      oper.set(lin+nr_eq1, col-dege2+nr_eq1) = mxd(lin,col) ;
	      }
	      oper.set_lu() ;
	      
	      // Filling the r.h.s
	      //------------------
	      for (int i=0; i<nr_eq1; i++) 
		  sec_membre.set(i) = (*S_eta.get_spectral_va().c_cf)(nzm1, k, j, i) ;
	      for (int i=0; i<nr_eq2; i++)
		  sec_membre.set(i+nr_eq1) = 0 ;
	      Tbl big_res = oper.inverse(sec_membre) ;
	      Tbl res_eta(nr) ;	 res_eta.set_etat_qcq() ;
	      Tbl res_vr(nr) ;  res_vr.set_etat_qcq() ;
		  
	      // Putting coefficients of h and v to individual arrays
	      //-----------------------------------------------------
	      for (int i=0; i<dege1; i++)
		  res_eta.set(i) = 0 ;
	      for (int i=dege1; i<nr0; i++)
		  res_eta.set(i) = big_res(i-dege1) ;
	      res_eta.set(nr0) = 0 ;
	      res_eta.set(nr0+1) = 0 ;
	      for (int i=0; i<dege2; i++)
		  res_vr.set(i) = 0 ;
	      for (int i=dege2; i<nr0; i++)
		  res_vr.set(i) = big_res(i-dege2+nr_eq1) ;
	      res_vr.set(nr0) = 0 ;
	      res_vr.set(nr0+1) = 0 ;

	      // Multiplication by r^2 
	      //-----------------------
	      Tbl res_v2(nr) ;  res_v2.set_etat_qcq() ;
	      Tbl res_e2(nr) ; res_e2.set_etat_qcq() ;
	      mult2_xm1_1d_cheb(nr, res_eta.t, res_e2.t) ; 
	      mult2_xm1_1d_cheb(nr, res_vr.t, res_v2.t) ; 

	      // Homogeneous solution (only 1/r^(l+2) in the CED)
	      Tbl sol_hom = solh(nr, l_q+1, 0., base_r) ;
	      for (int i=0 ; i<nr ; i++) {
		  sol_part_eta.set(nzm1, k, j, i) = alpha*alpha*res_e2(i) ;
		  sol_part_vr.set(nzm1, k, j, i) = alpha*alpha*res_v2(i) ;
		  solution_hom_un.set(nzm1, k, j, i) = sol_hom(i) ;
		  solution_hom_deux.set(nzm1, k, j, i) = 0. ; 
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
  int taille = 2*nzm1 ;
  Tbl sec_membre(taille) ; 
  Matrice systeme(taille, taille) ; 
  systeme.set_etat_qcq() ;
  int ligne ;  int colonne ;
  
  for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) {
	  base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	  if ((nullite_plm(j, nt, k, np, base) == 1)||(l_q != 1)) {
		
	      ligne = 0 ;
	      colonne = 0 ;
	      sec_membre.annule_hard() ;
	      for (int l=0; l<taille; l++) 
		  for (int c=0; c<taille; c++)
		      systeme.set(l,c) = 0 ;
	      //Nucleus 
	      nr = mg.get_nr(0) ;
	      alpha = mpaff->get_alpha()[0] ;
	      // value of x^(l-1) at 1 :
	      systeme.set(ligne, colonne) = 1. ;
	      for (int i=0 ; i<nr ; i++)
		  sec_membre.set(ligne) -= sol_part_eta(0, k, j, i) ;
	      ligne++ ;
	      // derivative of x^(l-1) at 1 :
	      systeme.set(ligne, colonne) = double(l_q-1)/alpha ; 
	      if (base_r == R_CHEBP) 
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= double(4*i*i)/alpha
			  *sol_part_eta(0, k, j, i) ;
	      else
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= double((2*i+1)*(2*i+1))/alpha
					    *sol_part_eta(0, k, j, i) ;
	      colonne++ ; 
      	      //shells
	      for (int zone=1 ; zone<nzm1 ; zone++) {
		  nr = mg.get_nr(zone) ;
		  alpha = mpaff->get_alpha()[zone] ;
		  double echelle = mpaff->get_beta()[zone]/alpha ;
		  ligne -- ;
		  //value of (x+echelle)^(l-1) at -1 :
		  systeme.set(ligne, colonne) = -pow(echelle-1., double(l_q-1)) ;
		  // value of 1/(x+echelle) ^(l+2) at -1 :
		  systeme.set(ligne, colonne+1) = -1/pow(echelle-1., double(l_q+2)) ;
		  for (int i=0 ; i<nr ; i++)
		      if (i%2 == 0)
			  sec_membre.set(ligne) += sol_part_eta(zone, k, j, i) ;
		      else sec_membre.set(ligne) -= sol_part_eta(zone, k, j, i) ;
		  ligne++ ;
		  // derivative of (x+echelle)^(l-1) at -1 :
		  systeme.set(ligne, colonne) = 
		      -double(l_q-1)*pow(echelle-1, double(l_q-2))/alpha ;
		  if (l_q == 1) systeme.set(ligne, colonne) = 0 ;
		  // derivative of 1/(x+echelle)^(l+2) at -1 :
		  systeme.set(ligne, colonne+1) = 
		      double(l_q+2)/pow(echelle-1, double(l_q+3))/alpha ;
		  for (int i=0 ; i<nr ; i++)
		      if (i%2 == 0) sec_membre.set(ligne) -= 
					i*i/alpha*sol_part_eta(zone, k, j, i) ;
		      else  sec_membre.set(ligne) += 
				i*i/alpha*sol_part_eta(zone, k, j, i) ;
		  ligne++ ;
			
		  //value of (x+echelle)^(l-1) at 1 :
		  systeme.set(ligne, colonne) = pow(echelle+1., double(l_q-1)) ;
		  // value of 1/(x+echelle)^(l+2) at 1 :
		  systeme.set(ligne, colonne+1) = 1./pow(echelle+1., double(l_q+2)) ;
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= sol_part_eta(zone, k, j, i) ;
		  ligne ++ ;
		  //derivative of (x+echelle)^(l-1) at 1 :
		  systeme.set(ligne, colonne) = 
		      double(l_q-1)*pow(echelle+1., double(l_q-2))/alpha ;
		  if (l_q == 1) systeme.set(ligne, colonne) = 0 ;
		  //derivative of 1/(echelle+x)^(l+2) at 1 :
		  systeme.set(ligne, colonne+1) = 
		      -double(l_q+2)/pow(echelle+1., double(l_q+3))/alpha ;
		  for (int i=0 ; i<nr ; i++)
		      sec_membre.set(ligne) -= i*i/alpha*sol_part_eta(zone, k, j, i) ;
		  colonne += 2 ;
	      }		    
	      //Compactified external domain
	      nr = mg.get_nr(nzm1) ;
	      alpha = mpaff->get_alpha()[nzm1] ;
	      ligne -- ;
	      //value of (x-1)^(l+2) at -1 :
	      systeme.set(ligne, colonne) = -pow(-2, double(l_q+2)) ;
	      for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0) sec_membre.set(ligne) += sol_part_eta(nz-1, k, j, i) ;
		  else sec_membre.set(ligne) -= sol_part_eta(nz-1, k, j, i) ;
	      //derivative of (x-1)^(l+2) at -1 :
	      systeme.set(ligne+1, colonne) = 
		  alpha*double(l_q+2)*pow(-2., double(l_q+3)) ;
	      for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0) sec_membre.set(ligne+1) -= 
			-4*alpha*i*i*sol_part_eta(nz-1, k, j, i) ;
		  else sec_membre.set(ligne+1) += 
			   -4*alpha*i*i*sol_part_eta(nz-1, k, j, i) ;
			
	      // Solution of the system giving the coefficients for the homogeneous 
	      // solutions for eta
	      //-------------------------------------------------------------------
	      if (taille > 2) systeme.set_band(2,2) ;
	      else systeme.set_band(1,1) ;
	      systeme.set_lu() ;
	      Tbl facteurs(systeme.inverse(sec_membre)) ;
	      int conte = 0 ;

	      // everything is put to the right place, the same combination of hom.
	      // solutions (with some l or -(l+1) factors) must be used for V^r
	      //-------------------------------------------------------------------
	      nr = mg.get_nr(0) ; //nucleus
	      for (int i=0 ; i<nr ; i++) {
		  cf_eta.set(0, k, j, i) = sol_part_eta(0, k, j, i)
		      +facteurs(conte)*solution_hom_un(0, k, j, i) ;
		  cf_vr.set(0, k, j, i) = sol_part_vr(0, k, j, i)
		      +double(l_q)*facteurs(conte)*solution_hom_un(0, k, j, i) ;
	      }
	      conte++ ;
	      for (int zone=1 ; zone<nzm1 ; zone++) { //shells
		  nr = mg.get_nr(zone) ;
		  for (int i=0 ; i<nr ; i++) {
		      cf_eta.set(zone, k, j, i) = 
			  sol_part_eta(zone, k, j, i)
			  +facteurs(conte)*solution_hom_un(zone, k, j, i) 
			  +facteurs(conte+1)*solution_hom_deux(zone, k, j, i) ;
		      cf_vr.set(zone, k, j, i) = sol_part_vr(zone, k, j, i)
			  +double(l_q)*facteurs(conte)*solution_hom_un(zone, k, j, i) 
			  -double(l_q+1)*facteurs(conte+1)*solution_hom_deux(zone, k, j, i) ;
		  }
		  conte+=2 ;
	      }
	      nr = mg.get_nr(nz-1) ; //compactified external domain
	      for (int i=0 ; i<nr ; i++) {
		  cf_eta.set(nzm1, k, j, i) = sol_part_eta(nzm1, k, j, i)
		      +facteurs(conte)*solution_hom_un(nzm1, k, j, i) ;
		  cf_vr.set(nzm1, k, j, i) = sol_part_vr(nzm1, k, j, i)
		      -double(l_q+1)*facteurs(conte)*solution_hom_un(nzm1, k, j, i) ;
	      }
	  } // End of nullite_plm  
      } //End of loop on theta
      
  vr.set_spectral_va().ylm_i() ;
  het.set_spectral_va().ylm_i() ;

  resu.set_vr_eta_mu(vr, het, mu_resu) ;

  return resu ;

}

