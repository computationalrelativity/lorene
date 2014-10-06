/*
 *  Method for vector Poisson equation inverting eqs. for V^r and eta as a block
 *  (with a boundary condition on the excised sphere).
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005  Jerome Novak
 *                       Francois Limousin
 *                       Jose Luis Jaramillo
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

char vector_poisson_boundary_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:21  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2005/06/09 08:00:09  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cassert>
#include <cstdlib>
#include <cmath>

// Lorene headers
#include "metric.h"
#include "diff.h"
#include "param_elliptic.h"
#include "proto.h"
#include "utilitaires.h"

void Vector::poisson_boundary(double lam, const Mtbl_cf& bound_vr, 
			    const Mtbl_cf& bound_eta, const Mtbl_cf& bound_mu, 
			      int num_front, double fact_dir, double fact_neu, 
			      Vector& resu) const {

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
    cout << "Not implemented yet !!" << endl ;
    abort() ;
    /*
      const Metric_flat& mets = mp->flat_met_spher() ;
      Vector_divfree sou(*mp, *triad, mets) ;
      for (int i=1; i<=3; i++) sou.set(i) = *cmp[i-1] ;
      resu = sou.poisson() ;
      return ;
    */
  }

  // Some working objects
  //---------------------
  const Mg3d& mg = *mpaff->get_mg() ;
  int nz = mg.get_nzone() ; int nzm1 = nz - 1;
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


  // The l_0 component is solved independently   // Understand this step 
  //------------------------------------------
  Scalar sou_l0 = (*cmp[0]) / (1. + lam) ;
  Param_elliptic param_l0(sou_l0) ;
  for (int l=0; l<nz; l++)
    param_l0.set_poisson_vect_r(l, true) ;
 
  //  Scalar vrl0 = sou_l0.sol_elliptic(param_l0) ;
  Scalar vrl0 = sou_l0.sol_elliptic_boundary(param_l0, bound_vr, fact_dir, fact_neu) ;

  // Build-up & inversion of the system for (eta, V^r) in each domain
  //-----------------------------------------------------------------

  // Shells
  //-------

  int nr ;
  int nt = mg.get_nt(0) ;
  int np = mg.get_np(0) ;
  int l_q = 0 ; int m_q = 0; int base_r = 0 ;
  double alpha, beta, ech ;

  assert(num_front+1 < nzm1) ; // Minimum one shell
  for (int zone=num_front+1 ; zone<nzm1 ; zone++) {                 //begins the loop on zones
      nr = mg.get_nr(zone) ;
      alpha = mpaff->get_alpha()[zone] ;
      beta = mpaff->get_beta()[zone] ;
      ech = beta / alpha ;
      assert (nr > 5) ;
      assert(nt == mg.get_nt(zone)) ;
      assert(np == mg.get_np(zone)) ;

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

	    // Building the operator   // Which is the eq. from the notes that it is actually implemented?  
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

	    //homogeneous solutions    //I do not understand!!! 
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
  alpha = mpaff->get_alpha()[nzm1] ;
  assert (nr > 4) ;
  int nr0 = nr - 2 ;  //two degrees of freedom less because of division by u^2

  // Loop on l and m
  //----------------
  for (int k=0 ; k<np+1 ; k++) {
    for (int j=0 ; j<nt ; j++) { 
      base.give_quant_numbers(nzm1, k, j, m_q, l_q, base_r) ;
      if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q != 0) ) {
	int dege1 = 1; //degeneracy of eq.1
	int dege2 = 0; //degeneracy of eq.2
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
	      = m2d2(lin,col) + 4*mxd(lin,col) 
	      + (2-(lam+1)*l_q*(l_q+1))*mid(lin,col) ;
	  for (int col=dege2; col<nr0; col++) 
	    oper.set(lin,col-dege2+nr_eq1) = 
	      -lam*mxd(lin,col) + 2*mid(lin,col) ;
	}
	for (int lin=0; lin<nr_eq2; lin++) {
	  for (int col=dege1; col<nr0; col++)
	    oper.set(lin+nr_eq1,col-dege1) 
	      = l_q*(l_q+1)*(lam*mxd(lin,col) + (3*lam+2)*mid(lin,col)) ;
	  for (int col=dege2; col<nr0; col++)
	    oper.set(lin+nr_eq1, col-dege2+nr_eq1) 
	      = (lam+1)*(m2d2(lin,col) + 4*mxd(lin,col)) 
	      - l_q*(l_q+1)*mid(lin,col) ;
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

	// Multiplication by u^2 
	//-----------------------
	Tbl res_v2(nr) ;  res_v2.set_etat_qcq() ;
	Tbl res_e2(nr) ; res_e2.set_etat_qcq() ;
	mult2_xm1_1d_cheb(nr, res_eta.t, res_e2.t) ; 
	mult2_xm1_1d_cheb(nr, res_vr.t, res_v2.t) ; 

	// Homogeneous solution (only 1/r^(l+2) and 1/r^l in the CED)
	Tbl sol_hom1 = solh(nr, l_q-1, 0., base_r) ;
	Tbl sol_hom2 = solh(nr, l_q+1, 0., base_r) ;
	for (int i=0 ; i<nr ; i++) {
	  sol_part_eta.set(nzm1, k, j, i) = alpha*alpha*res_e2(i) ;
	  sol_part_vr.set(nzm1, k, j, i) = alpha*alpha*res_v2(i) ;
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
  int taille = 4*(nzm1-num_front)-2  ; //## a verifier
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
		
	double f3_eta = lam*l_q + 3*lam + 2 ;
	double f4_eta = 2 + 2*lam - lam*l_q ;
	double f3_vr = (l_q+1)*(lam*l_q - 2) ;
	double f4_vr = l_q*(lam*l_q + lam + 2) ;
	ligne = 0 ;
	colonne = 0 ;
	sec_membre.annule_hard() ;
	for (int l=0; l<taille; l++) 
	  for (int c=0; c<taille; c++)
	    systeme.set(l,c) = 0 ;

	// First shell
	nr = mg.get_nr(num_front+1) ;
	alpha = mpaff->get_alpha()[num_front+1] ;
	double echelle = mpaff->get_beta()[num_front+1]/alpha ;
	// Conditions on eta (configuration space)
	//value and derivative of (x+echelle)^(l-1) at -1 
	systeme.set(ligne, colonne) = pow(echelle-1., double(l_q-1)) ;
	 
	// value and derivative of 1/(x+echelle) ^(l+2) at -1 
	systeme.set(ligne, colonne+1) = 1/pow(echelle-1., double(l_q+2)) ;

	//value and derivative of (x+echelle)^(l+1) at -1 
	systeme.set(ligne, colonne+2) = f3_eta*pow(echelle-1., double(l_q+1))  ;
	// value and derivative of 1/(x+echelle) ^l at -1 
	systeme.set(ligne, colonne+3) = f4_eta/pow(echelle-1., double(l_q))  ;
	for (int i=0 ; i<nr ; i++)
	  if (i%2 == 0)
	    sec_membre.set(ligne) -= sol_part_eta(num_front+1, k, j, i) ;
	  else sec_membre.set(ligne) += sol_part_eta(num_front+1, k, j, i) ;
	sec_membre.set(ligne) += bound_eta(num_front+1, k, j, 0) ;
	ligne++ ;

	// ... and their couterparts for V^r
	systeme.set(ligne, colonne) = fact_dir*l_q*pow(echelle-1., double(l_q-1)) + fact_neu*l_q*(l_q-1)*pow(echelle-1., double(l_q-2))/alpha ;
	systeme.set(ligne, colonne+1) = -fact_dir*(l_q+1)/pow(echelle-1., double(l_q+2)) + fact_neu*(l_q+1)*(l_q+2)/pow(echelle-1., double(l_q+3))/alpha ;
	systeme.set(ligne, colonne+2) = fact_dir*f3_vr*pow(echelle-1., double(l_q+1)) + fact_neu*f3_vr*(l_q+1)*pow(echelle-1., double(l_q))/alpha ;
	systeme.set(ligne, colonne+3) = fact_dir*f4_vr/pow(echelle-1., double(l_q)) - fact_neu*(f4_vr*l_q/pow(echelle-1., double(l_q+1)))/alpha ;
	for (int i=0 ; i<nr ; i++)
	  if (i%2 == 0)
	    sec_membre.set(ligne) -= fact_dir*sol_part_vr(num_front+1, k, j, i) - fact_neu*i*i/alpha*sol_part_vr(num_front+1, k, j, i) ;
	  else sec_membre.set(ligne) += fact_dir*sol_part_vr(num_front+1, k, j, i) - fact_neu*i*i/alpha*sol_part_vr(num_front+1, k, j, i)  ;
	sec_membre.set(ligne) += bound_vr(num_front+1, k, j, 0) ;
	
	ligne++ ;
	

	// Values at 1
	// eta
	//value of (x+echelle)^(l-1) at 1 
	systeme.set(ligne, colonne) = pow(echelle+1., double(l_q-1)) ;
	// value of 1/(x+echelle) ^(l+2) at 1 
	systeme.set(ligne, colonne+1) = 1./pow(echelle+1., double(l_q+2)) ;
	//value of (x+echelle)^(l+1) at 1 
	systeme.set(ligne, colonne+2) = f3_eta*pow(echelle+1., double(l_q+1));
	// value of 1/(x+echelle) ^l at 1 
	systeme.set(ligne, colonne+3) = f4_eta/pow(echelle+1., double(l_q)) ;
	for (int i=0 ; i<nr ; i++)
	  sec_membre.set(ligne) -= sol_part_eta(num_front+1, k, j, i) ;
	ligne++ ;
	// ... and their couterparts for V^r
	systeme.set(ligne, colonne) = l_q*pow(echelle+1., double(l_q-1)) ;
	systeme.set(ligne, colonne+1) 
	  = -double(l_q+1) / pow(echelle+1., double(l_q+2));
	systeme.set(ligne, colonne+2) = f3_vr*pow(echelle+1., double(l_q+1)) ;
	systeme.set(ligne, colonne+3) = f4_vr/pow(echelle+1., double(l_q));
	for (int i=0 ; i<nr ; i++)
	  sec_membre.set(ligne) -= sol_part_vr(num_front+1, k, j, i) ;
	ligne++ ;
	      
	//Derivatives at 1
	// eta
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
	  sec_membre.set(ligne) -= i*i/alpha*sol_part_eta(num_front+1, k, j, i) ;
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
	  sec_membre.set(ligne) -= i*i/alpha*sol_part_vr(num_front+1, k, j, i) ;
	      
	colonne += 4 ; // We pass to the next domain
      	    
	  
	// Next shells
	if (num_front+2<nzm1){
	  for (int zone=num_front+2 ; zone<nzm1 ; zone++) {
	    nr = mg.get_nr(zone) ;
	    alpha = mpaff->get_alpha()[zone] ;
	    echelle = mpaff->get_beta()[zone]/alpha ;
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
	if (taille > 2) systeme.set_band(5,5) ;
	else systeme.set_band(1,1) ;
	systeme.set_lu() ;
	Tbl facteurs(systeme.inverse(sec_membre)) ;
	int conte = 0 ;

	// everything is put to the right place, the same combination of hom.
	// solutions (with some l or -(l+1) factors) must be used for V^r
	//-------------------------------------------------------------------

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
	      -double(l_q+1)*facteurs(conte+1)*solution_hom_deux(zone, k, j, i)    // What is the origin of these factors?!
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
  
  Valeur temp_mu(mg.get_angu())  ;
  temp_mu = bound_mu ;
  const Valeur& limit_mu (temp_mu) ;

  resu.set_vr_eta_mu(vr, 0*het, mu().poisson_dirichlet(limit_mu, num_front)) ;

  return ; 
}


Vector Vector::poisson_dirichlet(double lam, const Valeur& bound_vr, 
			       const Valeur& bound_vt, const Valeur& bound_vp,
			       int num_front) const {

  Vector resu(*mp, CON, triad) ;
  resu = poisson_robin(lam, bound_vr, bound_vt, bound_vp, 1., 0., num_front) ;

  return resu ;

}

Vector Vector::poisson_neumann(double lam, const Valeur& bound_vr, 
			       const Valeur& bound_vt, const Valeur& bound_vp,
			       int num_front) const {

  Vector resu(*mp, CON, triad) ;
  resu = poisson_robin(lam, bound_vr, bound_vt, bound_vp, 0., 1., num_front) ;

  return resu ;

}

Vector Vector::poisson_robin(double lam, const Valeur& bound_vr, 
			     const Valeur& bound_vt, const Valeur& bound_vp,
			     double fact_dir, double fact_neu, 
			     int num_front) const {
 
 
  // Boundary condition for V^r    //Construction of a Mtbl_cf from Valeur with Ylm coefficients
  Valeur limit_vr (bound_vr) ; 

  limit_vr.coef() ;
  limit_vr.ylm() ; // Spherical harmonics transform.
  Mtbl_cf lim_vr (*(limit_vr.c_cf)) ;

  // bound_vt and bound_vp are only known at the boundary --> we fill 
  // all the zones extending the values at the boundary before calling to poisson_angu.
  Scalar temp_vt (*mp) ;
  Scalar temp_vp (*mp) ;
  temp_vt.annule_hard() ;
  temp_vp.annule_hard() ;
  int nz = mp->get_mg()->get_nzone() ;
  for (int l=0; l<nz; l++)
      for (int j=0; j<mp->get_mg()->get_nt(l); j++)
	  for (int k=0; k<mp->get_mg()->get_np(l); k++) {
	      temp_vt.set_grid_point(l, k, j, 0) = bound_vt(1, k, j, 0) ;
	      temp_vp.set_grid_point(l, k, j, 0) = bound_vp(1, k, j, 0) ;
	}
  temp_vt.set_spectral_va().set_base(bound_vt.base) ;     // We set the basis
  temp_vp.set_spectral_va().set_base(bound_vp.base) ;

  cout << "temp_vp" << endl << norme(temp_vp) << endl ;

  //Source for eta
  Scalar source_eta(*mp) ;
  Scalar vtstant (temp_vt) ;
  vtstant.div_tant() ;
  source_eta = temp_vt.dsdt() + vtstant + temp_vp.stdsdp() ;

  //Source for mu
  Scalar source_mu(*mp) ;
  Scalar vpstant (temp_vp) ;
  vpstant.div_tant() ;
  source_mu = temp_vp.dsdt() + vpstant - temp_vt.stdsdp() ;    //There was a wrong sign here

  Scalar temp_mu (source_mu.poisson_angu()) ;
  Scalar temp_eta (source_eta.poisson_angu()) ;

  // Boundary condition for mu
  Valeur limit_mu ((*mp).get_mg()->get_angu() )  ;
  int nnp = (*mp).get_mg()->get_np(1) ;
  int nnt = (*mp).get_mg()->get_nt(1) ;
  limit_mu= 1 ;
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      limit_mu.set(1, k, j, 0) = temp_mu.val_grid_point(1, k, j, 0) ;
  limit_mu.set_base(temp_mu.get_spectral_va().get_base()) ;

  limit_mu.coef() ;
  limit_mu.ylm() ; // Spherical harmonics transform.
  Mtbl_cf lim_mu (*(limit_mu.c_cf)) ;

  // Boundary condition for eta
  Valeur limit_eta ((*mp).get_mg()->get_angu() )  ;
  limit_eta = 1 ;
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      limit_eta.set(1, k, j, 0) = temp_eta.val_grid_point(1, k, j, 0) ;
  limit_eta.set_base(temp_eta.get_spectral_va().get_base()) ;

  limit_eta.coef() ;
  limit_eta.ylm() ; // Spherical harmonics transform.
  Mtbl_cf lim_eta (*(limit_eta.c_cf)) ;


  // Call to poisson_boundary(...)
  bool nullite = true ;
  for (int i=0; i<3; i++) {
    assert(cmp[i]->check_dzpuis(4)) ;
    if (cmp[i]->get_etat() != ETATZERO || bound_vr.get_etat() != ETATZERO || 
	bound_vt.get_etat() != ETATZERO || bound_vp.get_etat() != ETATZERO) 
      nullite = false ;
  }
  
  Vector resu(*mp, CON, triad) ;
  if (nullite)
    resu.set_etat_zero() ;
  else 
    poisson_boundary(lam, lim_vr, lim_eta, lim_mu, num_front, fact_dir, 
		     fact_neu, resu) ;

  return resu ;

}
