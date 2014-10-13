/*
 *  Solution of the l=0 part of the vector Poisson equation (only r-component)
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

char pois_vect_r0_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:29  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2007/01/23 17:08:46  j_novak
 * New function pois_vect_r0.C to solve the l=0 part of the vector Poisson
 * equation, which involves only the r-component.
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "metric.h"
#include "proto.h"
#include "diff.h"

/*
 * This function solves for the l=0 component of
 * 
 *         d2 f   2 df   2f
 *         ---- + - -- - --  = source 
 *          dr2   r dr   r2
 * 
 * and returns the soluton f. 
 * The input Scalar must have dzpuis = 4.
 */
namespace Lorene {
Scalar pois_vect_r0(const Scalar& source) {

    const Map& map0 = source.get_mp() ;
    const Map_af* map1 = dynamic_cast<const Map_af*>(&map0) ;
    assert(map1 != 0x0) ;
    const Map_af& map = *map1 ;

    const Mg3d& mgrid = *map.get_mg() ;
    int nz = mgrid.get_nzone() ;

    Scalar resu(map) ;
    if (source.get_etat() == ETATZERO) {
	resu = 0 ;
	return resu ;
    }

    resu.annule_hard() ;
    resu.std_spectral_base_odd() ;
    resu.set_spectral_va().ylm() ;	
    Mtbl_cf& sol_coef = (*resu.set_spectral_va().c_cf) ;
    const Base_val& base = source.get_spectral_base() ;
    assert(resu.get_spectral_base() == base) ;
    assert(source.check_dzpuis(4)) ;

    Mtbl_cf sol_part(mgrid, base) ; sol_part.annule_hard() ;
    Mtbl_cf sol_hom1(mgrid, base) ; sol_hom1.annule_hard() ;
    Mtbl_cf sol_hom2(mgrid, base) ; sol_hom2.annule_hard() ;

    { int lz = 0 ;
      int nr = mgrid.get_nr(lz) ;
      double alpha2 = map.get_alpha()[lz]*map.get_alpha()[lz] ;
      assert(mgrid.get_type_r(lz) == RARE) ;
      int base_r = R_CHEBI ;
      Matrice ope(nr,nr) ;
      ope.annule_hard() ;
      Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
      Diff_sxdsdx sdx(base_r, nr) ; const Matrice& msdx = sdx.get_matrice() ;
      Diff_sx2 sx2(base_r, nr) ; const Matrice& ms2 = sx2.get_matrice() ;
      
      for (int lin=0; lin<nr-1; lin++) 
	  for (int col=0; col<nr; col++)
	      ope.set(lin,col) = (mdx2(lin,col) + 2*msdx(lin,col) - 2*ms2(lin,col))/alpha2 ;
      
      ope.set(nr-1, 0) = 1 ; //for the true homogeneous solution
      for (int i=1; i<nr; i++)
	  ope.set(nr-1, i) = 0 ;

      Tbl rhs(nr) ;
      rhs.annule_hard() ;
      for (int i=0; i<nr-1; i++)
	  rhs.set(i) = (*source.get_spectral_va().c_cf)(lz, 0, 0, i) ;
      rhs.set(nr-1) = 0 ;
      ope.set_lu() ;
      Tbl sol = ope.inverse(rhs) ;
		       
      for (int i=0; i<nr; i++)
	  sol_part.set(lz, 0, 0, i) = sol(i) ;
      
      rhs.annule_hard() ;
      rhs.set(nr-1) = 1 ;
      sol = ope.inverse(rhs) ;
      for (int i=0; i<nr; i++)
	  sol_hom1.set(lz, 0, 0, i) = sol(i) ;
      
    }

    for (int lz=1; lz<nz-1; lz++) {
	int nr = mgrid.get_nr(lz) ;
	double alpha = map.get_alpha()[lz] ;
	double beta = map.get_beta()[lz] ;
	double ech = beta / alpha ;
	assert(mgrid.get_type_r(lz) == FIN) ;
	int base_r = R_CHEB ;
			
	Matrice ope(nr,nr) ;
	ope.annule_hard() ;
	
	Diff_dsdx dx(base_r, nr) ; const Matrice& mdx = dx.get_matrice() ;
	Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
	Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;
	Diff_xdsdx xdx(base_r, nr) ; const Matrice& mxdx = xdx.get_matrice() ;
	Diff_xdsdx2 xdx2(base_r, nr) ; const Matrice& mxdx2 = xdx2.get_matrice() ;
	Diff_x2dsdx2 x2dx2(base_r, nr) ; const Matrice& mx2dx2 = x2dx2.get_matrice() ;
	
	for (int lin=0; lin<nr-2; lin++) 
	    for (int col=0; col<nr; col++) 
		ope.set(lin, col) = mx2dx2(lin, col) + 2*ech*mxdx2(lin, col) + ech*ech*mdx2(lin, col) 
		    + 2*(mxdx(lin, col) + ech*mdx(lin, col)) - 2*mid(lin, col) ;
	
	ope.set(nr-2, 0) = 0 ;
	ope.set(nr-2, 1) = 1 ;
	for (int col=2; col<nr; col++) {
	    ope.set(nr-2, col) = 0 ;
	}
	ope.set(nr-1, 0) = 1 ;
	for (int col=1; col<nr; col++)
	    ope.set(nr-1, col) = 0 ;
	
	Tbl src(nr) ;
	src.set_etat_qcq() ;
	for (int i=0; i<nr; i++) 
	    src.set(i) = (*source.get_spectral_va().c_cf)(lz, 0, 0, i) ;
	Tbl xsrc = src ; multx_1d(nr, &xsrc.t, base_r) ;
	Tbl x2src = src ; multx2_1d(nr, &x2src.t, base_r) ;
	Tbl rhs(nr) ;
	rhs.set_etat_qcq() ;
	for (int i=0; i<nr-2; i++) 
	    rhs.set(i) = beta*beta*src(i) + 2*beta*alpha*xsrc(i) + alpha*alpha*x2src(i) ;
	rhs.set(nr-2) = 0 ;
	rhs.set(nr-1) = 0 ;
	ope.set_lu() ;
	Tbl sol = ope.inverse(rhs) ;
	
	for (int i=0; i<nr; i++)
	    sol_part.set(lz, 0, 0, i) = sol(i) ;
	
	rhs.annule_hard() ;
	rhs.set(nr-2) = 1 ;
	sol = ope.inverse(rhs) ;
	for (int i=0; i<nr; i++)
	    sol_hom1.set(lz, 0, 0, i) = sol(i) ;
	
	rhs.set(nr-2) = 0 ;
	rhs.set(nr-1) = 1 ;
	sol = ope.inverse(rhs) ;
	for (int i=0; i<nr; i++)
	    sol_hom2.set(lz, 0, 0, i) = sol(i) ;
	
    }

    { int lz = nz-1 ;
      int nr = mgrid.get_nr(lz) ;
      double alpha2 = map.get_alpha()[lz]*map.get_alpha()[lz] ;
      assert(mgrid.get_type_r(lz) == UNSURR) ;
      int base_r = R_CHEBU ;
			
      Matrice ope(nr,nr) ;
      ope.annule_hard() ;
      Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
      Diff_sx2 sx2(base_r, nr) ; const Matrice& ms2 = sx2.get_matrice() ;
      
      for (int lin=0; lin<nr-3; lin++) 
	  for (int col=0; col<nr; col++) 
	      ope.set(lin, col) = (mdx2(lin, col) - 2*ms2(lin, col))/alpha2 ;
      
      for (int i=0; i<nr; i++) {
	  ope.set(nr-3, i) = i*i ; //for the finite part (derivative = 0 at infty)
      }
      
      for (int i=0; i<nr; i++) {
	  ope.set(nr-2, i) = 1 ; //for the limit at inifinity
      }
      
      ope.set(nr-1, 0) = 1 ; //for the true homogeneous solution
      for (int i=1; i<nr; i++)
	  ope.set(nr-1, i) = 0 ;
      
      Tbl rhs(nr) ;
      rhs.annule_hard() ;
      for (int i=0; i<nr-3; i++)
	  rhs.set(i) = (*source.get_spectral_va().c_cf)(lz, 0, 0, i) ;
      rhs.set(nr-3) = 0 ;
      rhs.set(nr-2) = 0 ;
      rhs.set(nr-1) = 0 ;
      ope.set_lu() ;
      Tbl sol = ope.inverse(rhs) ;
      for (int i=0; i<nr; i++)
	  sol_part.set(lz, 0, 0, i) = sol(i) ;
      
      rhs.annule_hard() ;
      rhs.set(nr-1) = 1 ;
      sol = ope.inverse(rhs) ;
      for (int i=0; i<nr; i++)
	  sol_hom2.set(lz, 0, 0, i) = sol(i) ;
      
    }

    Mtbl_cf dpart = sol_part ; dpart.dsdx() ;
    Mtbl_cf dhom1 = sol_hom1 ; dhom1.dsdx() ;
    Mtbl_cf dhom2 = sol_hom2 ; dhom2.dsdx() ;
    
    Matrice systeme(2*(nz-1), 2*(nz-1)) ;
    systeme.annule_hard() ;
    Tbl rhs(2*(nz-1)) ;
    rhs.annule_hard() ;

    //Nucleus
    int lin = 0 ;
    int col = 0 ;
    double alpha = map.get_alpha()[0] ;
    systeme.set(lin,col) = sol_hom1.val_out_bound_jk(0, 0, 0) ;
    rhs.set(lin) -= sol_part.val_out_bound_jk(0, 0, 0) ;
    lin++ ;
    systeme.set(lin,col) = dhom1.val_out_bound_jk(0, 0, 0) / alpha ;
    rhs.set(lin) -= dpart.val_out_bound_jk(0, 0, 0) / alpha ;
    col++ ;
    
    //Shells
    for (int lz=1; lz<nz-1; lz++) {
	alpha = map.get_alpha()[lz] ;
	lin-- ;
	systeme.set(lin,col) -= sol_hom1.val_in_bound_jk(lz, 0, 0) ;
	systeme.set(lin,col+1) -= sol_hom2.val_in_bound_jk(lz, 0, 0) ;
	rhs.set(lin) += sol_part.val_in_bound_jk(lz, 0, 0) ;
	
	lin++ ;
	systeme.set(lin,col) -= dhom1.val_in_bound_jk(lz, 0, 0) / alpha ;
	systeme.set(lin,col+1) -= dhom2.val_in_bound_jk(lz, 0, 0) / alpha ;
	rhs.set(lin) += dpart.val_in_bound_jk(lz, 0, 0) / alpha;
	
	lin++ ;
	systeme.set(lin, col) += sol_hom1.val_out_bound_jk(lz, 0, 0) ;
	systeme.set(lin, col+1) += sol_hom2.val_out_bound_jk(lz, 0, 0) ;
	rhs.set(lin) -= sol_part.val_out_bound_jk(lz, 0, 0) ;
	
	lin++ ;
	systeme.set(lin, col) += dhom1.val_out_bound_jk(lz, 0, 0) / alpha ;
	systeme.set(lin, col+1) += dhom2.val_out_bound_jk(lz, 0, 0) / alpha ;
	rhs.set(lin) -= dpart.val_out_bound_jk(lz, 0, 0) / alpha ;
	col += 2 ;
    }
    
    //CED
    alpha = map.get_alpha()[nz-1] ;
    lin-- ;
    systeme.set(lin,col) -= sol_hom2.val_in_bound_jk(nz-1, 0, 0) ;
    rhs.set(lin) += sol_part.val_in_bound_jk(nz-1, 0, 0) ;
    
    lin++ ;
    systeme.set(lin,col) -= (-4*alpha)*dhom2.val_in_bound_jk(nz-1, 0, 0) ;
    rhs.set(lin) += (-4*alpha)*dpart.val_in_bound_jk(nz-1, 0, 0) ;
    
    systeme.set_lu() ;
    Tbl coef = systeme.inverse(rhs) ;
    int indice = 0 ;
    
    for (int i=0; i<mgrid.get_nr(0); i++)
	sol_coef.set(0, 0, 0, i) = sol_part(0, 0, 0, i) 
	    + coef(indice)*sol_hom1(0, 0, 0, i) ;
    indice++ ;
    for (int lz=1; lz<nz-1; lz++) {
	for (int i=0; i<mgrid.get_nr(lz); i++)
	    sol_coef.set(lz, 0, 0, i) = sol_part(lz, 0, 0, i)
		+coef(indice)*sol_hom1(lz, 0, 0, i) 
		+coef(indice+1)*sol_hom2(lz, 0, 0, i) ;
	indice += 2 ;
    }
    for (int i=0; i<mgrid.get_nr(nz-1); i++)
	sol_coef.set(nz-1, 0, 0, i) = sol_part(nz-1, 0, 0, i)
	    +coef(indice)*sol_hom2(nz-1, 0, 0, i) ;

    delete resu.set_spectral_va().c ;
    resu.set_spectral_va().c = 0x0 ;
    resu.set_dzpuis(0) ;
    resu.set_spectral_va().ylm_i() ;
    
    return resu ; 
}
}
