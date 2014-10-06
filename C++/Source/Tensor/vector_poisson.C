/*
 *  Methods for solving vector Poisson equation.
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

char vector_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.23  2014/10/06 15:13:21  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.22  2005/02/14 13:01:50  j_novak
 * p_eta and p_mu are members of the class Vector. Most of associated functions
 * have been moved from the class Vector_divfree to the class Vector.
 *
 * Revision 1.21  2005/01/10 15:36:09  j_novak
 * In method 5: added transformation back from the Ylm base.
 *
 * Revision 1.20  2004/12/23 10:23:06  j_novak
 * Improved method 5 in the case when some components are zero.
 * Changed Vector_divfree::poisson to deduce eta from the equation. 
 *
 * Revision 1.19  2004/08/24 09:14:50  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.18  2004/07/27 09:40:05  j_novak
 * Yet another method for solving vector Poisson eq. with spherical coordinates.
 *
 * Revision 1.17  2004/06/14 15:15:58  j_novak
 * New method (No.4) to solve the vector Poisson eq. The output is continuous
 * for all (spherical) components.
 *
 * Revision 1.16  2004/05/26 07:30:59  j_novak
 * Version 1.15 was not good.
 *
 * Revision 1.15  2004/05/25 15:15:47  f_limousin
 * Function Vector::poisson with parameters now returns a Vector (the
 * result) instead of void.
 *
 * Revision 1.12  2004/03/26 17:05:24  j_novak
 * Added new method n.3 using Tenseur::poisson_vect_oohara
 *
 * Revision 1.11  2004/03/11 08:48:45  f_limousin
 * Implement method Vector::poisson with parameters, only with method
 * 2 yet.
 *
 * Revision 1.10  2004/03/10 16:38:38  e_gourgoulhon
 * Modified the prototype of poisson with param. to let it
 * agree with declaration in vector.h.
 *
 * Revision 1.9  2004/03/03 09:07:03  j_novak
 * In Vector::poisson(double, int), the flat metric is taken from the mapping.
 *
 * Revision 1.8  2004/02/24 17:00:25  j_novak
 * Added a forgotten term.
 *
 * Revision 1.7  2004/02/24 09:46:20  j_novak
 * Correction to cope with SGI compiler's warnings.
 *
 * Revision 1.6  2004/02/22 15:47:46  j_novak
 * Added 2 more methods to solve the vector poisson equation. Method 1 is not
 * tested yet.
 *
 * Revision 1.5  2004/02/20 10:53:41  j_novak
 * Minor modifs.
 *
 * Revision 1.4  2004/02/16 17:40:14  j_novak
 * Added a version of poisson with the flat metric as argument (avoids
 * unnecessary calculations by decompose_div)
 *
 * Revision 1.3  2003/10/29 11:04:34  e_gourgoulhon
 * dec2_dpzuis() replaced by dec_dzpuis(2).
 * inc2_dpzuis() replaced by inc_dzpuis(2).
 *
 * Revision 1.2  2003/10/22 13:08:06  j_novak
 * Better handling of dzpuis flags
 *
 * Revision 1.1  2003/10/20 15:15:42  j_novak
 * New method Vector::poisson().
 *
 *
 * $Headers: $
 *
 */

//C headers
#include <cstdlib>
#include <cmath>

// Lorene headers
#include "metric.h"
#include "tenseur.h"
#include "param.h"
#include "param_elliptic.h"
#include "proto.h"

Vector Vector::poisson(double lambda, const Metric_flat& met_f, int method) 
  const {
 
  bool nullite = true ;
  for (int i=0; i<3; i++) {
    assert(cmp[i]->check_dzpuis(4)) ;
    if (cmp[i]->get_etat() != ETATZERO) nullite = false ;
  }
  assert ((method>=0) && (method<7)) ;

  Vector resu(*mp, CON, triad) ;
  if (nullite)
    resu.set_etat_zero() ;
  else {

    switch (method) {
    
    case 0 : {

      Scalar poten(*mp) ;
      if (fabs(lambda+1) < 1.e-8)
	poten.set_etat_zero() ;
      else {
	poten = (potential(met_f) / (lambda + 1)).poisson() ;
      }
      
      Vector grad = poten.derive_con(met_f) ;
      grad.dec_dzpuis(2) ;
      
      return ( div_free(met_f).poisson() + grad) ;
      break ;
    }
      
    case 1 : {
      
      Scalar divf(*mp) ;
      if (fabs(lambda+1) < 1.e-8)
	divf.set_etat_zero() ;
      else {
	divf = (potential(met_f) / (lambda + 1)) ;
      }
      
      int nz = mp->get_mg()->get_nzone() ;

      //-----------------------------------
      // Removal of the l=0 part of div(F)
      //-----------------------------------
      Scalar div_lnot0 = divf ;
      div_lnot0.div_r_dzpuis(4) ;
      Scalar source_r(*mp) ;
      Valeur& va_div = div_lnot0.set_spectral_va() ;
      if (div_lnot0.get_etat() != ETATZERO) {
	va_div.coef() ;
	va_div.ylm() ;
	for (int lz=0; lz<nz; lz++) {
	  int np = mp->get_mg()->get_np(lz) ;
	  int nt = mp->get_mg()->get_nt(lz) ;
	  int nr = mp->get_mg()->get_nr(lz) ;
	  if (va_div.c_cf->operator()(lz).get_etat() != ETATZERO)
	    for (int k=0; k<np+1; k++) 
	      for (int j=0; j<nt; j++) {
		int l_q, m_q, base_r ;
		if (nullite_plm(j, nt, k, np, va_div.base) == 1) {
		  donne_lm(nz, lz, j, k, va_div.base, m_q, l_q, base_r) ;
		  if (l_q == 0) 
		    for (int i=0; i<nr; i++) 
		      va_div.c_cf->set(lz, k, j, i) = 0. ;
		}
	      }
	}
	source_r.set_etat_qcq() ;
	source_r.set_spectral_va() = 2*(*va_div.c_cf) ; //2*div(F)
	source_r.set_spectral_va().ylm_i() ;
	source_r.set_dzpuis(4) ;
      }
      else 
	source_r.set_etat_zero() ;
	   
      //------------------------
      // Other source terms ....
      //------------------------
      source_r += *(cmp[0]) - lambda*divf.dsdr() ;
      Scalar f_r(*mp) ;
      if (source_r.get_etat() != ETATZERO) {

	  //------------------------
	  // The elliptic operator
	  //------------------------
	  Param_elliptic param_fr(source_r) ;
	  for (int lz=0; lz<nz; lz++) 
	      param_fr.set_poisson_vect_r(lz) ;
	  
	  f_r = source_r.sol_elliptic(param_fr) ;
      }
      else
	  f_r.set_etat_zero() ;
            
      divf.dec_dzpuis(1) ;
      Scalar source_eta = divf - f_r.dsdr() ;
      source_eta.mult_r_dzpuis(0) ;
      source_eta -= 2*f_r ;
      resu.set_vr_eta_mu(f_r, source_eta.poisson_angu(), mu().poisson());
      
      break ;
      
    }

    case 2 : {
      
      Tenseur source_p(*mp, 1, CON, mp->get_bvect_spher() ) ;
      source_p.set_etat_qcq() ;
      for (int i=0; i<3; i++) {
	source_p.set(i) = Cmp(*cmp[i]) ;
      }
      source_p.change_triad(mp->get_bvect_cart()) ;
      Tenseur vect_auxi (*mp, 1, CON, mp->get_bvect_cart()) ;
      vect_auxi.set_etat_qcq() ;
      Tenseur scal_auxi (*mp) ;
      scal_auxi.set_etat_qcq() ;
      
      Tenseur resu_p(source_p.poisson_vect(lambda, vect_auxi, scal_auxi)) ;
      resu_p.change_triad(mp->get_bvect_spher() ) ;
      
      for (int i=1; i<=3; i++) 
	resu.set(i) = resu_p(i-1) ;
      
      break ;
    }
      
    case 3 : {
      
      Tenseur source_p(*mp, 1, CON, mp->get_bvect_spher() ) ;
      source_p.set_etat_qcq() ;
      for (int i=0; i<3; i++) {
	source_p.set(i) = Cmp(*cmp[i]) ;
      }
      source_p.change_triad(mp->get_bvect_cart()) ;
      Tenseur scal_auxi (*mp) ;
      scal_auxi.set_etat_qcq() ;
      
      Tenseur resu_p(source_p.poisson_vect_oohara(lambda, scal_auxi)) ;
      resu_p.change_triad(mp->get_bvect_spher() ) ;
      
      for (int i=1; i<=3; i++) 
	resu.set(i) = resu_p(i-1) ;
      
      break ;
    }

    case 4 : {
      Scalar divf(*mp) ;
      if (fabs(lambda+1) < 1.e-8)
	divf.set_etat_zero() ;
      else {
	divf = (potential(met_f) / (lambda + 1)) ;
      }
      
      int nz = mp->get_mg()->get_nzone() ;

      //-----------------------------------
      // Removal of the l=0 part of div(F)
      //-----------------------------------
      Scalar div_tmp = divf ;
      div_tmp.div_r_dzpuis(4) ;
      Scalar source_r(*mp) ;
      Valeur& va_div = div_tmp.set_spectral_va() ;
      if (div_tmp.get_etat() != ETATZERO) {
	va_div.ylm() ;
	for (int lz=0; lz<nz; lz++) {
	  int np = mp->get_mg()->get_np(lz) ;
	  int nt = mp->get_mg()->get_nt(lz) ;
	  int nr = mp->get_mg()->get_nr(lz) ;
	  if (va_div.c_cf->operator()(lz).get_etat() != ETATZERO)
	    for (int k=0; k<np+1; k++) 
	      for (int j=0; j<nt; j++) {
		int l_q, m_q, base_r ;
		if (nullite_plm(j, nt, k, np, va_div.base) == 1) {
		  donne_lm(nz, lz, j, k, va_div.base, m_q, l_q, base_r) ;
		  if (l_q == 0) 
		    for (int i=0; i<nr; i++) 
		      va_div.c_cf->set(lz, k, j, i) = 0. ;
		}
	      }
	}
	source_r.set_etat_qcq() ;
	source_r.set_spectral_va() = 2*(*va_div.c_cf) ; //2*div(F)
	source_r.set_spectral_va().ylm_i() ;
	source_r.set_dzpuis(4) ;
      }
      else 
	source_r.set_etat_zero() ;
	   
      //------------------------
      // Other source terms ....
      //------------------------
      source_r += *(cmp[0]) - lambda*divf.dsdr() ;
      Scalar f_r(*mp) ;
      if (source_r.get_etat() != ETATZERO) {

	  //------------------------
	  // The elliptic operator
	  //------------------------
	  Param_elliptic param_fr(source_r) ;
	  for (int lz=0; lz<nz; lz++) 
	      param_fr.set_poisson_vect_r(lz) ;
	  
	  f_r = source_r.sol_elliptic(param_fr) ;
      }
      else
	  f_r.set_etat_zero() ;

      //--------------------------
      // Equation for eta 
      //--------------------------
      Scalar source_eta = *cmp[1] ;
      source_eta.mult_sint() ;
      source_eta = source_eta.dsdt() ;
      source_eta.div_sint() ;
      source_eta = (source_eta + cmp[2]->stdsdp()).poisson_angu() ;

      Scalar dfrsr = f_r.dsdr() ;
      dfrsr.div_r_dzpuis(4) ;
      Scalar frsr2 = f_r ;
      frsr2.div_r_dzpuis(2) ;
      frsr2.div_r_dzpuis(4) ;

      Valeur& va_dfr = dfrsr.set_spectral_va() ;
      Valeur& va_fsr = frsr2.set_spectral_va() ;
      va_dfr.ylm() ;
      va_fsr.ylm() ;
            
      //Since the operator depends on the domain, the source
      //must be transformed accordingly.
      Valeur& va_eta = source_eta.set_spectral_va() ;
      if (source_eta.get_etat() == ETATZERO) source_eta.annule_hard() ;
      va_eta.ylm() ;
      for (int lz=0; lz<nz; lz++) {
	bool ced = (mp->get_mg()->get_type_r(lz) == UNSURR) ;
	int np = mp->get_mg()->get_np(lz) ;
	int nt = mp->get_mg()->get_nt(lz) ;
	int nr = mp->get_mg()->get_nr(lz) ;
	if (va_eta.c_cf->operator()(lz).get_etat() != ETATZERO)
	  for (int k=0; k<np+1; k++) 
	    for (int j=0; j<nt; j++) {
	      int l_q, m_q, base_r ;
	      if (nullite_plm(j, nt, k, np, va_eta.base) == 1) {
		donne_lm(nz, lz, j, k, va_eta.base, m_q, l_q, base_r) ;
		if (l_q > 0) 
		  for (int i=0; i<nr; i++) {
 		    if (va_div.c_cf->operator()(lz).get_etat() != ETATZERO) 
 		      va_eta.c_cf->set(lz, k, j, i) 
			-= (lambda + 2. / double(ced ? -l_q : (l_q+1) )) 
 			* va_div.c_cf->operator()(lz, k, j, i) ;
 		    if (va_fsr.c_cf->operator()(lz).get_etat() != ETATZERO) {
 		      va_eta.c_cf->set(lz, k, j, i) 
 			+= 2. / double(ced ? -l_q : (l_q+1) ) 
 			* va_dfr.c_cf->operator()(lz, k, j, i) ;
 		      va_eta.c_cf->set(lz, k, j, i)
 			-= 2.*( ced ? double(l_q+2)/double(l_q)
			       : double(l_q-1)/double(l_q+1) ) 
 			* va_fsr.c_cf->set(lz, k, j, i) ;
 		    }
		  } //Loop on r
	      } //nullite_plm
	    } //Loop on theta
      } // Loop on domains
      if (va_eta.c != 0x0) {
	delete va_eta.c;
	va_eta.c = 0x0 ;
      }
      va_eta.ylm_i() ;

      //------------------------
      // The elliptic operator
      //------------------------
      Param_elliptic param_eta(source_eta) ; 
      for (int lz=0; lz<nz; lz++) 
	param_eta.set_poisson_vect_eta(lz) ;

      resu.set_vr_eta_mu(f_r, source_eta.sol_elliptic(param_eta), mu().poisson()) ;
      break ;
      
    }
      
    case 5 : {
      
      Scalar divf(*mp) ;
      if (fabs(lambda+1) < 1.e-8)
	divf.set_etat_zero() ;
      else {
	divf = (potential(met_f) / (lambda + 1)) ;
      }
      
      int nz = mp->get_mg()->get_nzone() ;

      //-----------------------------------
      // Removal of the l=0 part of div(F)
      //-----------------------------------
      Scalar div_lnot0 = divf ;
      div_lnot0.div_r_dzpuis(4) ;
      Scalar source_r(*mp) ;
      Valeur& va_div = div_lnot0.set_spectral_va() ;
      if (div_lnot0.get_etat() != ETATZERO) {
	va_div.coef() ;
	va_div.ylm() ;
	for (int lz=0; lz<nz; lz++) {
	  int np = mp->get_mg()->get_np(lz) ;
	  int nt = mp->get_mg()->get_nt(lz) ;
	  int nr = mp->get_mg()->get_nr(lz) ;
	  if (va_div.c_cf->operator()(lz).get_etat() != ETATZERO)
	    for (int k=0; k<np+1; k++) 
	      for (int j=0; j<nt; j++) {
		int l_q, m_q, base_r ;
		if (nullite_plm(j, nt, k, np, va_div.base) == 1) {
		  donne_lm(nz, lz, j, k, va_div.base, m_q, l_q, base_r) ;
		  if (l_q == 0) 
		    for (int i=0; i<nr; i++) 
		      va_div.c_cf->set(lz, k, j, i) = 0. ;
		}
	      }
	}
	source_r.set_etat_qcq() ;
	source_r.set_spectral_va() = 2*(*va_div.c_cf) ; //2*div(F)
	source_r.set_spectral_va().ylm_i() ;
	source_r.set_dzpuis(4) ;
      }
      else 
	source_r.set_etat_zero() ;
	   
      //------------------------
      // Other source terms ....
      //------------------------
      source_r += *(cmp[0]) - lambda*divf.dsdr() ;
      Scalar f_r(*mp) ;
      if (source_r.get_etat() != ETATZERO) {

	  //------------------------
	  // The elliptic operator
	  //------------------------
	  
	  Param_elliptic param_fr(source_r) ;
	  for (int lz=0; lz<nz; lz++) 
	      param_fr.set_poisson_vect_r(lz) ;
	  
	  f_r = source_r.sol_elliptic(param_fr) ;
      }
      else
	  f_r.set_etat_zero() ;
	  
      Scalar source_eta = - *(cmp[0]) ;
      source_eta.mult_r_dzpuis(3) ;
      source_eta -= (lambda+2.)*divf ;
      source_eta.dec_dzpuis() ;
      f_r.set_spectral_va().ylm() ;
      Scalar tmp = 2*f_r + f_r.lapang() ;
      tmp.div_r_dzpuis(2) ;
      source_eta += tmp ;
      tmp = (1.+lambda)*divf ;
      tmp.mult_r_dzpuis(0) ;
      tmp += f_r ;
      source_eta = source_eta.primr() ;
      f_r.set_spectral_va().ylm_i() ;
      resu.set_vr_eta_mu(f_r, (tmp+source_eta).poisson_angu(), mu().poisson()) ;      
      break ;
      
    }

    case 6 : {
	
	poisson_block(lambda, resu) ;
	break ;

    }
    default : {
      cout << "Vector::poisson : unexpected type of method !" << endl 
	   << "  method = " << method << endl ; 
      abort() ;
      break ; 
    }
      
    } // End of switch  

  } // End of non-null case

  return resu ;

}

Vector Vector::poisson(double lambda, int method) const {
 
  const Base_vect_spher* tspher = dynamic_cast<const Base_vect_spher*>(triad) ;
  const Base_vect_cart* tcart = dynamic_cast<const Base_vect_cart*>(triad) ;

  assert ((tspher != 0x0) || (tcart != 0x0)) ;
  const Metric_flat* met_f = 0x0 ;

  if (tspher != 0x0) {
    assert (tcart == 0x0) ;
    met_f = &(mp->flat_met_spher()) ;
  }

  if (tcart != 0x0) {
    assert (tspher == 0x0) ;
    met_f = &(mp->flat_met_cart()) ;
  }

  return ( poisson(lambda, *met_f, method) );
    
}

// Version with parameters
// -----------------------

Vector Vector::poisson(const double lambda, Param& par, int method) const {

    
    for (int i=0; i<3; i++)
	assert(cmp[i]->check_dzpuis(4)) ;

    assert ((method==0) || (method==2)) ;

    Vector resu(*mp, CON, triad) ;

    switch (method) {
	
	case 0 : {

    Metric_flat met_local(*mp, *triad) ;
    int nitermax = par.get_int(0) ; 
    int& niter = par.get_int_mod(0) ; 
    double  relax = par.get_double(0) ; 
    double precis = par.get_double(1) ;     
    Cmp& ss_phi = par.get_cmp_mod(0) ;
    Cmp& ss_khi = par.get_cmp_mod(1) ;
    Cmp& ss_mu = par.get_cmp_mod(2) ;

    Param par_phi ; 
    par_phi.add_int(nitermax, 0) ;
    par_phi.add_int_mod(niter, 0) ;
    par_phi.add_double(relax, 0) ;
    par_phi.add_double(precis, 1) ;
    par_phi.add_cmp_mod(ss_phi, 0) ;
    
    Scalar poten(*mp) ;
    if (fabs(lambda+1) < 1.e-8)
	poten.set_etat_zero() ;
    else {
	Scalar tmp = potential(met_local) / (lambda + 1) ;
	tmp.inc_dzpuis(2) ;
	tmp.poisson(par_phi, poten) ;
    }
    
    Vector grad = poten.derive_con(met_local) ;
    grad.dec_dzpuis(2) ;
    
    Param par_free ;
    par_free.add_int(nitermax, 0) ;
    par_free.add_int_mod(niter, 0) ;
    par_free.add_double(relax, 0) ;
    par_free.add_double(precis, 1) ;
    par_free.add_cmp_mod(ss_khi, 0) ;
    par_free.add_cmp_mod(ss_mu, 1) ;
   
    return  div_free(met_local).poisson(par_free) + grad ;
    break ;
  }



	case 2 : {
	
    Tenseur source_p(*mp, 1, CON, mp->get_bvect_spher() ) ;
    source_p.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      source_p.set(i) = Cmp(*cmp[i]) ;
    }
    source_p.change_triad(mp->get_bvect_cart()) ;

    Tenseur vect_auxi (*mp, 1, CON, mp->get_bvect_cart()) ;
    vect_auxi.set_etat_qcq() ;
    for (int i=0; i<3 ;i++){
	vect_auxi.set(i) = 0. ;
    }
    Tenseur scal_auxi (*mp) ;
    scal_auxi.set_etat_qcq() ;
    scal_auxi.set().annule_hard() ;
    scal_auxi.set_std_base() ;

    Tenseur resu_p(*mp, 1, CON, mp->get_bvect_cart() ) ;
    resu_p.set_etat_qcq() ;
    source_p.poisson_vect(lambda, par, resu_p, vect_auxi, scal_auxi) ;
    resu_p.change_triad(mp->get_bvect_spher() ) ;

     for (int i=1; i<=3; i++) 
       resu.set(i) = resu_p(i-1) ;

     break ;
  }

  default : {
    cout << "Vector::poisson : unexpected type of method !" << endl 
	 << "  method = " << method << endl ; 


    abort() ;
    break ; 
  }

    }// End of switch  
  return resu ;

}




