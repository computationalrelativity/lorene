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
 * Revision 1.13  2004/05/07 15:33:22  j_novak
 * Treated the case where all components are null.
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
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "metric.h"
#include "tenseur.h"

Vector Vector::poisson(double lambda, const Metric_flat& met_f, int method) 
  const {
 
  bool nullite = true ;
  for (int i=0; i<3; i++) {
    assert(cmp[i]->check_dzpuis(4)) ;
    if (cmp[i]->get_etat() != ETATZERO) nullite = false ;
  }
  assert ((method>=0) && (method<4)) ;

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
      
      Scalar source_r = *(cmp[0]) - lambda*divf.dsdr(); 
      source_r.mult_r_dzpuis(3) ;
      source_r += 2*divf ;
      Scalar khi = source_r.poisson() ; 
      Scalar f_r = khi ;
      f_r.div_r() ; 
      
      Scalar source_eta = divf ;
      source_eta.mult_r_dzpuis(2) ;
      source_eta -= khi.dsdr() ;
      source_eta.dec_dzpuis(2) ;
      source_eta -= f_r ;
      Scalar eta = source_eta.poisson_angu() ;
      
      Scalar mu = div_free(met_f).mu().poisson() ;
      
      resu.set(1) = f_r ;
      resu.set(2) = eta.dsdt() - mu.stdsdp() ;
      resu.set(3) = eta.stdsdp() + mu.dsdt() ;
      
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

void Vector::poisson(const double lambda, Param& par, Vector& shift,
		     int method) const {

    
  bool nullite = true ;
  for (int i=0; i<3; i++) {
    assert(cmp[i]->check_dzpuis(4)) ;
    if (cmp[i]->get_etat() != ETATZERO) nullite = false ;
  }
  
  assert ((method==2)) ;

  if (nullite) shift.set_etat_zero() ;

  else {
  
    switch (method) {
	
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
      scal_auxi.allocate_all() ;
      scal_auxi.set_std_base() ;

      Tenseur resu_p(*mp, 1, CON, mp->get_bvect_cart() ) ;
      resu_p.set_etat_qcq() ;
      source_p.poisson_vect(lambda, par, resu_p, vect_auxi, scal_auxi) ;
      resu_p.change_triad(mp->get_bvect_spher() ) ;

      for (int i=1; i<=3; i++) 
	shift.set(i) = resu_p(i-1) ;

      break ;
    }

    default : {
      cout << "Vector::poisson : unexpected type of method !" << endl 
	   << "  method = " << method << endl ; 
      
      
      abort() ;
      break ; 
      
      
    } // End of switch  
    }
  }
}




