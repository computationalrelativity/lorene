/*
 * Computations of Scalar partial derivatives for a Map_log mapping
 */

/*
 *   Copyright (c) 2004 Philippe Grandclement
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
 * Revision 1.4  2016/12/05 16:17:58  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2012/01/17 10:33:43  j_penner
 * added a derivative with respect to the computational coordinate xi
 *
 * Revision 1.1  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * 
 * $Header$
 *
 */
 
// Header Lorene
#include "map.h"
#include "tensor.h"

                    //---------------------------------------------------
                   // d/d\xi
                  //---------------------------------------------------
namespace Lorene {
void Map_log::dsdxi(const Scalar& uu, Scalar& resu) const {

  assert (uu.get_etat() != ETATNONDEF) ; 
  assert (uu.get_mp().get_mg() == mg) ; 

    
  if (uu.get_etat() == ETATZERO) {
    resu.set_etat_zero() ; 
  }
  else {   
    assert( uu.get_etat() == ETATQCQ ) ; 

    const Valeur& uuva = uu.get_spectral_va() ; 

    uuva.coef() ;    // (uu.va).c_cf is up to date
	
    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ;

    if ( uu.get_dzpuis() == 0 ) {
      resu = uuva.dsdx() ;     //  dsdx == d/d\xi
	
      if (mg->get_type_r(nzm1) == UNSURR) {
	resu.set_dzpuis(2) ;    // r^2 d/dr has been computed in the
	                        // external domain
      } 
    }
    else {
      assert(mg->get_type_r(nzm1) == UNSURR) ;

      int dzp = uu.get_dzpuis() ;
      
      resu = uuva.dsdx() ;
      resu.annule_domain(nzm1) ;  // zero in the CED
      
      // Special treatment in the CED
      Valeur tmp_ced = - uuva.dsdx() ; 
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED
      tmp_ced.mult_xm1_zec() ; 
      tmp_ced.set(nzm1) -= dzp * uuva(nzm1) ; 
      
      // Recombination shells + CED : 
      resu.set_spectral_va() += tmp_ced ;
      
      resu.set_dzpuis(dzp+1) ;         
    
    }
    resu.set_spectral_base( uuva.dsdx().get_base() ) ; // same basis as d/dxi
  }   
}    

                    //---------------------------------------------------
                   // Derivee standard par rapport au "vrai" rayon .....
                  //---------------------------------------------------
void Map_log::dsdr(const Scalar& uu, Scalar& resu) const {

  assert (uu.get_etat() != ETATNONDEF) ; 
  assert (uu.get_mp().get_mg() == mg) ; 

    
  if (uu.get_etat() == ETATZERO) {
    resu.set_etat_zero() ; 
  }
  else {   
    assert( uu.get_etat() == ETATQCQ ) ; 

    const Valeur& uuva = uu.get_spectral_va() ; 

    uuva.coef() ;    // (uu.va).c_cf is up to date
	
    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ;

    if ( uu.get_dzpuis() == 0 ) {
      resu = uuva.dsdx() * dxdr ;     //  dxdr = dxi/dR, - dxi/dU (ZEC)
	
      if (mg->get_type_r(nzm1) == UNSURR) {
	resu.set_dzpuis(2) ;    // r^2 d/dr has been computed in the
	                        // external domain
      } 
    }
    else {
      assert(mg->get_type_r(nzm1) == UNSURR) ;

      int dzp = uu.get_dzpuis() ;
      
      resu = uuva.dsdx() * dxdr ;
      resu.annule_domain(nzm1) ;  // zero in the CED
      
      // Special treatment in the CED
      Valeur tmp_ced = - uuva.dsdx() ; 
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED
      tmp_ced.mult_xm1_zec() ; 
      tmp_ced.set(nzm1) -= dzp * uuva(nzm1) ; 
      
      // Recombination shells + CED : 
      resu.set_spectral_va() += tmp_ced ;
      
      resu.set_dzpuis(dzp+1) ;         
    
    }
    resu.set_spectral_base( uuva.dsdx().get_base() ) ; // same basis as d/dxi
  }   
}    

                    //---------------------------------------------------
                   // Derivee par rapport au rayon numerique (r ou lnr)
                  //---------------------------------------------------
void Map_log::dsdradial (const Scalar& uu, Scalar& resu) const {

  assert (uu.get_etat() != ETATNONDEF) ; 
  assert (uu.get_mp().get_mg() == mg) ; 

    
  if (uu.get_etat() == ETATZERO) {
    resu.set_etat_zero() ; 
  }
  else {   
    assert( uu.get_etat() == ETATQCQ ) ; 

    const Valeur& uuva = uu.get_spectral_va() ; 

    uuva.coef() ;    // (uu.va).c_cf is up to date
	
    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ;

    if ( uu.get_dzpuis() == 0 ) {
      resu = uuva.dsdx() * dxdlnr ;     //  dxdr = dxi/dR, - dxi/dU (ZEC)
	
      if (mg->get_type_r(nzm1) == UNSURR) {
	resu.set_dzpuis(2) ;    // r^2 d/dr has been computed in the
	                        // external domain
      } 
    }
    else {
      assert(mg->get_type_r(nzm1) == UNSURR) ;

      int dzp = uu.get_dzpuis() ;
      
      resu = uuva.dsdx() * dxdlnr ;
      resu.annule_domain(nzm1) ;  // zero in the CED
      
      // Special treatment in the CED
      Valeur tmp_ced = - uuva.dsdx() ; 
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED
      tmp_ced.mult_xm1_zec() ; 
      tmp_ced.set(nzm1) -= dzp * uuva(nzm1) ; 
      
      // Recombination shells + CED : 
      resu.set_spectral_va() += tmp_ced ;
      
      resu.set_dzpuis(dzp+1) ;         
    
    }
    resu.set_spectral_base( uuva.dsdx().get_base() ) ; // same basis as d/dxi
  }   
}
}
