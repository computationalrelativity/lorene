/*
 * Computations of Scalar partial derivatives for a Map_star mapping
 */

/*
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

// Header Lorene
#include "map.h"
#include "valeur.h"
#include "tensor.h"




namespace Lorene{
			//---------------------//
			//        d/dr         //
			//---------------------//


void Map_star::dsdr(const Scalar& uu, Scalar& resu) const {

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
      int dzp = uu.get_dzpuis() ;
      
      resu = uuva.dsdx() * dxdr ;
      if (mg->get_type_r(nzm1) == UNSURR) {
	  resu.annule_domain(nzm1) ;  // zero in the CED
      
	  // Special treatment in the CED
	  Valeur tmp_ced = - uuva.dsdx() ; 
	  tmp_ced.annule(0, nz-2) ; // only non zero in the CED
	  tmp_ced.mult_xm1_zec() ; 
	  tmp_ced.set(nzm1) -= dzp * uuva(nzm1) ; 
	  
	  // Recombination shells + CED : 
	  resu.set_spectral_va() += tmp_ced ;
      }
      resu.set_dzpuis(dzp+1) ;         
    
    }

    resu.set_spectral_base( uuva.dsdx().get_base() ) ; // same basis as d/dxi

  }
    
}

			//-----------------------------------//
			//       1/sin(theta) d/dphi         //
			//-----------------------------------//
void Map_star::stdsdp(const Scalar& ci, Scalar& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp().get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
		resu.set_etat_zero() ; 
    }
    else {

		assert( ci.get_etat() == ETATQCQ ) ; 
    
    // Computation of 1/sin(theta) df/dphi'   ---> stdfdp
		// ----------------------------
		
		const Valeur& stdfdp = ci.get_spectral_va().stdsdp() ;
	

		// Computation of 1/(dR/dxi) 1/sin(theta) dR/dphi' df/dx   ----> adfdx
		// -------------------------------------------

		Valeur adfdx = ci.get_spectral_va().dsdx() ; 	// df/dx

		Base_val sauve_base = adfdx.get_base() ;
	 
		adfdx = adfdx * dxdr * stdrdp ;  // df/dx / (dR/dx) * 1/sin(th) dR/dphi' 
	
		adfdx.set_base( sauve_base ) ; 

		// Final result 
		// ------------

		resu = stdfdp - adfdx ;

    }
	
    resu.set_dzpuis( ci.get_dzpuis() ) ; 	// dzpuis unchanged

}


			//------------------------------------//
			//       1/(r sin(theta))  d/dphi     //
			//------------------------------------//

void Map_star::srstdsdp(const Scalar& uu, Scalar& resu) const {

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

    // Computation of 1/(R sin(theta')) df/dphi'   ---> srstdfdp
    // -----------------------------------------
    
    Valeur srstdfdp = uuva ; 
    
    srstdfdp = srstdfdp.dsdp() ;	// d/dphi
    srstdfdp = srstdfdp.ssint() ;	// 1/sin(theta)
    srstdfdp = srstdfdp.sx() ;	// 1/xi, Id, 1/(xi-1)
    
    Base_val sauve_base( srstdfdp.base ) ; 
    
    srstdfdp = srstdfdp * xsr ;	// xi/R, 1/R, (xi-1)/U
    
    srstdfdp.base = sauve_base ;   // The above operation does not change the basis
    
    // Computation of 1/(dR/dx) 1/(R sin(theta') dR/dphi' df/dx   --> bdfdx
    // --------------------------------------------------------
    Valeur bdfdx = uuva ; 

    bdfdx = bdfdx.dsdx()  ;	    // df/dx 
    
    sauve_base = bdfdx.base ; 
    bdfdx = bdfdx * dxdr * srstdrdp  ;  
    bdfdx.base = sauve_base ; 


    if (uu.get_dzpuis() == 0) {

	//Final result

	resu = srstdfdp - bdfdx ;
 
 	  
      if (mg->get_type_r(nz-1) == UNSURR) {
	resu.set_dzpuis(2) ;	    // r d/dtheta has been computed in
	// the external domain
      }
    }

    else {
      assert(mg->get_type_r(nzm1) == UNSURR) ;
          
      int dzp = uu.get_dzpuis() ;

      Valeur tmp = srstdfdp - bdfdx ;
      tmp.annule(nzm1) ; 

      // Special treatment in the CED

      Valeur tmp_ced = - bdfdx ;
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED

      tmp_ced = tmp_ced.mult_x() ;	// xi, Id, (xi-1)
      //s Base_val sauve_base( tmp_ced.get_base() ) ; 
      tmp_ced = tmp_ced / xsr ; // xi/R, 1/R, (xi-1)/U
      
      tmp_ced = tmp_ced + uuva.dsdp().ssint() ;
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED

      // Recombination shells + CED : 
      resu = tmp + tmp_ced ;
      
      resu.set_dzpuis(dzp+1) ;         
    }
  }    
}


			//------------------------//
			//       d/dtheta         //
			//------------------------//


void Map_star::dsdt(const Scalar& ci, Scalar& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp().get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
		resu.set_etat_zero() ; 
    }
    else {

		assert( ci.get_etat() == ETATQCQ ) ; 

    // Computation of df/dtheta'   ---> dfdt
		// ----------------------------

		const Valeur& dfdt = ci.get_spectral_va().dsdt() ;
	

		// Computation of 1/(dR/dxi) dR/dtheta' df/dx   ----> adfdx
		// -------------------------------------------

		Valeur adfdx = ci.get_spectral_va().dsdx() ; 	// df/dx

		Base_val sauve_base = adfdx.get_base() ; 
	 
		adfdx = adfdx * dxdr * drdt ;  // df/dx / (dR/dx) * dR/dtheta' 
	
		adfdx.set_base( sauve_base ) ; 

		// Final result 
		// ------------

		resu = dfdt - adfdx ;
		
    }

    resu.set_dzpuis( ci.get_dzpuis() ) ; 	// dzpuis unchanged
	
}

void Map_star::srdsdt(const Scalar& uu, Scalar& resu) const {
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
    
    // Computation of 1/R df/dtheta'   ---> srdfdt
    // ----------------------------
    Valeur srdfdt = uuva ; 
	
    srdfdt = srdfdt.dsdt() ;	// d/dtheta'
   
    srdfdt = srdfdt.sx() ;		// 1/xi, Id, 1/(xi-1)
    	
    Base_val sauve_base( srdfdt.base ) ; 
	
    srdfdt = srdfdt * xsr ;	// xi/R, 1/R, (xi-1)/U
	
    srdfdt.base = sauve_base ;   // The above operation does not change the basis
    // Computation of 1/(dR/dx) 1/R dR/dtheta' df/dx   ----> adfdx
    // ----------------------------------------------

    Valeur adfdx = uuva ; 

    adfdx = adfdx.dsdx()  ;	    // df/dx 
		    
    sauve_base = adfdx.base ; 
    adfdx = adfdx * dxdr * srdrdt ;  // 1/(dR/dx) 1/R dR/dtheta' df/dx
    adfdx.base = sauve_base ; 

    if (uu.get_dzpuis() == 0) {

      // Final result 
      // ------------

      resu = srdfdt - adfdx ;

      //s int nz = mg->get_nzone() ; 
      if (mg->get_type_r(nz-1) == UNSURR) {
	resu.set_dzpuis(2) ;	    // r^2 (1/r d/dtheta) has been computed in
	// the external domain
      }
    
    }

    else {
      assert(mg->get_type_r(nzm1) == UNSURR) ;
          
      int dzp = uu.get_dzpuis() ;

      Valeur tmp = srdfdt - adfdx ;
      tmp.annule(nzm1) ; 
     
      // Special treatment in the CED
      //-----------------------------
 
      Valeur tmp_ced = - adfdx ;
  
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED

      tmp_ced = tmp_ced.mult_x() ;	// xi, Id, (xi-1)
      //s Base_val sauve_base( tmp_ced.get_base() ) ; 
      tmp_ced = tmp_ced / xsr ; // xi/R, 1/R, (xi-1)/U
        
      tmp_ced = tmp_ced + uuva.dsdt() ;
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED
                 
      // Recombination shells + CED : 
      resu = tmp + tmp_ced ;
                
      resu.set_dzpuis(dzp+1) ;         
    }
            
  }
    
}

}