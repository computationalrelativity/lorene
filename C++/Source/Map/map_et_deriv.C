/*
 * Computations of Cmp partial derivatives for a Map_et mapping
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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


char map_et_deriv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2004/04/08 17:16:07  f_limousin
 * Add comments
 *
 * Revision 1.5  2004/04/08 16:39:23  f_limousin
 * Add the case dzpuis different of 0 for methods dsdr, srdsdt, srstdsdp
 * for Scalar's.
 *
 * Revision 1.4  2004/01/26 16:16:17  j_novak
 * Methods of gradient for Scalar s. The input can have any dzpuis.
 *
 * Revision 1.3  2003/10/20 19:45:53  e_gourgoulhon
 * check_dzpuis in dsdt and stdsdp.
 *
 * Revision 1.2  2003/10/15 10:37:43  e_gourgoulhon
 * Added new methods dsdt and stdsdp.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/02/25  09:01:28  eric
 * Remplacement de ci.get_dzpuis() == 0  par ci.check_dzpuis(0).
 * Suppression de l'affectation des dzpuis Mtbl/Mtnl_cf a la fin car
 *   c'est fait par Cmp::set_dzpuis.
 *
 * Revision 1.2  2000/01/26  13:09:52  eric
 * Reprototypage complet des routines de derivation:
 * le resultat est desormais suppose alloue a l'exterieur de la routine
 * et est passe en argument (Cmp& resu), si bien que le prototypage
 * complet devient:
 *            void DERIV(const Cmp& ci, Cmp& resu) const
 *
 * Revision 1.1  1999/12/17  12:59:29  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Header Lorene
#include "cmp.h"
#include "tensor.h"


			//---------------------//
			//        d/dr         //
			//---------------------//
			
			
void Map_et::dsdr(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {    
	assert( ci.get_etat() == ETATQCQ ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date
	
	resu = (ci.va).dsdx() * dxdr ;     //  dxi/dR, - dxi/dU (ZEC)
	
	(resu.va).base = (ci.va).dsdx().base ;	// same basis as d/dxi
	
	int nz = mg->get_nzone() ; 
	if (mg->get_type_r(nz-1) == UNSURR) {
	    resu.set_dzpuis(2) ;	    // r^2 d/dr has been computed in the
					    // external domain
	}

    }
    
}

void Map_et::dsdr(const Scalar& uu, Scalar& resu) const {

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
      Valeur tmp_ced = uuva.dsdx() ;
      Base_val sauve_base( tmp_ced.get_base() ) ; 
      tmp_ced = tmp_ced * dxdr ;
      tmp_ced.set_base(sauve_base) ;   // The above operation does not 
                                       //change the basis
      tmp_ced = tmp_ced.mult_x() ;	// xi, Id, (xi-1)
      tmp_ced = tmp_ced / xsr ; // xi/R, 1/R, (xi-1)/U
   
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED
      tmp_ced.set(nzm1) -= dzp * uuva(nzm1) ;
      
      // Recombination shells + CED : 
      resu.set_spectral_va() += tmp_ced ;
      
      resu.set_dzpuis(dzp+1) ;         
    
    }

    resu.set_spectral_base( uuva.dsdx().get_base() ) ; // same basis as d/dxi

  }

}


			//------------------------//
			//      1/r d/dtheta      //
			//------------------------//

void Map_et::srdsdt(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	// Computation of 1/R df/dtheta'   ---> srdfdt
	// ----------------------------
	Valeur srdfdt = ci.va ; 
	
	srdfdt = srdfdt.dsdt() ;	// d/dtheta'
	srdfdt = srdfdt.sx() ;		// 1/xi, Id, 1/(xi-1)
	
	Base_val sauve_base( srdfdt.base ) ; 
	
	srdfdt = srdfdt * xsr ;	// xi/R, 1/R, (xi-1)/U
	
	srdfdt.base = sauve_base ;   // The above operation does not change the basis

	// Computation of 1/(dR/dx) 1/R dR/dtheta' df/dx   ----> adfdx
	// ----------------------------------------------

	Valeur adfdx = ci.va ; 

	adfdx = adfdx.dsdx()  ;	    // df/dx 
		    
	sauve_base = adfdx.base ; 
	adfdx = adfdx * dxdr * srdrdt ;  // 1/(dR/dx) 1/R dR/dtheta' df/dx
	adfdx.base = sauve_base ; 

	// Final result 
	// ------------

	resu = srdfdt - adfdx ;

	int nz = mg->get_nzone() ; 
	if (mg->get_type_r(nz-1) == UNSURR) {
	    resu.set_dzpuis(2) ;	    // r d/dtheta has been computed in
					    // the external domain
	}

    }
    
}

void Map_et::srdsdt(const Scalar& uu, Scalar& resu) const {

  assert (uu.get_etat() != ETATNONDEF) ; 
  assert (uu.get_mp().get_mg() == mg) ; 
    
  if (uu.get_etat() == ETATZERO) {
    resu.set_etat_zero() ; 
  }
  else {

    assert( uu.get_etat() == ETATQCQ ) ; 

    const Valeur& uuva = uu.get_spectral_va() ;
    uuva.coef() ;   // (uu.va).c_cf is up to date

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

    int nz = mg->get_nzone() ; 
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
      Base_val sauve_base( tmp_ced.get_base() ) ; 
      tmp_ced = tmp_ced / xsr ; // xi/R, 1/R, (xi-1)/U
        
      tmp_ced = tmp_ced + uuva.dsdt() ;
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED
                 
      // Recombination shells + CED : 
      resu = tmp + tmp_ced ;
                
      resu.set_dzpuis(dzp+1) ;         
    }
            
  }

}


			//------------------------------------//
			//       1/(r sin(theta))  d/dphi     //
			//------------------------------------//

void Map_et::srstdsdp(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 

    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	// Computation of 1/(R sin(theta')) df/dphi'   ---> srstdfdp
	// -----------------------------------------

	Valeur srstdfdp = ci.va ; 
	
	srstdfdp = srstdfdp.dsdp() ;	// d/dphi
	srstdfdp = srstdfdp.ssint() ;	// 1/sin(theta)
	srstdfdp = srstdfdp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
	Base_val sauve_base( srstdfdp.base ) ; 
	
	srstdfdp = srstdfdp * xsr ;	// xi/R, 1/R, (xi-1)/U
	
	srstdfdp.base = sauve_base ;   // The above operation does not change the basis

	// Computation of 1/(dR/dx) 1/(R sin(theta') dR/dphi' df/dx   --> bdfdx
	// --------------------------------------------------------
	Valeur bdfdx = ci.va ; 

	bdfdx = bdfdx.dsdx()  ;	    // df/dx 
		    
	sauve_base = bdfdx.base ; 
	bdfdx = bdfdx * dxdr * srstdrdp  ;  
	bdfdx.base = sauve_base ; 

	// Final result 
	// ------------

	resu = srstdfdp - bdfdx ;

	int nz = mg->get_nzone() ; 
	if (mg->get_type_r(nz-1) == UNSURR) {
	    resu.set_dzpuis(2) ;	    // r/sin(theta) d/dphi has been 
					    // computed in the external domain
	}

    }
    
}

void Map_et::srstdsdp(const Scalar& uu, Scalar& resu) const {

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
      Base_val sauve_base( tmp_ced.get_base() ) ; 
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
			//        d/dtheta        //
			//------------------------//

void Map_et::dsdt(const Scalar& ci, Scalar& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp().get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
		resu.set_etat_zero() ; 
    }
    else {

	//  The relations are true for all dzpuis
	//	assert( ci.check_dzpuis(0) ) ; 
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
    
}

			//---------------------------------//
			//        1/sin(theta) d/dphi      //
			//---------------------------------//

void Map_et::stdsdp(const Scalar& ci, Scalar& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp().get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
		resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ ) ; 
	//  The relations are true for all dzpuis
	//  assert( ci.check_dzpuis(0) ) ; 

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
    
}



