/*
 * Computations of Cmp partial derivatives for a Map_af mapping
 */

/*
 *   Copyright (c) 1999-2004 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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
 * Revision 1.14  2016/12/05 16:17:56  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.13  2014/10/13 08:53:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2012/01/17 10:32:09  j_penner
 * added a derivative with respect to the computational coordinate xi
 *
 * Revision 1.11  2008/09/21 13:57:21  j_novak
 * Changed the test on the CED in the derivative.
 *
 * Revision 1.10  2004/12/17 13:35:02  m_forot
 * Add the case T_LEG
 *
 * Revision 1.9  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.8  2004/01/27 09:33:48  j_novak
 * New method Map_radial::div_r_zec
 *
 * Revision 1.7  2004/01/26 16:16:17  j_novak
 * Methods of gradient for Scalar s. The input can have any dzpuis.
 *
 * Revision 1.6  2004/01/22 16:13:00  e_gourgoulhon
 * Case dzpuis=2 treated in dsdr, srdsdt and srstdsdp (output: dzpuis =
 * 3).
 * Reorganization cases dzpuis = 0 and 4.
 *
 * Revision 1.5  2003/11/11 15:31:43  j_novak
 * Added a #ifnedef... to prevent warnings.
 *
 * Revision 1.4  2003/10/22 13:08:05  j_novak
 * Better handling of dzpuis flags
 *
 * Revision 1.3  2003/10/20 19:45:27  e_gourgoulhon
 * Treatment of dzpuis in dsdt and stdsdp.
 *
 * Revision 1.2  2003/10/15 10:34:07  e_gourgoulhon
 * Added new methods dsdt and stdsdp.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.12  2000/02/25  08:59:51  eric
 * Remplacement de ci.get_dzpuis() == 0  par ci.check_dzpuis(0).
 * Suppression de l'affectation des dzpuis Mtbl/Mtnl_cf a la fin car
 *   c'est fait par Cmp::set_dzpuis.
 *
 * Revision 2.11  2000/01/26  13:09:18  eric
 * Reprototypage complet des routines de derivation:
 * le resultat est desormais suppose alloue a l'exterieur de la routine
 * et est passe en argument (Cmp& resu), si bien que le prototypage
 * complet devient:
 *             void DERIV(const Cmp& ci, Cmp& resu) const
 *
 * Revision 2.10  1999/11/30  12:51:32  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.9  1999/11/26  14:23:55  eric
 * Traitement dzpuis des Cmp.
 *
 * Revision 2.8  1999/11/26  10:58:02  eric
 * Traitement dzpuis.
 *
 * Revision 2.7  1999/11/25  16:29:29  eric
 * Reorganisation complete du calcul des derivees partielles.
 *
 * Revision 2.6  1999/10/27  15:45:23  eric
 * Suppression du membre Cmp::c.
 *
 * Revision 2.5  1999/10/27  08:47:03  eric
 * Introduction de Cmp::va a la place de *(Cmp::c).
 *
 * Revision 2.4  1999/10/22  08:16:21  eric
 * const Map*.
 *
 * Revision 2.3  1999/10/14  14:27:17  eric
 * Methodes const.
 *
 * Revision 2.2  1999/10/13  15:54:40  eric
 * Mg3d* -> const Mg3d*
 *
 * Revision 2.1  1999/09/17  10:01:09  phil
 * correction pour deriv_x et deriv_y
 *
 * Revision 2.0  1999/09/14  16:37:06  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
// Header Lorene
#include "map.h"
#include "cmp.h"
#include "tensor.h"


			//-----------------------//
			//        d/d\xi         //
			//-----------------------//
			

namespace Lorene {
void Map_af::dsdxi(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 

    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {   
	assert( ci.get_etat() == ETATQCQ ) ; 
        

	(ci.va).coef() ;    // (ci.va).c_cf is up to date
	
	int nz = mg->get_nzone() ; 
	int nzm1 = nz - 1 ;

        switch( ci.get_dzpuis() ) {
        
            case 0 : {
	        resu = (ci.va).dsdx() ;     // dsdx == d/d\xi
	
	        if (mg->get_type_r(nzm1) == UNSURR) {
	            resu.set_dzpuis(2) ;    // r^2 d/dr has been computed in the // SOMETHING IS WRONG HERE
					    // external domain
	        } 
                break ; 
            }
            
            case 2 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;
                
	        Valeur tmp((ci.va).dsdx() ) ;
                tmp.annule(nzm1) ;  // zero in the CED

                // Special treatment in the CED
                Valeur tmp_ced = - (ci.va).dsdx() ; 
                tmp_ced.annule(0, nz-2) ; // only non zero in the CED
                tmp_ced.mult_xm1_zec() ; 
                tmp_ced.set(nzm1) -= 2. * ci.va(nzm1) ; 
                
                // Recombination shells + CED : 
                resu = tmp + tmp_ced ;
                
                resu.set_dzpuis(3) ;         
                break ; 
            }
            
            case 4 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;

	        Valeur tmp(ci.va.dsdx()) ;
	        Valeur tmp2 = tmp ;
	        tmp2.base = (ci.va).dsdx().base ;
	        tmp.annule(nzm1) ; // not in the CED
	        tmp2.annule(0, nz-2) ; // special treatment of the CED
	        tmp2.mult_xm1_zec() ;
	        tmp2 = tmp2 / xsr ;
	        tmp2.set(nzm1) -= 4*ci.va(nzm1) ;
	        tmp2.base = ci.va.base ; //Just for the CED
	        tmp2.mult_xm1_zec() ;

	        resu = tmp + tmp2 / xsr  ; // do not know what this is, but for now I can get away with it, no CED
	        resu.set_dzpuis(4) ;
                break ; 
            }
            
            default : {
                cerr << "Map_af::dsdxi: unexpected value of input dzpuis !\n"
                    << "  ci.get_dzpuis() = " << ci.get_dzpuis() << endl ; 
                abort() ; 
                break ; 
            }
            
        }

	  
	(resu.va).set_base( (ci.va).dsdx().get_base() ) ; // same basis as d/dxi

    }
    
}

void Map_af::dsdxi(const Scalar& uu, Scalar& resu) const {

  assert (uu.get_etat() != ETATNONDEF) ; 
  assert (uu.get_mp().get_mg() == mg) ; 

  Mtbl unity = r/r;
    
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
      resu = uuva.dsdx() * unity ;     //  dxds == d/d\xi, unity is used to set the correct formated output
	
      if (mg->get_type_r(nzm1) == UNSURR) {
	resu.set_dzpuis(2) ;    // r^2 d/dr has been computed in the
	                        // external domain
      } 
    }
    else {
      int dzp = uu.get_dzpuis() ;
      
      resu = uuva.dsdx() ; 
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

			//---------------------//
			//        d/dr         //
			//---------------------//
			

void Map_af::dsdr(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 

    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {   
	assert( ci.get_etat() == ETATQCQ ) ; 
        

	(ci.va).coef() ;    // (ci.va).c_cf is up to date
	
	int nz = mg->get_nzone() ; 
	int nzm1 = nz - 1 ;

        switch( ci.get_dzpuis() ) {
        
            case 0 : {
	        resu = (ci.va).dsdx() * dxdr ;     //  dxdr = dxi/dR, - dxi/dU (ZEC)
	
	        if (mg->get_type_r(nzm1) == UNSURR) {
	            resu.set_dzpuis(2) ;    // r^2 d/dr has been computed in the
					    // external domain
	        } 
                break ; 
            }
            
            case 2 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;
                
	        Valeur tmp((ci.va).dsdx() * dxdr) ;
                tmp.annule(nzm1) ;  // zero in the CED

                // Special treatment in the CED
                Valeur tmp_ced = - (ci.va).dsdx() ; 
                tmp_ced.annule(0, nz-2) ; // only non zero in the CED
                tmp_ced.mult_xm1_zec() ; 
                tmp_ced.set(nzm1) -= 2. * ci.va(nzm1) ; 
                
                // Recombination shells + CED : 
                resu = tmp + tmp_ced ;
                
                resu.set_dzpuis(3) ;         
                break ; 
            }
            
            case 4 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;

	        Valeur tmp(ci.va.dsdx() * dxdr) ;
	        Valeur tmp2 = tmp ;
	        tmp2.base = (ci.va).dsdx().base ;
	        tmp.annule(nzm1) ; // not in the CED
	        tmp2.annule(0, nz-2) ; // special treatment of the CED
	        tmp2.mult_xm1_zec() ;
	        tmp2 = tmp2 / xsr ;
	        tmp2.set(nzm1) -= 4*ci.va(nzm1) ;
	        tmp2.base = ci.va.base ; //Just for the CED
	        tmp2.mult_xm1_zec() ;

	        resu = tmp + tmp2 / xsr  ; 
	        resu.set_dzpuis(4) ;
                break ; 
            }
            
            default : {
                cerr << "Map_af::dsdr: unexpected value of input dzpuis !\n"
                    << "  ci.get_dzpuis() = " << ci.get_dzpuis() << endl ; 
                abort() ; 
                break ; 
            }
            
        }

	  
	(resu.va).set_base( (ci.va).dsdx().get_base() ) ; // same basis as d/dxi

    }
    
}

void Map_af::dsdr(const Scalar& uu, Scalar& resu) const {

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
// VARIABLE RADIALE

void Map_af::dsdradial(const Scalar& uu, Scalar& resu) const {

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

			//------------------------//
			//      1/r d/dtheta      //
			//------------------------//

void Map_af::srdsdt(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	Valeur tmp = ci.va ; 
	
	tmp = tmp.dsdt() ;	// d/dtheta

	int nz = mg->get_nzone() ; 
	int nzm1 = nz - 1 ;

        switch( ci.get_dzpuis() ) {
        
            case 0 : {
	        tmp = tmp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
	        Base_val sauve_base( tmp.get_base() ) ; 

	        tmp = tmp * xsr ;	// xi/R, 1/R, (xi-1)/U

	        tmp.set_base(sauve_base) ;   // The above operation does not the basis
	        resu = tmp ;
	  
	        if (mg->get_type_r(nz-1) == UNSURR) {
	            resu.set_dzpuis(2) ;	    // r d/dtheta has been computed in
					            // the external domain
	        }
                break ; 
            }
            
            case 2 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;
                
                // Special treatment in the CED
                Valeur tmp_ced = tmp ;    // d/dtheta 
                tmp_ced.annule(0, nz-2) ; // only non zero in the CED

                tmp.annule(nzm1) ; 
	        tmp = tmp.sx() ;	// 1/xi, Id
	        Base_val sauve_base( tmp.get_base() ) ; 
	        tmp = tmp * xsr ;	// xi/R, 1/R
	        tmp.set_base(sauve_base) ;   // The above operation does not the basis
                
                // Recombination shells + CED : 
                resu = tmp + tmp_ced ;
                
                resu.set_dzpuis(3) ;         
                break ; 
            }
            
            case 4 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;

                // Special treatment in the CED
	        Valeur tmp_ced = tmp ;  // d/dtheta
                tmp_ced.annule(0, nz-2) ; // only non zero in the CED
	        tmp_ced.mult_xm1_zec() ;

                tmp.annule(nzm1) ; 
	        tmp = tmp.sx() ;	// 1/xi, Id
	        Base_val sauve_base( tmp.get_base() ) ; 
	        tmp = tmp * xsr ;	// xi/R, 1/R

                // Recombination shells + CED : 
	        resu = tmp + tmp_ced / xsr ;

	        resu.va.set_base( sauve_base ) ;
	        resu.set_dzpuis(4) ;
                break ; 
            }
            
            default : {
                cerr << "Map_af::srdsdt: unexpected value of input dzpuis !\n"
                    << "  ci.get_dzpuis() = " << ci.get_dzpuis() << endl ; 
                abort() ; 
                break ; 
            }
            
        }

    }
    
}

void Map_af::srdsdt(const Scalar& uu, Scalar& resu) const {
  assert (uu.get_etat() != ETATNONDEF) ; 
  assert (uu.get_mp().get_mg() == mg) ; 
    
  if (uu.get_etat() == ETATZERO) {
    resu.set_etat_zero() ; 
  }
  else {

    assert( uu.get_etat() == ETATQCQ ) ; 

    const Valeur& uuva = uu.get_spectral_va() ;
    uuva.coef() ;    // (uu.va).c_cf is up to date

    Valeur tmp  = uuva.dsdt() ;	// d/dtheta

    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ;

    if (uu.get_dzpuis() == 0) {
      tmp = tmp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
      Base_val sauve_base( tmp.get_base() ) ; 

      tmp = tmp * xsr ;	// xi/R, 1/R, (xi-1)/U

      tmp.set_base(sauve_base) ;   // The above operation does not change the basis
      resu = tmp ;
	  
      if (mg->get_type_r(nz-1) == UNSURR) {
	resu.set_dzpuis(2) ;	    // r d/dtheta has been computed in
	// the external domain
      }
    }

    else {
      assert(mg->get_type_r(nzm1) == UNSURR) ;
          
      int dzp = uu.get_dzpuis() ;
      // Special treatment in the CED
      Valeur tmp_ced = tmp ;    // d/dtheta 
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED

      tmp.annule(nzm1) ; 
      tmp = tmp.sx() ;	// 1/xi, Id
      Base_val sauve_base( tmp.get_base() ) ; 
      tmp = tmp * xsr ;	// xi/R, 1/R
      tmp.set_base(sauve_base) ;   // The above operation does not change the basis
                
      // Recombination shells + CED : 
      resu = tmp + tmp_ced ;
                
      resu.set_dzpuis(dzp+1) ;         
    }
            
  }
    
}


			//------------------------------------//
			//       1/(r sin(theta))  d/dphi     //
			//------------------------------------//

void Map_af::srstdsdp(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	Valeur tmp = ci.va ; 
	
	tmp = tmp.dsdp() ;	// d/dphi
	tmp = tmp.ssint() ;	// 1/sin(theta)

	int nz = mg->get_nzone() ; 
	int nzm1 = nz - 1 ;

        switch( ci.get_dzpuis() ) {
        
            case 0 : {
	        tmp = tmp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
	        Base_val sauve_base( tmp.get_base() ) ; 

	        tmp = tmp * xsr ;	// xi/R, 1/R, (xi-1)/U

	        tmp.set_base(sauve_base) ;   // The above operation does not the basis
	        resu = tmp ;
	  
	        if (mg->get_type_r(nz-1) == UNSURR) {
	            resu.set_dzpuis(2) ;	    // r d/dtheta has been computed in
					            // the external domain
	        }
                break ; 
            }
            
            case 2 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;
                
                // Special treatment in the CED
                Valeur tmp_ced = tmp ;    // 1/sin(theta) d/dphi 
                tmp_ced.annule(0, nz-2) ; // only non zero in the CED

                tmp.annule(nzm1) ; 
	        tmp = tmp.sx() ;	// 1/xi, Id
	        Base_val sauve_base( tmp.get_base() ) ; 
	        tmp = tmp * xsr ;	// xi/R, 1/R
	        tmp.set_base(sauve_base) ;   // The above operation does not the basis
                
                // Recombination shells + CED : 
                resu = tmp + tmp_ced ;
                
                resu.set_dzpuis(3) ;         
                break ; 
            }
            
            case 4 : {
	        assert(mg->get_type_r(nzm1) == UNSURR) ;

                // Special treatment in the CED
	        Valeur tmp_ced = tmp ;  // 1/sin(theta) d/dphi 
                tmp_ced.annule(0, nz-2) ; // only non zero in the CED
	        tmp_ced.mult_xm1_zec() ;

                tmp.annule(nzm1) ; 
	        tmp = tmp.sx() ;	// 1/xi, Id
	        Base_val sauve_base( tmp.get_base() ) ; 
	        tmp = tmp * xsr ;	// xi/R, 1/R

                // Recombination shells + CED : 
	        resu = tmp + tmp_ced / xsr ;

	        resu.va.set_base( sauve_base ) ;
	        resu.set_dzpuis(4) ;
                break ; 
            }
            
            default : {
                cerr << "Map_af::srstdsdp: unexpected value of input dzpuis !\n"
                    << "  ci.get_dzpuis() = " << ci.get_dzpuis() << endl ; 
                abort() ; 
                break ; 
            }
            
        }

    }

    
}

void Map_af::srstdsdp(const Scalar& uu, Scalar& resu) const {

  assert (uu.get_etat() != ETATNONDEF) ; 
  assert (uu.get_mp().get_mg() == mg) ; 
    
  if (uu.get_etat() == ETATZERO) {
    resu.set_etat_zero() ; 
  }
  else {

    assert( uu.get_etat() == ETATQCQ ) ; 

    const Valeur& uuva = uu.get_spectral_va() ;
    uuva.coef() ;    // (uu.va).c_cf is up to date

    Valeur tmp = uuva.dsdp() ;
	
    tmp = tmp.ssint() ;	// 1/sin(theta)

    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ;

    if (uu.get_dzpuis() == 0) {

      tmp = tmp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
      Base_val sauve_base( tmp.get_base() ) ; 

      tmp = tmp * xsr ;	// xi/R, 1/R, (xi-1)/U

      tmp.set_base(sauve_base) ;   // The above operation does not change the basis
      resu = tmp ;
	  
      if (mg->get_type_r(nz-1) == UNSURR) {
	resu.set_dzpuis(2) ;	    // r d/dtheta has been computed in
	// the external domain
      }
    }

    else {
      assert(mg->get_type_r(nzm1) == UNSURR) ;
          
      int dzp = uu.get_dzpuis() ;

      // Special treatment in the CED
      Valeur tmp_ced = tmp ;    // 1/sin(theta) d/dphi 
      tmp_ced.annule(0, nz-2) ; // only non zero in the CED

      tmp.annule(nzm1) ; 
      tmp = tmp.sx() ;	// 1/xi, Id
      Base_val sauve_base( tmp.get_base() ) ; 
      tmp = tmp * xsr ;	// xi/R, 1/R
      tmp.set_base(sauve_base) ;   // The above operation does not change the basis
                
      // Recombination shells + CED : 
      resu = tmp + tmp_ced ;
                
      resu.set_dzpuis(dzp+1) ;         
    }
  }    
}


			//------------------------//
			//       d/dtheta         //
			//------------------------//


void Map_af::dsdt(const Scalar& ci, Scalar& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp().get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
		resu.set_etat_zero() ; 
    }
    else {

		assert( ci.get_etat() == ETATQCQ ) ; 

		resu = ci.get_spectral_va().dsdt() ; 	// d/dtheta
		
    }

    resu.set_dzpuis( ci.get_dzpuis() ) ; 	// dzpuis unchanged
	
}


			//-----------------------------------//
			//       1/sin(theta) d/dphi         //
			//-----------------------------------//


void Map_af::stdsdp(const Scalar& ci, Scalar& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp().get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
		resu.set_etat_zero() ; 
    }
    else {

		assert( ci.get_etat() == ETATQCQ ) ; 
	
		resu = ci.get_spectral_va().stdsdp() ; 	// 1/sin(theta) d/dphi
		
    }
	
    resu.set_dzpuis( ci.get_dzpuis() ) ; 	// dzpuis unchanged

}








}
