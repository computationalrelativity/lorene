/*
 * Computes the Laplacian of a scalar field (represented by a Cmp) when
 *  the mapping belongs to the Map_et class
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char map_et_lap_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/15 16:03:37  j_novak
 * Added the angular Laplace operator for Scalar.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.4  2000/02/25  09:02:18  eric
 * Remplacement de (uu.get_dzpuis() == 0) par uu.check_dzpuis(0).
 *
 * Revision 1.3  2000/01/26  13:11:07  eric
 * Reprototypage complet :
 * le resultat est desormais suppose alloue a l'exterieur de la routine
 * et est passe en argument (Cmp& resu).
 *  .
 *
 * Revision 1.2  2000/01/14  14:55:05  eric
 * Suppression de l'assert(sauve_base == vresu.base)
 *  car sauve_base == vresu.base n'est pas necessairement vrai (cela
 *  depend de l'histoire du Cmp uu, notamment de si uu.p_dsdx est
 *  a jour).
 *
 * Revision 1.1  1999/12/20  17:27:30  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header Lorene
#include "map.h"
#include "cmp.h" 
#include "tensor.h"
#include "type_parite.h" 

void Map_et::laplacien(const Cmp& uu, int zec_mult_r, Cmp& resu) const {
 
    assert (uu.get_etat() != ETATNONDEF) ; 
    assert (uu.get_mp()->get_mg() == mg) ; 

    // Particular case of null value:

    if (uu.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
	return ; 
    }

    assert( uu.get_etat() == ETATQCQ ) ; 
    assert( uu.check_dzpuis(0) ) ; 
    
    int nz = get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ;        

    // Indicator of existence of a compactified external domain
    bool zec = false ; 		
    if (mg->get_type_r(nzm1) == UNSURR) {
	zec = true ;
    }
     
    if ( zec && (zec_mult_r != 4) ) {
	cout << "Map_et::laplacien : the case zec_mult_r = " << 
		zec_mult_r << " is not implemented !" << endl ; 
	abort() ; 
    }
 
    //--------------------
    // First operations
    //--------------------
    
    Valeur duudx = (uu.va).dsdx() ;	    // d/dx 
    Valeur d2uudx2 = (uu.va).d2sdx2() ;	    // d^2/dx^2

    const Valeur& d2uudtdx = duudx.dsdt() ;	    // d^2/dxdtheta

    const Valeur& std2uudpdx = duudx.stdsdp() ;    // 1/sin(theta) d^2/dxdphi   

    //------------------
    // Angular Laplacian 
    //------------------
    
    Valeur sxlapang = uu.va ; 

    sxlapang.ylm() ; 
    
    sxlapang = sxlapang.lapang() ;    
    
    sxlapang = sxlapang.sx() ;	  // 1/x    in the nucleus
				  // Id	    in the shells
				  // 1/(x-1) in the ZEC
    
    //------------------------------------
    // (2 dx/dR d/dx + x/R 1/x Lap_ang)/x
    //------------------------------------

    Valeur varduudx = duudx ; 

    if (zec) {
	varduudx.annule(nzm1) ;	    // term in d/dx set to zero in the ZEC
    }

    resu.set_etat_qcq() ; 
    
    Valeur& vresu = resu.va ; 

    Base_val sauve_base = varduudx.base ; 
    
    vresu = double(2) * dxdr * varduudx + xsr * sxlapang ; 

    vresu.set_base(sauve_base) ; 

    vresu = vresu.sx() ; 

    //--------------
    // Final result
    //--------------

    Mtbl unjj = double(1) + srdrdt*srdrdt + srstdrdp*srstdrdp ;

    sauve_base = d2uudx2.base ; 
//    assert(sauve_base == vresu.base) ;  // this is not necessary true
    
    vresu = dxdr*dxdr * unjj * d2uudx2 + xsr * vresu 
		- double(2) * dxdr * (	sr2drdt * d2uudtdx 
			      + sr2stdrdp * std2uudpdx ) ;

    vresu += - dxdr * ( lapr_tp + dxdr * ( 
		dxdr* unjj * d2rdx2 
		- double(2) * ( sr2drdt * d2rdtdx  + sr2stdrdp * sstd2rdpdx ) ) 
				 ) * duudx ;		    

    vresu.set_base(sauve_base) ; 

    if (zec == 1) {
	resu.set_dzpuis(zec_mult_r) ; 
    }

}

void Map_et::lapang(const Scalar& , Scalar& ) const {
 
  cout << "Map_et::lapang : not implemented yet!" << endl ;
  abort() ;

}


