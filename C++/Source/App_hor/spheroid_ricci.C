/*
 *  Method for calculating the Ricci tensor for a Spheroid.
 *
 *    (see file spheroid.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006  Jose-Luis Jaramillo, Jerome Novak & Nicolas Vasset
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

char spheroid_ricci_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2006/05/26 13:49:10  j_novak
 * Method for computing the Ricci tensor on the Spheroid.
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "spheroid.h"

const Sym_tensor& Spheroid::ricci() const {

    if (p_ricci == 0x0) {
	const Map& map_2d = h_surf.get_mp() ;
	p_ricci = 
	    new Sym_tensor(map_2d, COV, map_2d.get_bvect_spher()) ;

	for (int i=1; i<=3; i++)
	    p_ricci->set(1,i) = 0 ;

	Vector ei(map_2d, CON, map_2d.get_bvect_spher()) ;
	ei.set_etat_zero() ;
	ei.set(2) = 1 ;
	ei.set(2).std_spectral_base() ;
	ei.set(2).mult_cost() ;
	ei.set(2).mult_sint() ;

	Tensor dei = ei.derive_cov(qab) ;
	dei = contract(proj, 1, contract(proj, 0, dei, 1), 1) ;
	Tensor d2ei = dei.derive_cov(qab) ;
	d2ei = contract(proj, 1, contract(proj, 0, contract(proj, 0, d2ei, 2), 2), 2) ;
	Vector r1 = d2ei.trace(0,2) - d2ei.trace(0,1) ;
	for (int i=1; i<=3; i++) {
	    r1.set(i).div_sint() ;
	    r1.set(i).div_cost() ;
	}
	for (int i=1; i<=3; i++)
	    p_ricci->set(2,i) = r1(i) ;
	
	ei.set_etat_zero() ;
	ei.set(3) = 1 ;
	ei.set(3).std_spectral_base() ;
	ei.set(3).mult_sint() ;

	dei = ei.derive_cov(qab) ;
	dei = contract(proj, 1, contract(proj, 0, dei, 1), 1) ;
	d2ei = dei.derive_cov(qab) ;
	d2ei = contract(proj, 1, contract(proj, 0, contract(proj, 0, d2ei, 2), 2), 2) ;
	Vector r2 = d2ei.trace(0,2) - d2ei.trace(0,1) ;
	for (int i=1; i<=3; i++) {
	    r2.set(i).div_sint() ;
	}
	for (int i=1; i<=3; i++)
	    p_ricci->set(3,i) = r2(i) ;
    }

    return *p_ricci ;
}
 
 
