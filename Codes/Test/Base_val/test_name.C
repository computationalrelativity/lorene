/*
 * Code for testing Base_val::name_* methods
 *
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon. 
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

char test_name_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/10/19 20:03:28  e_gourgoulhon
 * First version
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>

// Lorene headers
#include "base_val.h"

int main() {

	int nz = 3 ; 
	Base_val base(nz) ; 
	
	base.set_base_r(0, R_CHEBU) ; 
	
	base.set_base_t(T_LEG_II) ; 

	base.set_base_p(P_COSSIN_I) ; 
	
	char name[8] ; 
	
	int np = 6 ;
	int nt = 13 ; 
	int nr = 17 ;
	
	for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
			base.name_theta(0, k, j, name) ; 
			cout << "k=" << k << ", j=" << j << " : " << name << endl ; 
		}
	}

	cout << endl ; 
	for (int k=0; k<np+1; k++) {
		for (int i=0; i<nr; i++) {
			base.name_r(0, k, 0, i, name) ; 
			cout << "k=" << k << ", i=" << i << " : " << name << endl ; 
		}
	}


	return EXIT_SUCCESS ; 
}
