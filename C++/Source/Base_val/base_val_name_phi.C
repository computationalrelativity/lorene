/*
 *  Method Base_val::name_phi
 *
 *	(see file base_val.h for documentation). 
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

char base_val_name_phi_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2012/01/17 14:44:27  j_penner
 * Modified phi variables to only use 16 integers in arrays
 *
 * Revision 1.1  2003/10/19 19:49:40  e_gourgoulhon
 * First version
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <string.h>
#include <stdlib.h>

// Lorene headers
#include "base_val.h"

// Local prototypes
void basename_p_unknown(int, char*) ; 
void basename_p_cossin(int, char*) ; 
void basename_p_cossin_p(int, char*) ; 
void basename_p_cossin_i(int, char*) ; 

			//----------------------------//
			//      Base_val method       //
			//----------------------------//

void Base_val::name_phi(int l, int k, char* name) const {

	// Array of actual base name functions
    static void(*vbasename_p[MAX_BASE_2])(int, char*) ;  

    static bool first_call = true ;

    // Initializations at first call
    // -----------------------------
    if ( first_call ) {

		first_call = false ;

		for (int i=0 ; i<MAX_BASE_2 ; i++) {
	    	vbasename_p[i] = basename_p_unknown ;
		}

		vbasename_p[P_COSSIN >> TRA_P] = basename_p_cossin ;
		vbasename_p[P_COSSIN_P >> TRA_P] = basename_p_cossin_p ;
		vbasename_p[P_COSSIN_I >> TRA_P] = basename_p_cossin_i ;

    }
	
	// Call to the function adapted to the basis in domain l
	//------------------------------------------------------
	
	assert( (l>=0) && (l<nzone) ) ; 

    int base_p = ( b[l] & MSQ_P ) >> TRA_P ;
	
	vbasename_p[base_p](k, name) ; 

}
	
	
			//-------------------------------//
            //  individual basis functions   //
			//-------------------------------//
	
void basename_p_unknown(int, char*) {
	cout << "Base_val::name_phi : unknwon basis !" << endl ; 
	abort() ; 
} 


void basename_p_cossin(int k, char* name) {

	assert( k>=0 ) ; 

	if (k%2 == 0) {
		strcpy(name, "cos") ; 
	}		
	else {
		if (k == 1) {
			strcpy(name, "unused") ; 
			return ;
		}
		else { 
			strcpy(name, "sin") ; 
		}
	}
		
	int m = k / 2 ; 
		
	char cm[4] ;
	assert( m < 1000) ; 
	sprintf(cm, "%d", m) ; 
	strcat(name, cm) ; 
	strcat(name, "p") ; 

}	

	
	
void basename_p_cossin_p(int k, char* name) {

	assert( k>=0 ) ; 

	if (k%2 == 0) {
		strcpy(name, "cos") ; 
	}		
	else {
		if (k == 1) {
			strcpy(name, "unused") ; 
			return ;
		}
		else { 
			strcpy(name, "sin") ; 
		}
	}
		
	int m = 2 * (k / 2) ; 
		
	char cm[4] ;
	assert( m < 1000) ; 
	sprintf(cm, "%d", m) ; 
	strcat(name, cm) ; 
	strcat(name, "p") ; 

} 


void basename_p_cossin_i(int k, char* name) {

	assert( k>=0 ) ; 

	if (k == 0) {
		strcpy(name, "cos1p") ; 
		return  ; 		
	}

	if (k%2 == 0) {
		strcpy(name, "sin") ; 
	}		
	else {
		if (k == 1) {
			strcpy(name, "unused") ; 
			return ;
		}
		else { 
			strcpy(name, "cos") ; 
		}
	}
		
	int m = 2 * ((k-1) / 2) + 1 ; 
		
	char cm[4] ;
	assert( m < 1000) ; 
	sprintf(cm, "%d", m) ; 
	strcat(name, cm) ; 
	strcat(name, "p") ; 

} 


	
	
	
	
	
	
	
	
	
	
	
	
	
