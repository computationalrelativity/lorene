/*
 *  Write binary data into a file according to the Big Endian convention
 *
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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

char fwrite_be_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/12/04 21:32:39  e_gourgoulhon
 * Functions similar to the stdio fread/fwrite except that they ensure
 * the big endian convention, whatever the system convention is.
 *
 * Revision 1.1  2001/11/23 15:09:09  e_gourgoulhon
 * Templates for new source files
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdio.h>
#include <assert.h>

			//-------------------------//
			//	int version 	   //
			//-------------------------//
			

int fwrite_be(const int* aa, int size, int nb, FILE* fich) {

	assert(size == 4) ;
	
	// Determines whether the default storage is big endian
	//  or large endians
	
	int itest = 1 ;
	bool little_endian = ( *( (char*) &itest ) == 1) ;
	
	if (little_endian) {

		int size_tot = 4 * nb ;

		char bytes_big[size_tot] ;
		char* pbig =  bytes_big ;
		char* plit = (char*) aa ;
		
		for (int j=0; j< nb; j++) {
		
			for (int i=0; i<4; i++) {
				pbig[i] = plit[3-i] ;
			}
		
			plit += 4 ; 	// next item
			pbig += 4 ;
			
		}
		
		return fwrite(bytes_big, 1, size_tot, fich) / 4 ;		
		
	}
	else {  // Big endian case: nothing to do:
	
		return fwrite(aa, size, nb, fich) ;
	}
		
}


			//-------------------------//
			//	double version 	   //
			//-------------------------//
			

int fwrite_be(const double* aa, int size, int nb, FILE* fich) {

	assert(size == 8) ;
	
	// Determines whether the default storage is big endian
	//  or large endians
	
	int itest = 1 ;
	bool little_endian = ( *( (char*) &itest ) == 1) ;
	
	if (little_endian) {

		int size_tot = 8 * nb ;

		char bytes_big[size_tot] ;
		char* pbig =  bytes_big ;
		char* plit = (char*) aa ;
		
		for (int j=0; j< nb; j++) {
		
			for (int i=0; i<8; i++) {
				pbig[i] = plit[7-i] ;
			}
		
			plit += 8 ; 	// next item
			pbig += 8 ;
			
		}
		
		return fwrite(bytes_big, 1, size_tot, fich) / 8 ;		
		
	}
	else {  // Big endian case: nothing to do:
	
		return fwrite(aa, size, nb, fich) ;
	}
		
}
