/*
 * Determines whether a given number of points N is allowed by the
 *   Fast Fourier Transform algorithm, i.e. if
 *	
 *	    N = 2^p 3^q 5^r  and N >= 4, p>=1
 *
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


char admissible_fft_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/12/21 17:06:01  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  1999/11/24  16:06:52  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */
 
bool admissible_fft(int n) {
     
    if (n < 4) {
	return false ; 
    }

     // Division by 2
     //--------------
     
    int reste = n % 2 ; 
    if (reste != 0) {
	return false ; 
    }
     
    int k = n/2 ; 
     
    while ( k % 2 == 0 ) {
	k = k / 2 ;  
    }

    if (k == 1) return true ;	    // n = 2^p 

    // Division by 3
    //--------------
     
    while ( k % 3 == 0 ) {
	k = k / 3 ;  
    }
     
    if (k == 1) return true ;	    // n = 2^p * 3^q 
     
    // Division by 5
    //--------------
     
    while ( k % 5 == 0 ) {
	k = k / 5 ;  
    }
     
    if (k == 1) return true ;	    // n = 2^p * 3^q * 5^r 

    return false ; 

 }
