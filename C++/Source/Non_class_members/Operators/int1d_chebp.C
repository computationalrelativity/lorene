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


char int1d_chebp_C[] = "$Header$" ;

/*
 *  Calcul de l'integrale
 * 
 *	    int_0^1 f(x) dx					    (1)
 *
 *  pour une fonction f(x) paire donnee par ses coefficients de Tchebyshev
 *
 *	    f(x) = som_{i=0}^{nr-1} c_i T_{2i}(x)		    (2)
 *
 * Entree:
 * ------
 *  int nr  :		    Nombre de coefficients de Tchebyshev dans le 
 *			    developpement (2)
 *  const double* cf	:   Tableau des nr coefficients c_i de la fonction
 *			    definis par (2). Le stokage doit etre le suivant
 *				cf[i] = c_i   0 <= i <= nr - 1
 *			    L'espace memoire correspondant au pointeur cf doit
 *			    etre de taille au moins nr et doit avoir ete 
 *			    alloue avant l'appel a la routine
 *
 * Sortie (valeur de retour) :
 * ------
 *  double int1d_chebp	:   Valeur de l'integrale (1) 
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  1999/02/22  15:37:46  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


//*****************************************************************************

double int1d_chebp(int nr, const double* cf){
    
    double som = - cf[0] ;
    const double* cc = cf + 1 ;
    
    for (int i=1; i<nr ; i++) {
	som += *cc / (4*i*i - 1) ; 
	cc++ ;
    }
    
    return -som ; 
    
}
