/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char mult2_xm1_1d_cheb_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.1  1999/10/11  14:27:55  phil
 * vire (double(0.5))
 *
 * Revision 2.0  1999/04/26  16:28:15  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
 
/*
 * Operateur (x-1)^2 Id applique a une fonction f(x) developpee en 
 * polynomes de Tchebychev (echantillonnage fin: x ds. [-1, 1]) :
 *
 *	    f(x) = som_{i=0}^{nr-1} c_i T_i(x)		    (1)
 *
 *
 * Entree:
 * ------
 *  int nr  :		    Nombre de coefficients de Tchebyshev dans le 
 *			    developpement (1)
 *
 *  const double* cf	:   Tableau des nr coefficients c_i de la fonction f(x)
 *			    definis par (1). Le stokage doit etre le suivant
 *				cf[i] = c_i   0 <= i <= nr - 1
 *			    L'espace memoire correspondant au pointeur cf doit
 *			    etre de taille au moins nr et doit avoir ete 
 *			    alloue avant l'appel a la routine.
 * Sortie :
 * -------
 *  double* cresu	:  Tableau des nr coefficients de la fonction 
 *			    (x-1)^2 f(x).
 *			   L'espace memoire correspondant au pointeur cresu doit
 *			   etre de taille au moins nr et doit avoir ete 
 *			   alloue avant l'appel a la routine.
 *
 */
 
 #include <assert.h>

//*****************************************************************************

void mult2_xm1_1d_cheb(int nr, const double* cf,  double* cresu) {

    double aim2 = 0.25 ; 
    double aim1 = -1. ; 
    double ai	= 1.5 ; 
    double aip1 = -1. ; 
    double aip2 = 0.25 ; 
    
    assert(nr>=3) ; 
    
// Coefficient i=0 du resultat :

    cresu[0] = ai*cf[0] + aip1*cf[1] + aip2*cf[2] ; 
    
// Coefficient i=1 du resultat :

    cresu[1] = double(-2)*cf[0] + double(1.75)*cf[1] + aip1*cf[2] + aip2*cf[3] ; 
         
// Coefficient i=2 du resultat :

    cresu[2] = double(0.5)*cf[0] + aim1*cf[1] + ai*cf[2] + aip1*cf[3]
		    + aip2*cf[4] ; 
    
// Coefficients 3 <= i <= nr-3 du resultat :

    int i ; 
    for (i=3; i<nr-2; i++) {
	cresu[i] = aim2*cf[i-2] + aim1*cf[i-1] + ai*cf[i] + aip1*cf[i+1]
		    + aip2*cf[i+2] ; 
    }
    
// Coefficient i=nr-2 du resultat :

    cresu[nr-2] = aim2*cf[nr-4] + aim1*cf[nr-3] + ai*cf[nr-2] + aip1*cf[nr-1] ;
    
// Coefficient i=nr-1 du resultat :  
 
    cresu[nr-1] = aim2*cf[nr-3] + aim1*cf[nr-2] + ai*cf[nr-1] ;
     
}



