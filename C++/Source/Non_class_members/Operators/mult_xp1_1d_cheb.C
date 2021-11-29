/*
 *   Copyright (c) 2021 Gaël Servignat
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
 * Operateur (x+1) Id applique a une fonction f(x) developpee en 
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
 *			    (x+1) f(x).
 *			   L'espace memoire correspondant au pointeur cresu doit
 *			   etre de taille au moins nr et doit avoir ete 
 *			   alloue avant l'appel a la routine.
 *
 */
 
 
 #include <cassert>

namespace Lorene {

//*****************************************************************************

void mult_xp1_1d_cheb(int nr, const double* cf,  double* cresu) {

    double aim1 = 0.5 ; 
    double ai	= 1. ; 
    double aip1 = 0.5 ; 
    
    assert(nr>=3) ; 
    
// Coefficient i=0 du resultat :

    cresu[0] = ai*cf[0] + aip1*cf[1] ; 
    
// Coefficient i=1 du resultat :

    cresu[1] = cf[0] + ai*cf[1] + aip1*cf[2] ; 
         
    
// Coefficients 2 <= i <= nr-2 du resultat :

    int i ; 
    for (i=2; i<nr-1; i++) {
	cresu[i] = aim1*cf[i-1] + ai*cf[i] + aip1*cf[i+1] ; 
    }
    
// Coefficient i=nr-1 du resultat :  
 
    cresu[nr-1] = aim1*cf[nr-2] + ai*cf[nr-1] ;
     
}



}