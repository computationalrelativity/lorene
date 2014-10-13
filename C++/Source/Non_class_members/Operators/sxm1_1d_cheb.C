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


char sxm1_1d_cheb_C[] = "$Header$" ;

/*
 * Operateur [f(x) - f(1)]/(x-1) applique a une fonction f(x) developpee en 
 * polynomes de Tchebychev (echantillonnage fin: x ds. [-1, 1]) :
 *
 *	    f(x) = som_{i=0}^{nr-1} c_i T_i(x)		    (1)
 *
 *
 * Entree:
 * ------
 *  int nr  :		    Nombre de coefficients de Tchebyshev dans le 
 *			    developpement (2)
 *
 * Entree/Sortie :
 * -------------
 *  double* cf	:   entree: Tableau des nr coefficients c_i de la fonction f(x)
 *			    definis par (1). Le stokage doit etre le suivant
 *				cf[i] = c_i   0 <= i <= nr - 1
 *			    L'espace memoire correspondant au pointeur cf doit
 *			    etre de taille au moins nr et doit avoir ete 
 *			    alloue avant l'appel a la routine
 *		    sortie: Tableau des nr coefficients c_i de la fonction 
 *			    f(x)/(x-1). On a cf[nr-1] = 0. 
 *
 */


/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:27  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/04/26  14:59:56  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

namespace Lorene {


//*****************************************************************************

void sxm1_1d_cheb(int nr, double* cf) {

//-------------------------------------------------------
// Formulation effectuant f(x)-f(1)/(x-1) (le coef. c_0 n'intervient donc pas)
//-------------------------------------------------------

    int i, j ; 

// Coefficient i=0 du resultat :

    double som = cf[1] ;
    for (j=2; j<nr; j++) {
	som += j * cf[j] ;
    }
    cf[0] = som ;
    
// Coefficients 1 <= i <= nr-2 du resultat :

    for (i=1; i<nr-1; i++) {
	som = cf[i+1] ;
	for (j=i+2; j<nr; j++) {
	    som += (j-i) * cf[j] ;
	}
	cf[i] = 2 * som ;
    }
    
// Coefficient i=nr-1 du resultat :
    cf[nr-1] = 0 ;
  

/*
//-------------------------------------------------------
// Formulation privilegiant c_{0} au detriment de c_{N-1}
//-------------------------------------------------------


    // Coefficient i=0 du resultat :
    // ---------------------------
    int nrm1 = nr - 1 ; 
    double som = nrm1*cf[0] ;
    double cfim1 = cf[0] ;	// pour ne pas perdre le coef c_0 de l'entree
    int i ; 
    for (i=1; i<nrm1; i++) {
	som += (nrm1-i)*cf[i] ;
    } 
    cf[0] = - som ; 
    
    // Coefficient i=1 du resultat :
    // ---------------------------
    som = cfim1 ;	// coef c_0 de l'entree
    cfim1 = cf[1] ;	// coef c_1 de l'entree 
    cf[1] = 2 * (cf[0] + som) ;  
    som += cfim1 ;	// a ce stade som = c_0 + c_1 de l'entree

    // Coefficients 2 <= i <= nr-2 du resultat :
    // ----------============-----------------
    for (i=2; i<nrm1; i++) {
	cfim1 = cf[i] ;		    // coef c_i de l'entree
	cf[i] = cf[i-1] + 2*som ; 
	som += cfim1 ;		    // som = c_0 + c_1 + ... + c_i de l'entree
    }
    
    // Coefficient i=nr-1 du resultat :
    // ------------------------------
    cf[nrm1] = 0 ;

*/

}



}
