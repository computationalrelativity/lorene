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


 

/*
 *  Calcul de l'integrale
 * 
 *	    int_-1^1 f(x) dx					    (1)
 *
 *  pour une fonction f(x)  donnee par ses coefficients de Tchebyshev
 *
 *	    f(x) = som_{i=0}^{nr-1} c_i T_{i}(x)		    (2)
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
 *  double int1d_cheb	:   Valeur de l'integrale (1) 
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2016/12/05 16:18:07  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:23  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2005/02/16 15:27:55  m_forot
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

namespace Lorene {

//*****************************************************************************

double int1d_cheb(int nr, const double* cf){
    
    double som = - cf[0] ;
    const double* cc = cf + 2 ;
    
    for (int i=2; i<nr ; i=i+2) {
      som += *cc / (i*i - 1) ;
      cc = cc + 2 ;
    }
    
    return -2*som ; 
    
}
}
