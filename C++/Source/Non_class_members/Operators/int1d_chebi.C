/*
 *   Copyright (c) 2005 Jerome Novak
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


char int1d_chebi_C[] = "$Header$" ;

/*
 *  Calcul de l'integrale
 * 
 *	    int_0^1 f(x) dx					    (1)
 *
 *  pour une fonction f(x) paire donnee par ses coefficients de Tchebyshev
 *
 *	    f(x) = som_{i=0}^{nr-1} c_i T_{2i+1}(x)		    (2)
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
 *  double int1d_chebi	:   Valeur de l'integrale (1) 
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:23  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2005/11/02 15:08:18  j_novak
 * Minor change to prevent warning message.
 *
 * Revision 1.2  2005/05/13 13:22:33  j_novak
 * *** empty log message ***
 *
 * Revision 1.1  2005/05/13 08:51:02  j_novak
 * Added the function int1d_chebi.
 *
 *
 * $Header$
 *
 */

namespace Lorene {

//*****************************************************************************

double int1d_chebi(int nr, const double* cf){
    
    double som = 0. ;
    const double* cc = cf ;

    for (int i=0; i<nr-2 ; i+=2) {
	som += (cc[0] - cc[1] ) / double(2*i + 2) ; 
	cc += 2 ;
    }
    
    if (nr%2 == 0) som += (*cc) / double(2*nr - 2) ;

    return som ; 
    
}
}
