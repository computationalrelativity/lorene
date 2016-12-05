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
 * Fournit la matrice de passage pour la transformation des coefficients du
 * developpement en fonctions associees de Legendre 
 * P_l^m(cos(theta)) impaires (i.e. telles que l-m est impair)
 * dans les coefficients du developpement 
 * en cos((2*j+1)*theta) [m pair] / sin( 2*j * theta) [m impair]
 * 
 * Cette routine n'effectue le calcul de la matrice que si celui-ci n'a pas 
 * deja ete fait, sinon elle renvoie le pointeur sur une valeur precedemment 
 * calculee.
 * 
 * Entree:
 * -------
 *  int np  :	Nombre de degres de liberte en phi 
 *  int nt  :	Nombre de degres de liberte en theta
 *
 * Sortie (valeur de retour) :
 * ---------------------------
 *  double* mat_legi_cossinci : pointeur sur le tableau contenant l'ensemble
 *			        (pour les np/2+1 valeurs de m) des 
 *				matrices de passage.
 *				La dimension du tableau est (np/2+1)*nt^2
 *				Le stokage est le suivant: 
 *
 *		mat_legi_cossinci[ nt*nt* m + nt*j + l] = B_{mjl}
 *			    
 *				ou B_{mjl} est defini par
 *
 * pour m pair :  
 *   P_{2l+1}^m( cos(theta) ) = som_{j=0}^{nt-2} B_{mjl} cos((2*j+1)*theta)
 *
 * pour m impair :  
 *   P_{2l}^m( cos(theta) ) = som_{j=1}^{nt-2} B_{mjl} sin(2*j*theta) 
 *
 *  ou P_n^m(x) represente la fonction de Legendre associee de degre n et 
 *     d'ordre m normalisee de facon a ce que
 *
 *	 int_0^pi [ P_n^m(cos(theta)) ]^2  sin(theta) dtheta = 1
 *
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:18:02  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:13  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:16:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/02/18 13:14:15  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.2  2002/10/16 14:36:55  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/22  15:34:21  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C
#include <cstdlib>
#include <cmath>

// Prototypage
#include "headcpp.h"
#include "proto.h"

namespace Lorene {
//******************************************************************************

double* mat_legi_cossinci(int np, int nt) {

#define NMAX	30		// Nombre maximun de couples(np,nt) differents 
static	double*	tab[NMAX] ;	// Tableau des pointeurs sur les tableaux 
static	int	nb_dejafait = 0 ;	// Nombre de tableaux deja initialises 
static	int	np_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait
static	int	nt_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait

int i, indice,  j,  j2,  m,  l ;

    {
    // Les matrices B_{mjl} pour ce couple (np,nt) ont-elles deja ete calculees ? 
    indice = -1 ;
    for ( i=0 ; i < nb_dejafait ; i++ ) {
	if ( (np_dejafait[i] == np) && (nt_dejafait[i] == nt) ) indice = i ;
	}


// Si le calcul n'a pas deja ete fait, il faut le faire : 
    if (indice == -1) {		   
	if ( nb_dejafait >= NMAX ) {
	    cout << "mat_legi_cossinci: nb_dejafait >= NMAX : "
		<< nb_dejafait << " <-> " << NMAX << endl ;
	    abort () ;	
	    exit(-1) ;
	    }
	indice = nb_dejafait ; 
	nb_dejafait++ ; 
	np_dejafait[indice] = np ;
	nt_dejafait[indice] = nt ;

	tab[indice] = new double[(np/2+1)*nt*nt] ; //(double *) malloc( sizeof(double) * (np/2+1)*nt*nt ) ;

//-----------------------
// Preparation du calcul 
//-----------------------

// Sur-echantillonnage :
	int nt2 = 2*nt - 1 ; 

	int deg[3] ;
	deg[0] = 1 ;
	deg[1] = nt2 ;
	deg[2] = 1 ;

// Tableaux de travail
	double* yy = new double[nt2] ; //(double*)( malloc( nt2*sizeof(double) ) ) ;


//-------------------
// Boucle sur m
//-------------------

	for (m=0; m < np/2+1 ; m++) {

// Recherche des fonctions de Legendre associees d'ordre m :

	    double* leg = legendre_norm(m, nt) ;	
	
	    if (m%2==0) {
// Cas m pair
//-----------
		for (l=m/2; l<nt-1; l++) {	// boucle sur les P_{2l+1}^m
	    
		    int ll = 2*l+1 ;	// degre des fonctions de Legendre
		    
		    for (j2=0; j2<nt2; j2++) {
			    yy[j2] = leg[nt2* (ll-m) + j2] ;
		    }

//....... transformation en cos((2*j+1)*theta) : 

		    cftcosi(deg, deg, yy, deg, yy) ; 

//....... le resultat fournit les elements de matrice : 
		    for (j=0; j<nt; j++) {
			tab[indice][ nt*nt* m + nt*j + l] = yy[j] ; 
		    }
	    
		}  // fin de la boucle sur l (indice de P_{2l}^m)
	    
	    
	    }   // fin du cas m pair
	    else {
		
// Cas m impair
//-------------

		for (l=(m+1)/2; l<nt-1; l++) {	// boucle sur les P_{2l}^m
	    
		    int ll = 2*l ;	// degre des fonctions de Legendre
		    
			
		    for (j2=0; j2<nt2; j2++) {
			    yy[j2] = leg[nt2* (ll-m) + j2] ;
		    }

//....... transformation en sin( 2j*theta) : 

		    cftsinp(deg, deg, yy, deg, yy) ; 

//....... le resultat fournit les elements de matrice : 
		    for (j=0; j<nt-1; j++) {
			tab[indice][ nt*nt* m + nt*j + l] = yy[j] ; 
		    }
	    
		}  // fin de la boucle sur l (indice de P_{2l+1}^m)
	    

	    } // fin du cas m impair
	    
	    delete [] leg ;
	   	    
	}  // fin de la boucle sur m

// Liberation espace memoire
// -------------------------

	delete [] yy ;

    } // fin du cas ou le calcul etait necessaire

    }	// Fin de zone critique
    
    return tab[indice] ;

}


}
