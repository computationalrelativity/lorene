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


char mat_cossinci_legi_C[] = "$Header$" ;

/*
 * Fournit la matrice de passage pour la transformation des coefficients du
 * developpement en cos((2*j+1)*theta) [m pair] / sin( 2*j * theta) [m impair]
 * dans les coefficients du developpement en fonctions associees de Legendre 
 * P_l^m(cos(theta)) impaires (i.e. telles que l-m est impair). 
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
 *  double* mat_cossinci_legi : pointeur sur le tableau contenant l'ensemble
 *			        (pour les np/2+1 valeurs de m) des 
 *				matrices de passage.
 *				La dimension du tableau est (np/2+1)*nt^2
 *				Le stokage est le suivant: 
 *
 *		mat_cossinci_legi[ nt*nt* m + nt*l + j] = A_{mlj}
 *			    
 *				ou A_{mlj} est defini par
 *
 * pour m pair :  
 *   cos((2*j+1)*theta) = som_{l=m/2}^{nt-2} A_{mlj} P_{2l+1}^m( cos(theta) )
 *						pour 0 <= j <= nt-2
 *
 * pour m impair :  
 *   sin(2*j*theta) = som_{l=(m+1)/2}^{nt-2} A_{mlj} P_{2l}^m( cos(theta) )
 *						pour 1 <= j <= nt-2
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
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  1999/02/22  15:35:09  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// headers du C++
#include <iostream.h>

// headers du C
#include <stdlib.h>
#include <math.h>

// Prototypage
#include "proto.h"

// Variable de loch
int loch_mat_cossinci_legi = 0 ;

//******************************************************************************

double* mat_cossinci_legi(int np, int nt) {

#define NMAX	30		// Nombre maximun de couples(np,nt) differents 
static	double*	tab[NMAX] ;	// Tableau des pointeurs sur les tableaux 
static	int	nb_dejafait = 0 ;	// Nombre de tableaux deja initialises 
static	int	np_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait
static	int	nt_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait

int i, indice,  j,  j2,  m,  l ;
	
    #pragma critical (loch_mat_cossinci_legi)
    { 

    // Les matrices A_{mlj} pour ce couple (np,nt) ont-elles deja ete calculees ? 
    indice = -1 ;
    for ( i=0 ; i < nb_dejafait ; i++ ) {
	if ( (np_dejafait[i] == np) && (nt_dejafait[i] == nt) ) indice = i ;
	}


    // Si le calcul n'a pas deja ete fait, il faut le faire : 
    if (indice == -1) {		   
	if ( nb_dejafait >= NMAX ) {
	    cout << "mat_cossinci_legi: nb_dejafait >= NMAX : "
		<< nb_dejafait << " <-> " << NMAX << endl ;
	    abort () ;	
	    exit(-1) ;
	    }
	indice = nb_dejafait ; 
	nb_dejafait++ ; 
	np_dejafait[indice] = np ;
	nt_dejafait[indice] = nt ;

	tab[indice] = (double *) malloc( sizeof(double) * (np/2+1)*nt*nt ) ;

//-----------------------
// Preparation du calcul 
//-----------------------

// Sur-echantillonnage pour calculer les produits scalaires sans aliasing:
	int nt2 = 2*nt - 1 ; 
	int nt2m1 = nt2 - 1 ; 

	int deg[3] ;
	deg[0] = 1 ;
	deg[1] = 1 ;
	deg[2] = nt2 ;

// Tableaux de travail
	double* yy = (double*)( malloc( nt2*sizeof(double) ) ) ;
	double* cost = (double*)( malloc( nt*nt2*sizeof(double) ) ) ;
	double* sint = (double*)( malloc( nt*nt2*sizeof(double) ) ) ;
   
// Calcul des cos(2*j*theta) / sin( (2*j+1)*theta ) aux points de collocation
//  de l'echantillonnage double : 

	double dt = M_PI / double(2*(nt2-1)) ;
	for (j=0; j<nt-1; j++) {
	    for (j2=0; j2<nt2; j2++) {
		double theta = j2*dt ;
		cost[nt2*j + j2] = cos( (2*j+1) * theta ) ;
		sint[nt2*j + j2] = sin( 2*j * theta ) ;
	    }
	}	


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
		    
		    for (j=0; j<nt-1; j++) { // boucle sur les cos((2j+1) theta)
			
//... produit scalaire de cos((2j+1) theta) par P_{2l+1}^m(cos(theta)) 

			for (j2=0; j2<nt2; j2++) {
			    yy[nt2m1-j2] = cost[nt2*j + j2] *
					   leg[nt2* (ll-m) + j2] ;
			}

//....... on passe en Tchebyshev vis-a-vis de x=cos(theta) pour	calculer
//	    l'integrale (routine int1d_chebp) : 		
			cfrchebp(deg, deg, yy, deg, yy) ;
			tab[indice][ nt*nt* m + nt*l + j] = 
					    2.*int1d_chebp(nt2, yy) ;

		    }	// fin de la boucle sur j  (indice de cos((2j+1) theta) )
	    
		}  // fin de la boucle sur l (indice de P_{2l+1}^m)
	    
	    
	    }   // fin du cas m pair
	    else {
		
// Cas m impair
//-------------

		for (l=(m+1)/2; l<nt-1; l++) {	// boucle sur les P_{2l}^m
	    
		    int ll = 2*l ;	// degre des fonctions de Legendre
		    
		    for (j=0; j<nt-1; j++) {  // boucle sur les sin(2j theta)
			
//... produit scalaire de sin((2j+1) theta) par P_{2l+1}^m(cos(theta)) 

			for (j2=0; j2<nt2; j2++) {
			    yy[nt2m1-j2] = sint[nt2*j + j2] *
					   leg[nt2* (ll-m) + j2] ;
			}

//....... on passe en Tchebyshev vis-a-vis de x=cos(theta) pour	calculer
//	    l'integrale (routine int1d_chebp) : 		
			cfrchebp(deg, deg, yy, deg, yy) ;
			tab[indice][ nt*nt* m + nt*l + j] = 
					    2.*int1d_chebp(nt2, yy) ;

		    }	// fin de la boucle sur j  (indice de  sin(2j theta) )
	    
		}  // fin de la boucle sur l (indice de P_{2l}^m)
	    

	    } // fin du cas m impair
	    
	free(leg) ;
	   	    
	}  // fin de la boucle sur m

// Liberation espace memoire
// -------------------------

	free(yy) ;
	free(cost) ; 
	free(sint) ; 

    } // fin du cas ou le calcul etait necessaire
    
    } // Fin de zone critique
    
    return tab[indice] ;

}


