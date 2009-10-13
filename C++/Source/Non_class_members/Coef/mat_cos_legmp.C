/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *                 2009 Jerome Novak
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


char mat_cos_legmp_C[] = "$Header$" ;

/*
 * Fournit la matrice de passage pour la transformation des coefficients du
 * developpement en cos(j*theta) 
 * dans les coefficients du developpement en fonctions associees de Legendre 
 * P_l^m(cos(theta)) telles que m est pair.
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
 *  double* mat_cos_legmp :     pointeur sur le tableau contenant l'ensemble
 *			        (pour les np/2+1 valeurs de m: m=0,2,...,np) des 
 *				matrices de passage.
 *				La dimension du tableau est (np/2+1)*nt^2
 *				Le stokage est le suivant: 
 *
 *		mat_cos_legmp[ nt*nt* m/2 + nt*l + j] = A_{mlj}
 *			    
 *				ou A_{mlj} est defini par
 *
 *   cos(j*theta) = som_{l=m}^{nt-1} A_{mlj} P_l^m( cos(theta) )
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
 * Revision 1.1  2009/10/13 13:49:36  j_novak
 * New base T_LEG_MP.
 *
 *
 *
 * $Header$
 *
 */

// headers du C
#include <stdlib.h>
#include <math.h>

// Prototypage
#include "headcpp.h"
#include "proto.h"

//******************************************************************************

double* mat_cos_legmp(int np, int nt) {

#define NMAX	30		// Nombre maximun de couples(np,nt) differents 
    static	double*	tab[NMAX] ;	// Tableau des pointeurs sur les tableaux 
    static	int	nb_dejafait = 0 ;      // Nombre de tableaux deja initialises 
    static	int	np_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
                                               // calcul a deja ete fait
    static	int	nt_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
                                               // calcul a deja ete fait

    int i, indice,  j,  j2,  m,  l ;
    
	
    // Les matrices A_{mlj} pour ce couple (np,nt) ont-elles deja ete calculees ? 
    indice = -1 ;
    for ( i=0 ; i < nb_dejafait ; i++ ) {
	if ( (np_dejafait[i] == np) && (nt_dejafait[i] == nt) ) indice = i ;
    }


// Si le calcul n'a pas deja ete fait, il faut le faire : 
    if (indice == -1) {		   
	if ( nb_dejafait >= NMAX ) {
	    cout << "mat_cos_legmp: nb_dejafait >= NMAX : "
		 << nb_dejafait << " <-> " << NMAX << endl ;
	    abort () ;	
	    exit(-1) ;
	}
	indice = nb_dejafait ; 
	nb_dejafait++ ; 
	np_dejafait[indice] = np ;
	nt_dejafait[indice] = nt ;
	
	tab[indice] = new double[(np/2+1)*nt*nt] ;
	
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
	double* yy = new double[nt2] ; 
	double* cost = new double[nt*nt2] ; 
   
// Calcul des cos(j*theta)  aux points de collocation
//  de l'echantillonnage double : 

	double dt = M_PI / double(nt2-1) ;
	for (j=0; j<nt; j++) {
	    for (j2=0; j2<nt2; j2++) {
		double theta = j2*dt ;
		cost[nt2*j + j2] = cos( j * theta ) ;
	    }
	}	


//-------------------
// Boucle sur m
//-------------------

	int m_max = np ; 
	if (np == 1) m_max = 0 ; 

	for (m=0; m <= m_max ; m+=2) {

// Recherche des fonctions de Legendre associees d'ordre m :

	    double* leg = legendre_norm(m, nt) ;	
	
	    for (l=m; l<nt; l++) {	// boucle sur les P_l^m
	    
		  int parite = 1 - 2*((l-m)%2) ;  // parite du P_l^m par rapport au plan theta=pi/2
		for (j=0; j<nt; j++) {	// boucle sur les cos(j theta)
			
//... produit scalaire de cos(j theta) par P_l^m(cos(theta)) 

		    for (j2=0; j2<nt; j2++) {
			yy[nt2m1-j2] = cost[nt2*j + j2] *
					   leg[nt2* (l-m) + 2*j2] ;
		    }

		    for (j2 = nt; j2<nt2; j2++) {
			yy[nt2m1-j2] = cost[nt2*j + j2] *
			    parite * leg[nt2* (l-m) + 2*nt2 -2 -2*j2] ;
		    }

//....... on passe en Tchebyshev vis-a-vis de x=cos(theta) pour	calculer
//	    l'integrale (routine int1d_cheb) : 		
		    cfrcheb(deg, deg, yy, deg, yy) ;
		    tab[indice][ nt*nt* m/2 + nt*l + j] = 
					    int1d_cheb(nt2, yy) ;

		}	// fin de la boucle sur j  (indice de cos(j theta) )
	    
	    }  // fin de la boucle sur l (indice de P_l^m)
	    	    
	delete [] leg ; 
	   	    
	}  // fin de la boucle sur m

// Liberation espace memoire
// -------------------------

	delete [] yy ;
	delete [] cost ;

    } // fin du cas ou le calcul etait necessaire
    
    return tab[indice] ;

}


