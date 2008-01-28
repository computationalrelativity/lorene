/*
 * Test of scalar Poisson equation for a spherical non-compact source
 * 
 */
 
/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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

char test_poisson_spherique_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2008/01/28 16:08:47  j_novak
 * Changed return type of main from void to int
 *
 * Revision 1.2  2003/01/09 11:07:55  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/02/15  14:54:42  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */

// LORENE

#include "type_parite.h"
#include "nbr_spx.h"
#include "grilles.h"
#include "map.h"
#include "valeur.h"
#include "proto.h"
#include "coord.h"
#include "cmp.h"
#include "graphique.h"

//standard
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main() { 
    int symetrie = SYM ;
    int puis_zec = 2 ;
    
    // Construction de la grille ...
    int nz = 3 ;
    double R = (nz-1) ;
    
    //Construction du mapping :
    double* bornes = new double[nz+1] ;
    for (int i=0 ; i<nz ; i++)
	bornes[i] = i ;
    bornes[nz] = __infinity ;
    
    // echantillonnage en phi :
    int* np = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	np[l] = 1 ;
    int type_p = symetrie ;
    
    // echantillonnage en theta :
    int* nt = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	nt[l] = 1 ;
    int type_t = SYM ;
    
    int nbrer ;
    cin >> nbrer;
    
    // echantillonage en r :
    int* nr = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	nr[l] = nbrer ;
	

    int* type_r = new int[nz] ;
    type_r[0] = RARE ;
    for (int l=1 ; l<nz-1 ; l++)
	type_r[l] = FIN ;
    type_r[nz-1] = UNSURR ;
    
    Mg3d grille (nz, nr, type_r, nt, type_t, np, type_p) ;
    
   
    Map_af mapping(grille, bornes) ;

    
	    // Construction de la source .
    Coord& r = mapping.r ;
   
    // Composante x
    Valeur val_x(&grille) ;
    val_x = R-pow(r, 2.)/R ;
    
    // ZEC ...
    Valeur val_x_zec(&grille) ;
    val_x_zec = pow(R, 5.)/pow(r, 4-puis_zec) ;
    *(val_x.c->t[nz-1]) = *(val_x_zec.c->t[nz-1]) ; 
    
    Cmp source(mapping) ;
    source =  val_x ;
    source.set_dzpuis(puis_zec) ;
    source.std_base_scal() ;
    
    
    // On construit la solution analytique
    
    // Composante x
    Valeur sol_x(&grille) ;
    sol_x = R/6.*pow(r, 2.)-pow(r, 4.)/20./R-3./4.*pow(R, 3.) ;
    
    // ZEC ...
    Valeur sol_x_zec(&grille) ;
    sol_x_zec = pow(R, 5.)/2./pow(r, 2.)-17./15.*pow(R, 4.)/r ;
    for (int k=0 ; k<np[nz-1] ; k++)
	for (int j=0 ; j<nt[nz-1] ; j++)
	    sol_x_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
    *(sol_x.c->t[nz-1]) = *(sol_x_zec.c->t[nz-1]) ; 
    
    
    // On construit le vecteur ...
    Cmp soluce(mapping) ;
    soluce = sol_x ;
    soluce.set_dzpuis(0) ;
    soluce.std_base_scal() ;
    
    Cmp sol_poisson(source.poisson()) ;
   
    cout << nr[0] << endl ;
    cout <<"Erreur : " << endl ;
    cout << diffrelmax(sol_poisson, soluce) << endl ;
   
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    delete [] bornes ;
}
