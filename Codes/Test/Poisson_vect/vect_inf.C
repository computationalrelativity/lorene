/*
 * Test for the vectorial Poisson equation for a discontinuous non-compact source
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

char vect_inf_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/01/09 11:07:56  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/02/15  14:57:18  phil
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
#include "utilitaires.h"
#include "coord.h"
#include "cmp.h"
#include "graphique.h"
#include "proto.h"
#include "tenseur.h"
#include "base_vect.h"
#include "param.h"

//standard
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void main() { 

    int symetrie = NONSYM ;
    double lambda = 1./3. ;
    int n=4 ;
    
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
	np[l] = 8 ;
    int type_p = symetrie ;
    
    // echantillonnage en theta :
    int* nt = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	nt[l] = 9 ;
    int type_t = SYM ;
    int* nr = new int[nz] ;
    int* type_r = new int[nz] ;
    int nbrer ;
    
    cin >> nbrer ;

	// echantillonage en r :
	
	for (int l=0 ; l<nz ; l++)
	    nr[l] = nbrer ;
	
	type_r[0] = RARE ;
	for (int l=1 ; l<nz-1 ; l++)
	    type_r[l] = FIN ;
	type_r[nz-1] = UNSURR ;
    
	Mg3d grille (nz, nr, type_r, nt, type_t, np, type_p) ;
    
   
	Map_af mapping(grille, bornes) ;

    
	    // Construction de la source .
	    
	Coord& x = mapping.x ;
	Coord& y = mapping.y ;
	Coord& z = mapping.z ;
	Coord& r = mapping.r ;
	
	Mtbl a_cont (&grille) ;
	a_cont = -(4.+n)/2./pow(R, 6.+n) ;
	Mtbl b_cont (&grille) ;
	b_cont = (6.+n)/2./pow(R, 4.+n) ;
	
	// Composante x
	Valeur val_x(&grille) ;
	val_x = x*((54+18*lambda)*a_cont*pow(r, 4.)+
		    (28+12*lambda)*b_cont*pow(r, 2.))
		+ lambda*pow(x, 3.)*(24*a_cont*pow(r, 2.)+8*b_cont) ;
	// ZEC ...
	Valeur val_x_zec(&grille) ;
	val_x_zec = n*(n-3-3*lambda)*x/pow(r, n)+
		    n*(n+2)*lambda*pow(x, 3.)/pow(r, n+2.)  ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		val_x_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(val_x.c->t[nz-1]) = *(val_x_zec.c->t[nz-1]) ; 
	
	// Composante y
	Valeur val_y(&grille) ;
	val_y = lambda*y*(6*a_cont*pow(r, 4.)+4*b_cont*pow(r, 2.)) 
		+lambda*pow(x, 2.)*y*(24*a_cont*pow(r, 2.)+8*b_cont);
	// ZEC ...
	Valeur val_y_zec(&grille) ;
	val_y_zec = -lambda*n*y/pow(r, n)+lambda*n*(n+2)*pow(x, 2.)*y/pow(r, n+2.)  ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		val_y_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(val_y.c->t[nz-1]) = *(val_y_zec.c->t[nz-1]) ; 
	
	// Composante z
	Valeur val_z(&grille) ;
	val_z = lambda*z*(6*a_cont*pow(r, 4.)+4*b_cont*pow(r, 2.)) 
		+lambda*pow(x, 2.)*z*(24*a_cont*pow(r, 2.)+8*b_cont);
	// ZEC ...
	Valeur val_z_zec(&grille) ;
	val_z_zec = -lambda*n*z/pow(r, n)+lambda*n*(n+2)*pow(x, 2.)*z/pow(r, n+2.)  ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		val_z_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(val_z.c->t[nz-1]) = *(val_z_zec.c->t[nz-1]) ; 
	
	// On construit le vecteur ...
	Tenseur source (mapping, 1, CON, mapping.get_bvect_cart()) ;
	source.set_etat_qcq() ;
	source.set(0) = val_x ;
	source.set(1) = val_y ;
	source.set(2) = val_z ;
	for (int i=0 ; i<3 ; i++)
	    source.set(i).set_dzpuis(2) ;
	source.set_std_base() ;
	
	
	// On construit la solution analytique
	// Composante x
	Valeur sol_x(&grille) ;
	sol_x = x*(a_cont*pow(r, 6.)+b_cont*pow(r, 4.)) ;
    
	// ZEC ...
	Valeur sol_x_zec(&grille) ;
	sol_x_zec = x/pow(r, n) ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		sol_x_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(sol_x.c->t[nz-1]) = *(sol_x_zec.c->t[nz-1]) ; 

	// On construit le vecteur ...
	Tenseur soluce(mapping, 1, CON, mapping.get_bvect_cart()) ;
	soluce.set_etat_qcq() ;
	soluce.set(0) = sol_x ;
	soluce.set(1) = 0 ;
	soluce.set(2) = 0 ;
	for (int i=0 ; i<3 ; i++)
	    soluce.set(i).set_dzpuis(0) ;
	soluce.set_std_base() ;
	
	Tenseur vect_auxi (mapping, 1, CON, mapping.get_bvect_cart()) ;
	vect_auxi.set_etat_qcq() ;
	Tenseur scal_auxi (mapping) ;
	scal_auxi.set_etat_qcq() ;
	
	Tenseur shibata (source.poisson_vect(lambda, vect_auxi, scal_auxi)) ;
	
	cout << "Shibata" << endl ;
	for (int i=0 ; i<3 ; i++)
	    cout << "Composante " << i << " : " <<
		diffrelmax(shibata(i), soluce(i)) << endl ;
	cout <<"-------------------------------------------" << endl ;
	
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    delete [] bornes ;
}    
