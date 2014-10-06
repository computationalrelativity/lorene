/*
 * Test of the vectorial Poisson euqation with a "vectorial" Gibbs-like 
 * phenomenon
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

char vect_gibbs_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:12:55  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:56  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/02/15  14:56:27  phil
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
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>

void main() { 

    int symetrie = NONSYM ;
    double lambda = 1./3. ;
    
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
	
	// Composante z
	Valeur val_z(&grille) ;
	val_z = z/pow(R, 7.) ;
	// ZEC ...
	Valeur val_z_zec(&grille) ;
	val_z_zec = z/pow(r, 5.) ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		val_z_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(val_z.c->t[nz-1]) = *(val_z_zec.c->t[nz-1]) ; 
	val_z.c->dzpuis = 2 ;
		
	// On construit le vecteur ...
	Tenseur source (mapping, 1, CON, mapping.get_bvect_cart()) ;
	source.set_etat_qcq() ;
	source.set(0) = 0 ;
	source.set(1) = 0 ;
	source.set(2) = val_z ;
	for (int i=0 ; i<3 ; i++)
	    source.set(i).set_dzpuis(2) ;
	source.set_std_base() ;
	
	
	// Composante x
	Valeur sol_x(&grille) ;
	sol_x = -lambda/2./(lambda+1)*(2.*x*pow(z, 2.)/35./pow(R, 7.)
				    + pow(r, 2.)*x/35./pow(R, 7.)
				    -7.*x/75./pow(R, 5.)) ;
	 ;
	// ZEC ...
	Valeur sol_x_zec(&grille) ;
	sol_x_zec = -lambda/2./(lambda+1.)*(pow(z, 2.)*x/pow(r, 7.)*(-9./14.+log(R)-log(r))
		    +x/pow(r, 5.)*(59./350.+(log(r)-log(R))/5.)
		    -7.*x/30./pow(R, 2.)/pow(r, 3.)
		    +7./10.*pow(z, 2.)*x/pow(R, 2.)/pow(r, 5.));
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		sol_x_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(sol_x.c->t[nz-1]) = *(sol_x_zec.c->t[nz-1]) ; 
	
	// Composante y
	Valeur sol_y(&grille) ;
	sol_y = -lambda/2./(lambda+1)*(2.*y*pow(z, 2.)/35./pow(R, 7.)
				    + pow(r, 2.)*y/35./pow(R, 7.)
				    -7.*y/75./pow(R, 5.)) ;
	// ZEC ...
	Valeur sol_y_zec(&grille) ;
	sol_y_zec = -lambda/2./(lambda+1.)*(pow(z, 2.)*y/pow(r, 7.)*(-9./14.+log(R)-log(r))
		    +y/pow(r, 5.)*(59./350.+(log(r)-log(R))/5.)
		    -7.*y/30./pow(R, 2.)/pow(r, 3.)
		    +7./10.*pow(z, 2.)*y/pow(R, 2.)/pow(r, 5.)) ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		sol_y_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(sol_y.c->t[nz-1]) = *(sol_y_zec.c->t[nz-1]) ; 
	
	// Composante z
	Valeur sol_z(&grille) ;
	sol_z = (lambda+2.)/2./(lambda+1.)*z*(pow(r, 2.)/10./pow(R, 7.)
					-7./30./pow(R, 5.))
	    -lambda/2./(lambda+1.)*(2.*pow(z, 3.)/35./pow(R, 7.)
				    -z*pow(r, 2.)/70./pow(R, 7.)
				    -7.*z/150./pow(R, 5.)) ;
	// ZEC ...
	Valeur sol_z_zec(&grille) ;
	sol_z_zec =  (lambda+2.)/2./(lambda+1.)*z* (1./10./pow(r, 5.)
						-7./30./pow(R, 2.)/pow(r, 3.))
		-lambda/2./(lambda+1.)*(pow(z, 3.)/pow(r, 7.)*(-9./14.
						+log(R)-log(r))
					+7./10.*pow(z, 3.)/pow(R, 2.)/pow(r, 5.)
				+z/pow(r, 5.)*(71./175.+3./5.*(log(r)-log(R)))
				-7./15.*z/pow(R, 2.)/pow(r, 3.)) ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		sol_z_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(sol_z.c->t[nz-1]) = *(sol_z_zec.c->t[nz-1]) ; 
	
	// On construit le vecteur ...
	Tenseur soluce (mapping, 1, CON, mapping.get_bvect_cart()) ;
	soluce.set_etat_qcq() ;
	soluce.set(0) = sol_x ;
	soluce.set(1) = sol_y ;
	soluce.set(2) = sol_z ;
	for (int i=0 ; i<3 ; i++)
	    soluce.set(i).set_dzpuis(0) ;
	soluce.set_std_base() ;
	
	
	Tenseur vect_auxi (mapping, 1, CON, mapping.get_bvect_cart()) ;
	vect_auxi.set_etat_qcq() ;
	Tenseur scal_auxi (mapping) ;
	scal_auxi.set_etat_qcq() ;
	Tenseur chi_auxi (mapping) ;
	chi_auxi.set_etat_qcq() ;
	
	
	Tenseur shibata (source.poisson_vect(lambda, vect_auxi, scal_auxi)) ;
	Tenseur oohara (source.poisson_vect_oohara(lambda, chi_auxi)) ;
	
	cout << "Points en r : " << nbrer << endl ;
	cout << "Shibata vs analytique" << endl ;
	for (int i=0 ; i<3 ; i++)
	    cout << diffrelmax(shibata(i), soluce(i)) << endl ;
    	
	cout << "Oohara vs analytique" << endl ;
	for (int i=0 ; i<3 ; i++)
	    cout << diffrelmax(oohara(i), soluce(i)) << endl ;
    	
	
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    delete [] bornes ;
}    
