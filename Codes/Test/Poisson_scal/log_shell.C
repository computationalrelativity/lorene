/*
 * Test of scalar Poisson equation for a solution with ln functions in a shell
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

char log_shell_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:12:54  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:55  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/02/15  14:53:51  phil
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
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>

void main() { 
    int symetrie = NONSYM ;
    
    // Construction de la grille ...
    int nz = 3 ;
    double R_2 = nz-1 ;
    double R_1 = nz-2 ;
    
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
    Coord& z = mapping.z ;
    
    // Composante x
    Valeur val_x(&grille) ;
    val_x = 1./pow(r, 3.) + 3*pow(z, 2.)/pow(r, 2.)-1. ;
    
    Cmp source(mapping) ;
    source =  val_x ;
    source.annule(0) ;
    source.annule(nz-1) ;
    source.std_base_scal() ;
   
    // On construit la solution analytique
    
    Mtbl a_cont(&grille) ;
    a_cont = 1./R_2-1./R_1 ;
    Mtbl b_cont(&grille) ;
    b_cont = (log(R_1)-log(R_2))/5. ;
    Mtbl c_cont(&grille) ;
    c_cont = 1./R_2 ;
    Mtbl d_cont(&grille) ;
    d_cont = log(R_1)-1. ;
    Mtbl e_cont(&grille) ;
    e_cont = -log(R_2)/5.-1./25. ;
    Mtbl f_cont(&grille) ;
    f_cont = pow(R_1, 5.)/25.;
    Mtbl g_cont(&grille) ;
    g_cont = log(R_1)-log(R_2) ;
    Mtbl h_cont(&grille) ;
    h_cont = (pow(R_1, 5.)-pow(R_2, 5.))/25. ;
    
    
    Valeur sol(&grille) ;
    sol = a_cont+ b_cont*(3*pow(z, 2.)-pow(r, 2.)) ;
 
    // Shell ...
    Valeur sol_shell(&grille) ;
    sol_shell = (3*pow(z, 2.)/pow(r, 2.)-1.)*
		(pow(r, 2.)*log(r)/5.+e_cont*pow(r, 2.)+f_cont/pow(r, 3.))
		-log(r)/r + c_cont +d_cont/r ;
    *(sol.c->t[nz-2]) = *(sol_shell.c->t[nz-2]) ; 

    // ZEC ...
    Valeur sol_zec(&grille) ;
    sol_zec = g_cont/r + h_cont*3*pow(z, 2.)/pow(r, 5.)-h_cont/pow(r, 3.) ;
    for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		sol_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
    *(sol.c->t[nz-1]) = *(sol_zec.c->t[nz-1]) ; 
  
    // On construit le vecteur ...
    Cmp soluce(mapping) ;
    soluce = sol ;
    soluce.set_dzpuis(0) ;
    soluce.std_base_scal() ;
    
    Cmp sol_poisson( source.poisson()) ;
    cout << nbrer << endl ;
    cout <<"Erreur : " << endl ;
    cout << diffrelmax(sol_poisson, soluce) << endl ;
   
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    delete [] bornes ;
}    

