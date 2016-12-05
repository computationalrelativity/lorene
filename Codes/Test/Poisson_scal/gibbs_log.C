/*
 * Test of scalar Poisson equation for a source that accounts for a Gibbs-like
 *  phenomenon
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:29  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/06 15:12:54  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:54  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/02/15  14:53:03  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */

// LORENE

#include "itbl.h"
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

Tbl donne_pl (int l) {
    assert (l>=0) ;
    
    Tbl plp1(l+1) ;
    Tbl pl(l+1) ;
    Tbl plm1(l+1) ;
    
    double a, b ;
    //Initialisation : 
    plm1.annule_hard() ;
    plm1.set(0) = 1 ;
    if (l==0)
	return plm1 ;
    assert (l != 0) ;
    
    pl.annule_hard() ;
    pl.set(1) = sqrt(3.) ;
    if (l==1)
	return pl ;
    
    plp1.set_etat_qcq() ;
    
    assert (l>=2) ;
    for (int conte = 1 ; conte<l ; conte++) {
	
	a = (conte+1)/sqrt((2*conte+1)*(2*conte+3)) ;
	b = conte/sqrt((2*conte+1)*(2*conte-1)) ;
	
	plp1.set(0) = -b*plm1(0) ;
	for (int i=1 ; i<l+1 ; i++)
	    plp1.set(i) = pl(i-1)-b*plm1(i) ;
	
	plm1 = pl ;
	pl = plp1/a ;
    }
    return pl ;
}

void main() {
    
    // nbres de points en r possibles ...
    Itbl r_pos (10) ;
    r_pos.set_etat_qcq() ;
    r_pos.set(0) = 7 ;
    r_pos.set(1) = 13 ;
    r_pos.set(2) = 17 ;
    r_pos.set(3) = 19 ;
    r_pos.set(4) = 21 ;
    r_pos.set(5) = 25 ;
    r_pos.set(6) = 33 ;
    r_pos.set(7) = 49 ;
    r_pos.set(8) = 65 ;
    r_pos.set(9) = 129 ;
    
    int puis_zec = 4 ;
    
    int l ;
    cin >> l ;
    assert (l>=0) ;
    if (l==0)
	assert (puis_zec == 2) ;
    
    
    
    // Construction de la grille ...
    int nz = 3 ;
    double R = (nz-1) ;
    
    //Construction du mapping :
    double* bornes = new double[nz+1] ;
    for (int i=0 ; i<nz ; i++)
	bornes[i] = i ;
    bornes[nz] = __infinity ;
    
    //On determine le nombre de points minimal en theta :
    int ntmin = (l%2==0) ? l/2+1 : (l+3)/2;
    assert (ntmin <= 129) ;
    int nbret = 0 ;
    int conte = 0 ;
    while (ntmin > r_pos(conte))
	conte++ ;
    nbret = r_pos(conte) ;
    
    int symetrie = (l%2 == 0) ? SYM : NONSYM ;
    int nbrep = (l%2 == 0) ? 1 : 4 ;
    // echantillonnage en phi :
    int* np = new int [nz] ;
    for (int lz=0 ; lz<nz ; lz++)
	np[lz] = nbrep ;
    int type_p = symetrie ;
    
    // echantillonnage en theta :
    int* nt = new int [nz] ;
    for (int lz=0 ; lz<nz ; lz++)
	nt[lz] = nbret ;
    int type_t = SYM ;
    
    // echantillonage en r :
    int* nr = new int [nz] ;
    int* type_r = new int[nz] ;
    type_r[0] = RARE ;
    for (int lz=1 ; lz<nz-1 ; lz++)
	type_r[lz] = FIN ;
    type_r[nz-1] = UNSURR ;
     
    // C'est parti pour la grande boucle (vas y Richard ...)
    int nrmin = 2*nbret ;
    assert (nrmin <=129) ;
    conte = 0 ;
    while (nrmin > r_pos(conte))
	conte++ ;

    int borne = (l%2 == 0) ? 10 : 9 ;
    
    for (int boucle = conte ; boucle<borne ; boucle ++) {
    
	for (int lz=0 ; lz<nz ; lz++)
	    nr[lz] = r_pos(boucle) ;
    

 
	Mg3d grille (nz, nr, type_r, nt, type_t, np, type_p) ;
    
   
	Map_af mapping(grille, bornes) ;
	
	Base_val** bases ;
	bases = grille.std_base_vect_cart() ;
	
	Base_val base(nz) ;
	
	if (l%2 ==0)
	    base = *bases[0] ;
	else
	    base = *bases[2] ;
	    
	for (int i=0 ; i<3 ; i++)
	    delete bases[i] ;
	delete [] bases ;
	
	//Les coord :
	Coord& r = mapping.r ;
	Coord& z = mapping.z ;
	Coord& cost = mapping.cost ;
    
	//Tableau contenant les coefficients du Pl considere :
	Tbl pl (donne_pl(l)) ;
    
	//Construction des Mtbl contenant Pl dans la zec et r^l*Pl dans ZI :
	Mtbl pl_zec (grille) ;
	pl_zec.annule_hard() ;
	for (int i=l ; i>=0 ; i-=2)
	    pl_zec += pl(i)*pow(cost, double(i)) ;
	pl_zec.annule(0, nz-2) ;
    
	Mtbl pl_zin (grille) ;
	pl_zin.annule_hard() ;
	for (int i=l ; i>=0 ; i-=2)
	    pl_zin += pl(i)*pow(z, double(i))*pow(r, double(l-i)) ;
	pl_zin.annule(nz-1, nz-1) ;
    
	// Construction de la source.
	Valeur val_so(&grille) ;
	val_so = -pl_zec/pow(r, double(l+3-puis_zec)) ;
	val_so.annule(0, nz-2) ;
	Cmp source(mapping) ;
	source =  val_so ;
	source.set_dzpuis(puis_zec) ;
	source.va.set_base(base) ;
    
	// On construit la solution analytique
	Valeur val_sol(&grille) ;
	val_sol = pl_zin/pow(2*l+1, 2)/pow(R, double(2*l+1)) ;

	Valeur val_sol_zec(&grille) ;
	val_sol_zec = pl_zec*(log(r)+1./(2*l+1)-log(R))/(2*l+1)/pow(r, double(l+1)) ;
	for (int k=0 ; k<np[nz-1] ; k++)
	    for (int j=0 ; j<nt[nz-1] ; j++)
		val_sol_zec.set(nz-1, k, j, nr[nz-1]-1) = 0. ;
	*(val_sol.c->t[nz-1]) = *(val_sol_zec.c->t[nz-1]) ; 

	// On construit le vecteur ...
	Cmp soluce(mapping) ;
	soluce = val_sol ;
	soluce.va.set_base(base) ;
	soluce.set_dzpuis(0) ;
    
	Cmp sol_poisson(source.poisson()) ;

        cout << nr[0] << " " ;
	
	Tbl erreur(diffrel(soluce, sol_poisson)) ;
	for (int i=0 ; i<nz ; i++)
	    cout << erreur(i) << " " ;
	cout << endl ;
    }
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    delete [] bornes ;
}

