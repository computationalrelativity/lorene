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


char comb_lin_cpt_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:28  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:16:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:37:11  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/11/22  19:29:50  eric
 * Changement de nom_C en comb_lin_cpt_C
 * Nettoyage des includes
 *
 * Revision 2.0  2000/03/16  16:25:18  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

// Headers Lorene
#include "matrice.h"

/*
 * Gestion des CL permettant de mettre a bande les operateurs associes a poisson 
 * compact. Version pour les sources et les matrices des operateurs.
 */


// Version Matrice --> Matrice
namespace Lorene {
Matrice _cl_cpt_pas_prevu (const Matrice &source, int) {
    cout << "Combinaison lineaire pas prevu..." << endl ;
    cout << "Source : " << source << endl ;
    abort() ;
    return source;
}


		//-------------------
	       //--  R_CHEBP   -----
	      //-------------------


Matrice _cl_cpt_r_chebp (const Matrice &source, int) {
    
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;

    Matrice barre(source) ;
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	    barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j))/(i+1) ;
	if (i==0) dirac = 0 ;
    }

    Matrice res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = barre(i, j)-barre(i+2, j) ;
    
    res.set_band(4, 1) ;
    res.set_lu() ;
    return res ;
}

		//-------------------
	       //--  R_CHEBI   -----
	      //-------------------


Matrice _cl_cpt_r_chebi (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
   
    Matrice barre(source) ;
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n ; j++)
	    barre.set(i, j) = (source(i, j)-source(i+1, j))/(i+1) ;
    
    Matrice res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = barre(i, j)-barre(i+2, j) ;
    
    if (l==1)
	res.set_band(3, 0) ;
    else
	res.set_band(3, 1) ;
    res.set_lu() ;
    return res ;
    
}
		

		//-------------------------
	       //- La routine a appeler ---
	      //---------------------------

Matrice combinaison_cpt (const Matrice &source, int l, int base_r) {
    
		// Routines de derivation
    static Matrice (*combinaison_cpt[MAX_BASE])
			(const Matrice &, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    combinaison_cpt[i] = _cl_cpt_pas_prevu ;
	}
		// Les routines existantes
	combinaison_cpt[R_CHEBP >> TRA_R] = _cl_cpt_r_chebp ;
	combinaison_cpt[R_CHEBI >> TRA_R] = _cl_cpt_r_chebi ;
    }
    
    Matrice res(combinaison_cpt[base_r](source, l)) ;
    return res ;
}

	    //--------------------------------------------------------------
	    //		    Version Tbl
	    //--------------------------------------------------------------



Tbl _cl_cpt_pas_prevu(const Tbl& tb) {
    cout << "combinaison_nul_pas_prevu " << endl ;
    cout << "tb : " << tb << endl ;
    abort() ;
    return tb ;
}

	
		//-------------------
	       //--  R_CHEBP   -----
	      //-------------------
Tbl _cl_cpt_r_chebp(const Tbl& tb) {
    
    assert (tb.get_etat() != ETATNONDEF) ;
    int n=tb.get_dim(0) ;
    
    Tbl barre(tb) ;
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	barre.set(i) = ((1+dirac)*tb(i)-tb(i+2))/(i+1) ;
	if (i==0) dirac = 0 ;
    }

    Tbl res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	res.set(i) = barre(i)-barre(i+2) ;
		
    return res ;
}


	
		//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
Tbl _cl_cpt_r_chebi(const Tbl& tb) {
   
    assert (tb.get_etat() != ETATNONDEF) ;
    int n=tb.get_dim(0) ;
    
    Tbl barre(tb) ;
    for (int i=0 ; i<n-2 ; i++)
	barre.set(i) = (tb(i)-tb(i+1))/(i+1) ;

    Tbl res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	res.set(i) = barre(i)-barre(i+2) ;
		
    return res ;
}  


		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl combinaison_cpt (const Tbl &source, int base_r) {
    
		// Routines de derivation
    static Tbl (*combinaison_cpt[MAX_BASE])(const Tbl&) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    combinaison_cpt[i] = _cl_cpt_pas_prevu ;
	}
		// Les routines existantes
	combinaison_cpt[R_CHEBP >> TRA_R] = _cl_cpt_r_chebp ;
	combinaison_cpt[R_CHEBI >> TRA_R] = _cl_cpt_r_chebi ;
    }
    
    Tbl res(combinaison_cpt[base_r](source)) ;
    return res ;
}

}
