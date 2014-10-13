/*
 *   Copyright (c) 2001 Philippe Grandclement
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


char separation_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:52:41  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:12:59  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/10/03 15:58:44  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2001/04/02  12:16:20  phil
 * *** empty log message ***
 *
 * Revision 2.5  2001/03/30  13:48:17  phil
 * on appelle raccord externe
 *
 * Revision 2.4  2001/03/22  10:40:30  phil
 * modification prototypage
 *
 * Revision 2.3  2001/03/02  10:19:05  phil
 * modification parametrage pour affichage
 *
 * Revision 2.2  2001/02/28  13:39:34  phil
 * modif cas etat_zero
 *
 * Revision 2.1  2001/02/28  13:23:00  phil
 * modif etat initial
 *
 * Revision 2.0  2001/02/28  11:24:34  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>

// Lorene
#include "cmp.h"
#include "proto.h"

namespace Lorene {
void separation (const Cmp& c1, const Cmp& c2, Cmp& res1, Cmp& res2, int decrois, 
    int puiss, int lmax, double precision, const double relax, const int itemax, const int flag) {
    
    assert (c1.get_etat() != ETATNONDEF) ;
    assert (c2.get_etat() != ETATNONDEF) ;
    
    if ((c1.get_etat() == ETATZERO) && (c2.get_etat() == ETATZERO)) {
	res1.set_etat_zero() ;
	res2.set_etat_zero() ;
	return ;
    }
    else {
    
    res1 = c1 ;
    if (res1.get_etat() == ETATZERO) {
	    res1.annule_hard() ;
	    res1.std_base_scal() ;
	    }
	    
    res1.raccord_externe (decrois, puiss, lmax) ;
    for (int i=0 ; i<decrois ; i++)
	res1.dec_dzpuis() ;
	
    res2 = c2 ;
    if (res2.get_etat() == ETATZERO) {
	    res2.annule_hard() ;
	    res2.std_base_scal() ;
	    }
    res2.raccord_externe (decrois, puiss, lmax) ;
    for (int i=0 ; i<decrois ; i++)
	res2.dec_dzpuis() ;
	
    int indic = 1 ;
    int conte = 0 ;
    // On commence la boucle pour separer :
    while (indic == 1) {
	Cmp old_un (res1) ;
	Cmp old_deux (res2) ;

	// On fait les modifications :
	Mtbl xa_mtbl_un (c1.get_mp()->xa) ;
	Mtbl ya_mtbl_un (c1.get_mp()->ya) ;
	Mtbl za_mtbl_un (c1.get_mp()->za) ;
	
	Mtbl xa_mtbl_deux (c2.get_mp()->xa) ;
	Mtbl ya_mtbl_deux (c2.get_mp()->ya) ;
	Mtbl za_mtbl_deux (c2.get_mp()->za) ;
    
	double xabs, yabs, zabs, air, theta, phi ;
	int np, nt, nr ;
	
	// On modifie le Cmp 1
	int nz_un = c1.get_mp()->get_mg()->get_nzone() ;
	for (int l=1 ; l<nz_un-1 ; l++) {
	    
	    np = c1.get_mp()->get_mg()->get_np(l) ;
	    nt = c1.get_mp()->get_mg()->get_nt(l) ;
	    nr = c1.get_mp()->get_mg()->get_nr(l) ;
	    
	    for (int k=0 ; k<np ; k++)
		for (int j=0 ; j<nt ; j++)
		    for (int i=0 ; i<nr ; i++) {
			xabs = xa_mtbl_un (l, k, j, i) ;
			yabs = ya_mtbl_un (l, k, j, i) ;
			zabs = za_mtbl_un (l, k, j, i) ;
			
			c2.get_mp()->convert_absolute(xabs, yabs, zabs, air, theta, phi) ;
			res1.set(l, k, j, i) = 
		(1-relax)*res1.set(l, k, j, i) +
	relax*(c1(l, k, j, i) - old_deux.val_point(air, theta, phi)) ;	
	    }
	}

	// On modifie le trou 2 
	int nz_deux = c2.get_mp()->get_mg()->get_nzone() ;
	for (int l=1 ; l<nz_deux-1 ; l++) {
	    
	    np = c2.get_mp()->get_mg()->get_np(l) ;
	    nt = c2.get_mp()->get_mg()->get_nt(l) ;
	    nr = c2.get_mp()->get_mg()->get_nr(l) ;
	    
	    for (int k=0 ; k<np ; k++)
		for (int j=0 ; j<nt ; j++)
		    for (int i=0 ; i<nr ; i++) {
			
			xabs = xa_mtbl_deux (l, k, j, i) ;
			yabs = ya_mtbl_deux (l, k, j, i) ;
			zabs = za_mtbl_deux (l, k, j, i) ;
			
			c1.get_mp()->convert_absolute(xabs, yabs, zabs, air, theta, phi) ;
			res2.set(l, k, j, i) = 
		(1-relax)*res2.set(l, k, j, i) +
	relax*(c2(l, k, j, i) - old_un.val_point(air, theta, phi)) ;	
	    }
	}
	
	// les coefficients ne sont plus a jour :
	res1.va.set_etat_c_qcq() ;
	res2.va.set_etat_c_qcq() ;
	// On raccord dans la zec :
	res1.raccord_externe (decrois, puiss, lmax) ;
	for (int i=0 ; i<decrois ; i++)
	    res1.dec_dzpuis() ;
	    
	res1.va.coef_i() ;
	res2.raccord_externe (decrois, puiss, lmax) ;
	for (int i=0 ; i<decrois ; i++)
	    res2.dec_dzpuis() ;
	res2.va.coef_i() ;
	
	// On regarde si on a converge :
	
	double erreur = 0 ;
	
	Tbl diff_un (diffrelmax(res1, old_un)) ;
	for (int i=1 ; i<nz_un-1 ; i++)
	    if (diff_un(i)>erreur)
		erreur = diff_un(i) ;
	
	Tbl diff_deux (diffrelmax(res2, old_deux)) ;
	for (int i=1 ; i<nz_deux-1 ; i++)
	    if (diff_deux(i)>erreur)
		erreur = diff_deux(i) ;
	
	if (flag == 1) 
	    cout << "Pas " << conte << " : erreur = " << erreur << endl ;
	if (erreur<=precision)
	    indic = -1 ;
	
	conte ++ ;
	if (conte > itemax)
	    indic = -1 ;
	}
    }
}
}
