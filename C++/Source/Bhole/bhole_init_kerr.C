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


char bhole_init_kerr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:12:58  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:36:32  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/12/14  10:45:20  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.0  2000/10/20  09:18:56  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
//standard
#include <cstdlib>
#include <cmath>

// Lorene
#include "tenseur.h"
#include "bhole.h"

void Bhole::init_kerr (double masse, double moment) {
    
    // On verifie si le rayon est bien calcule 
    assert (rayon == sqrt (masse*masse-moment*moment)/2.) ;
    
    // Valeur de omega :
    omega = moment/2/masse/(masse+sqrt(masse*masse-moment*moment)) ;
        
    // Calcul de R :
    Mtbl grand_r (mp.get_mg()) ;
    grand_r = mp.r + (masse*masse-moment*moment)/4/mp.r + masse ;
   
    // Calcul de sigma :
    Mtbl sigma (mp.get_mg()) ;
    sigma = moment*moment*mp.cost*mp.cost + grand_r*grand_r ;
    
    // Calcul de grand_a :
    Cmp grand_a (mp) ;
    grand_a = 1 + 2*masse/mp.r + 
     (3*masse*masse+moment*moment*mp.cost*mp.cost)/2/mp.r/mp.r
     + (2*masse*rayon*rayon)/pow(mp.r, 3.) + pow(rayon/mp.r, 4.) ;
    grand_a.set_val_inf(1) ;
    grand_a.std_base_scal() ;
    grand_a.raccord(1) ;
   
    // Calcul de n_phi :
    Cmp n_phi(mp) ;
    n_phi = (2*moment*masse*grand_r) / (sigma*(grand_r*grand_r+moment*moment)
	+ 2*moment*moment*masse*grand_r*mp.sint*mp.sint) ;
    n_phi.annule(0) ;
    n_phi.set_val_inf (0) ;
    n_phi.std_base_scal() ;
    
    // Calcul de N :
    Cmp carre (mp) ;
    carre = 1-(2*masse*grand_r)/sigma + (4*moment*moment*masse*masse
	*grand_r*grand_r*mp.sint*mp.sint)/
	(sigma*sigma*(grand_r*grand_r+moment*moment)+2*moment*moment*sigma*masse*
	grand_r*mp.sint*mp.sint) ;
    carre.set_val_inf(1) ;
    carre.set_val_hor(0, 1) ;
    carre.std_base_scal() ;
    carre.annule(0) ;
    
    n_auto.set_etat_qcq() ;
    n_auto.set() = sqrt(carre) ;
    n_auto.set().set_dzpuis(0) ;
    n_auto.set_std_base() ;
    n_auto.set().raccord(1) ;
    
    // Calcul de psi :
    psi_auto.set_etat_qcq() ;
    psi_auto.set() = pow(grand_a, 0.25) ;
    psi_auto.set().set_dzpuis(0) ;
    psi_auto.set_std_base() ;
    psi_auto.set().raccord(1) ;
    
    // Calcul du shift :
    shift_auto.set_etat_qcq() ;
    shift_auto.set_std_base() ;
    Valeur auxi (mp.get_mg()) ;
    auxi = n_phi.va.mult_st().mult_sp() ;
    shift_auto.set(0) = auxi ;
    auxi = -n_phi.va.mult_st().mult_cp() ;
    shift_auto.set(1) = auxi ;
    shift_auto.set(2).set_etat_zero() ;
    
    shift_auto.inc_dzpuis() ;
    
    for (int i=0 ; i<2 ; i++) {
	shift_auto.set(i).mult_r() ;
	shift_auto.set(i).raccord(1) ;
	assert (shift_auto(i).check_dzpuis (0)) ;
	}
}
