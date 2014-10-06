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


char bhole_solve_phi_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:12:58  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:36:33  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2001/04/26  12:06:44  phil
 * *** empty log message ***
 *
 * Revision 2.2  2001/04/06  08:56:49  phil
 * *** empty log message ***
 *
 * Revision 2.1  2001/04/05  13:42:46  phil
 * *** empty log message ***
 *
 * Revision 2.0  2001/04/05  13:35:14  phil
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
#include "nbr_spx.h"
#include "tenseur.h"
#include "bhole.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"



void Bhole::init_bhole_phi () {
    
    Cmp auxi(mp) ;
    
    auxi = 1./2.-2*rayon/mp.r ;
    auxi.annule(0);
    auxi.set_dzpuis(0) ;
    n_auto = auxi;
    n_comp = 0 ; n_tot = 0;
    n_auto.set_std_base() ;
    n_auto.set().raccord(1) ;
    
    auxi = log (1+rayon/mp.r) ;
    auxi.annule(0);
    auxi.set_dzpuis(0) ;
    psi_auto = auxi;
    psi_comp = 0 ; psi_tot = 0;
    psi_auto.set_std_base() ;
    psi_auto.set().raccord(1) ;
    
    grad_n_tot = n_auto.gradient() ;
    grad_psi_tot = psi_auto.gradient() ;
    shift_auto.set_etat_zero() ;

    taij_auto.set_etat_zero();
    taij_comp.set_etat_zero();
    
    taij_tot.set_etat_zero() ;
    tkij_auto.set_etat_zero() ;
    tkij_tot.set_etat_zero();
    decouple.set_etat_zero() ;
}

void Bhole_binaire::solve_phi (double precision, double relax) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "-----------------------------------------------" << endl ;
    cout << "Resolution PSI" << endl ;
    
    Tenseur psi_un_old (hole1.psi_auto) ;
    Tenseur psi_deux_old (hole2.psi_auto) ;
    
    // Les sources totales, raccordees dans les zec
    Cmp source_un (-flat_scalar_prod(hole1.grad_psi_tot, 
			    hole1.psi_auto.gradient())()) ;
    source_un.std_base_scal() ;
    
    Cmp source_deux (-flat_scalar_prod(hole2.grad_psi_tot, 
			    hole2.psi_auto.gradient())())  ;
    source_deux.std_base_scal() ;
    
    // Les valeurs limites :
    Valeur lim_un (hole1.mp.get_mg()->get_angu()) ;
    lim_un = -0.5/hole1.rayon ;
    lim_un.std_base_scal() ;
    
    Valeur lim_deux (hole2.mp.get_mg()->get_angu()) ;
    lim_deux = -0.5/hole2.rayon ;
    lim_deux.std_base_scal() ;
     
    //On resout
    neumann_binaire (source_un, source_deux, lim_un, lim_deux, 
	hole1.psi_auto.set(), hole2.psi_auto.set(), 0, precision) ;
     
    hole1.psi_auto.set().raccord(1) ;
    hole2.psi_auto.set().raccord(1) ;
     
    //On verifie qu on a bien resolu :
    cout << diffrelmax (source_un, hole1.psi_auto().laplacien(4)) << endl ;
    cout << diffrelmax (source_deux, hole2.psi_auto().laplacien(4)) << endl ;
  
    // La relaxation :
    hole1.psi_auto.set() = relax*hole1.psi_auto() + (1-relax)*psi_un_old() ;
    hole2.psi_auto.set() = relax*hole2.psi_auto() + (1-relax)*psi_deux_old() ;
    
    hole1.fait_psi_comp (hole2) ;
    hole2.fait_psi_comp (hole1) ;
}

void Bhole_binaire::init_phi() {
    set_omega (0) ;
    hole1.init_bhole_phi() ;
    hole2.init_bhole_phi() ;
    
    hole1.fait_psi_comp(hole2) ;
    hole2.fait_psi_comp(hole1) ;
}


