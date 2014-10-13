/*
 *   Copyright (c) 2003 Keisuke Taniguchi
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char pseudo_misner_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:52:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:02  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2007/04/24 20:13:53  f_limousin
 * Implementation of Dirichlet and Neumann BC for the lapse
 *
 * Revision 1.1  2005/09/09 09:00:02  p_grandclement
 * add pseudo_misner
 *

 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "bhole.h"
#include "nbr_spx.h"
#include "et_bin_nsbh.h"
#include "etoile.h"
#include "param.h"
#include "bin_ns_bh.h"

#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
void Bin_ns_bh::pseudo_misner (int& ite, int itemax, double relax, 
			       double precis, int bound_nn, double lim_nn) {
			    
    using namespace Unites ;
 
    // Parameters for the function Map_et::poisson for n_auto
    // ------------------------------------------------------
    Cmp source_n_prev (star.get_mp()) ;
    source_n_prev.set_etat_zero() ;
    
    Param par_poisson1 ;
    par_poisson1.add_int(itemax,  0) ;  // maximum number of iterations
    par_poisson1.add_double(relax,  0) ; // relaxation parameter
    par_poisson1.add_double(precis, 1) ; // required precision
    par_poisson1.add_int_mod(ite, 0) ; // number of iterations actually used
    par_poisson1.add_cmp_mod(source_n_prev) ;
    
    // Parameters for the function Map_et::poisson for confpsi_auto
    // ------------------------------------------------------------
    Cmp source_psi_prev (star.get_mp()) ;
    source_psi_prev.set_etat_zero() ;
    
    Param par_poisson2 ;
    par_poisson2.add_int(itemax,  0) ;  // maximum number of iterations
    par_poisson2.add_double(relax,  0) ; // relaxation parameter
    par_poisson2.add_double(precis, 1) ; // required precision
    par_poisson2.add_int_mod(ite, 0) ; // number of iterations actually used
    par_poisson2.add_cmp_mod(source_psi_prev) ;
  

    bool loop = true ;
    
    //=========================================================================
    // 			Start of iteration
    //=========================================================================
    double erreur ;
    int itere =  1 ;
    
    while (loop) {

    	Tenseur n_auto_old (star.n_auto()) ;
	Tenseur psi_auto_old (star.confpsi_auto()) ;
	Tenseur n_auto_hole (hole.n_auto()) ;
	
	Tenseur confpsi_q = pow(star.confpsi, 4.) ;
	Tenseur confpsi_c = pow(star.confpsi, 5.) ;
		
	// Lapse star
	Tenseur source_n (qpig * star.nnn * confpsi_q * (star.ener_euler + star.s_euler)
		- 2.*flat_scalar_prod((star.d_confpsi_auto+star.d_confpsi_comp), star.d_n_auto) / star.confpsi) ;
	source_n.set_std_base() ;
	source_n().poisson(par_poisson1, star.n_auto.set()) ;
	star.n_auto.set() = star.n_auto() + 0.5 ; 
	star.n_auto.set() = relax * star.n_auto() + (1-relax)*n_auto_old() ;
        erreur = max(diffrelmax(star.n_auto(), n_auto_old())) ;
	
	// Psi star
	Tenseur source_psi (-0.5 * qpig * confpsi_c * star.ener_euler) ;
	source_psi.set_std_base() ;
	source_psi().poisson(par_poisson2, star.confpsi_auto.set()) ;
	star.confpsi_auto.set() = star.confpsi_auto() + 0.5 ;
	star.confpsi_auto.set() = relax*star.confpsi_auto() + (1-relax)*psi_auto_old() ;
	
	// Trou noir :
	hole.update_metric (star) ;
	hole.solve_lapse_with_ns (relax, bound_nn, lim_nn) ;
	//erreur = max(diffrelmax(hole.n_auto(), n_auto_hole())) ;
	
	hole.solve_psi_with_ns (relax) ;
		
	star.update_metric (hole) ;
	star.update_metric_der_comp(hole) ;
		
	star.equation_of_state() ;
        star.kinematics(get_omega(), get_x_axe()) ;
	star.fait_d_psi() ;
        star.hydro_euler() ;
	 
	cout << "Step " << itere << " " << erreur << endl ;
	
	if ((itere==itemax) || (erreur<precis))
	    loop = false ;
	itere ++ ;
    }
    
  ite = itere ;
}
}
