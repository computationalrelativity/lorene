/*
 * Reads a black hole with an isolated horizon configuration 
 *
 */

/*
 *   Copyright (c) 2005 Francois Limousin
 *                      Jose Luis Jaramillo
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

char lit_ih_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/06 15:09:45  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2005/09/12 12:34:09  f_limousin
 * Compilation Warning - Change of convention for the angular velocity
 * Add Berlin boundary condition in the case of binary horizons.
 *
 * Revision 1.3  2005/03/09 10:35:28  f_limousin
 * Drawing of K_{ij}s^is^j and new initialisation of A^{ij} thanks to
 * the function update_aa().
 *
 * Revision 1.2  2005/03/04 18:24:09  jl_jaramillo
 * Extrinsic curvature graphical output.
 *
 * Revision 1.1  2005/03/03 10:19:47  f_limousin
 * First version
 *
 * 
 * $Header$
 *
 */

//standard
#include <cstdlib>
#include <cmath>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "param.h"
#include "coord.h"
#include "scalar.h"
#include "cmp.h"
#include "tensor.h"
#include "tenseur.h"
#include "isol_hor.h"
#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"


int main(int argc, char** argv) {

  using namespace Unites ;
  
    if (argc <2) {
	cout <<" Passer nom du ficher en arguments SVP !" << endl ;
	abort() ;
    }
  
    char* name_fich = argv[1] ;
  
    FILE* fich = fopen(name_fich, "r") ;    
    Mg3d grid (fich) ;
    Map_af mp (grid, fich) ;
    Isol_hor isolhor (mp, fich, true) ;
    fclose(fich) ;
    
    isolhor.update_aa() ;

    // PLOTS 
    // -------

    // Graphic output of the different fields
    //---------------------------------------

    Vector beta (isolhor.beta()) ;
    Scalar bb (contract(beta, 0, isolhor.tgam().radial_vect()
			.down(0, isolhor.tgam()), 0)
	       *isolhor.psi()*isolhor.psi()) ;

    des_meridian(isolhor.nn(), 1.00001, 10, "nn", 0) ; 
    arrete() ;
    des_profile(isolhor.psi(), 1.00001, 10, 1., 1., "psi") ;
    des_profile(bb, 1.00001, 10, M_PI/2., 0., "b in Pi/2, 0") ;
    des_profile(bb, 1.00001, 10, M_PI/2., M_PI, "b in Pi/2, Pi") ;
    des_profile(isolhor.beta()(3), 1.00001, 10, 
		M_PI/2., 0., "beta_phi en pi/2") ;

    // Definition of the surface

    Cmp surface (mp) ;
    surface = pow(mp.r, 2.)-pow(isolhor.get_radius(), 2.) ;
    surface.annule(grid.get_nzone()-1) ;
    surface.std_base_scal() ;
    
    // Shift
    
    Tenseur shift (mp, 1, CON, mp.get_bvect_cart()) ;
    shift.set_etat_qcq() ;
    beta.change_triad(mp.get_bvect_cart()) ;
    Cmp shift_x (beta(1)) ;
    Cmp shift_y (beta(2)) ;
    Cmp shift_z (beta(3)) ;
    shift.set(0) = shift_x ;
    shift.set(1) = shift_y ;
    shift.set(2) = shift_z ;

    des_coupe_vect_z(shift, 0, -1, 0.7, 3, "shift", &surface) ;
    des_coupe_vect_x(shift, 0, -1, 0.7, 3, "shift", &surface) ;

    Vector beta_rad (mp, CON, mp.get_bvect_spher()) ;
    beta_rad = isolhor.b_tilde() * isolhor.tgam().radial_vect() ;
    beta_rad.change_triad(mp.get_bvect_cart()) ;
    shift_x = beta_rad(1) ;
    shift_y = beta_rad(2) ;
    shift_z = beta_rad(3) ;
    shift.set(0) = shift_x ;
    shift.set(1) = shift_y ;
    shift.set(2) = shift_z ;

    des_coupe_vect_z(shift, 0, -1, 0.7, 3, "radial shift ", &surface) ;
    des_coupe_vect_x(shift, 0, -1, 0.7, 3, "radial shift", &surface) ;

    Vector beta_ang (mp, CON, mp.get_bvect_spher()) ;
    beta_ang = beta - beta_rad ;
    beta_ang.change_triad(mp.get_bvect_cart()) ;
    shift_x = beta_ang(1) ;
    shift_y = beta_ang(2) ;
    shift_z = beta_ang(3) ;
    shift.set(0) = shift_x ;
    shift.set(1) = shift_y ;
    shift.set(2) = shift_z ;

    des_coupe_vect_z(shift, 0, -1, 0.7, 3, "shift angular", &surface) ;
    des_coupe_vect_x(shift, 0, -1, 0.7, 3, "shift angular", &surface) ;

    // \gamma^{rr} and \gamma^{rp}

    Cmp gam_rr (isolhor.gam_uu()(1,1)) ;
    Cmp gam_rp (isolhor.gam_uu()(3,1)) ;

    des_coupe_z(gam_rr, 0, -5, 5, -5, 5, "gam_rr", &surface) ;
    des_coupe_x(gam_rr, 0, -5, 5, -5, 5, "gam_rr", &surface) ;
    des_coupe_z(gam_rp, 0, -5, 5, -5, 5, "gam_rp", &surface) ;
    des_coupe_x(gam_rp, 0, -5, 5, -5, 5, "gam_rp", &surface) ;
    
    // K^{rr}, K^{rt} and K^{rp}
    
    Vector pp(contract(isolhor.k_dd(), 0,isolhor.tgam().radial_vect() , 0)/isolhor.psi()/isolhor.psi()) ;
    Cmp ksp (pp(3)) ;
    des_coupe_z(ksp, 0, -5, 5, -5, 5, "K^sp", &surface) ;
    
    pp.change_triad(mp.get_bvect_cart()) ;
    Cmp ksx (pp(1)) ;
    Cmp ksy (pp(2)) ;
    Cmp ksz (pp(3)) ;
    
    des_coupe_z(ksx, 0, -5, 5, -5, 5, "K^sx", &surface) ;
    des_coupe_z(ksy, 0, -5, 5, -5, 5, "K^sy", &surface) ;
    des_coupe_z(ksz, 0, -5, 5, -5, 5, "K^sz", &surface) ; 

    // Lapse

    Cmp lapse (isolhor.nn()) ;
    des_coupe_z(lapse, 0, -5, 5, -5, 5, "lapse", &surface) ;
    des_coupe_x(lapse, 0, -5, 5, -5, 5, "lapse", &surface) ;
    
    // Psi
    
    Cmp psi (isolhor.psi()) ;
    des_coupe_z(psi, 0, -5, 5, -5, 5, "psi", &surface) ;
    des_coupe_x(psi, 0, -5, 5, -5, 5, "psi", &surface) ;
    



    return 1 ; 
}
