/*
 *  Method of class single_hor to compute boundary conditions
 *
 *    (see file isol_hor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
 *                       Francois Limousin
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

char single_bound_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:01  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/04/13 15:28:35  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cassert>

// Lorene headers
#include "time_slice.h"
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "graphique.h"
#include "utilitaires.h"
#include "param.h"



namespace Lorene {
const Valeur Single_hor::boundary_psi_app_hor()const {

  Metric flat (mp.flat_met_spher()) ;
  Vector temp (mp, CON, mp.get_bvect_spher()) ;
  temp.set(1) = 1.;
  temp.set(2) = 0.;
  temp.set(3) = 0.;
  temp.std_spectral_base() ;

  Scalar tmp = psi * psi * psi * trK 
    - contract(get_k_dd(),0, 1, tgam.radial_vect() * tgam.radial_vect(), 0, 1) 
    / psi
    - psi * tgam.radial_vect().divergence(ff) 
    - 4 * ( tgam.radial_vect()(2) * psi.derive_cov(ff)(2) 
	    + tgam.radial_vect()(3) * psi.derive_cov(ff)(3) ) ;
  
  tmp = tmp / (4 * tgam.radial_vect()(1)) ;

  // in this case you don't have to substract any value
 
  Valeur psi_bound (mp.get_mg()->get_angu() )  ;
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;
			
  psi_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  psi_bound.std_base_scal() ;
  
  return psi_bound ;

}

const Valeur Single_hor::boundary_nn_Neu(double cc)const  {
  
  double rho = 1. ; // 1 is the standart case;

  Scalar tmp = - cc * nn ;
  //  Scalar tmp = - nn()/psi()*psi().dsdr() ;

  // in this case you don't have to substract any value
  tmp += (rho - 1) * tgam.radial_vect()(1)  * dn(1)  ;
  tmp = tmp / (rho * tgam.radial_vect()(1))  ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur nn_bound (mp.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}


const Valeur Single_hor::boundary_nn_Dir(double cc)const {

  Scalar rho(mp);
  rho = 0. ;  // 0 is the standard case

  Scalar tmp(mp) ;
  tmp = cc ; 
  
 
  //tmp = 1./(2*psi()) ;
  //  tmp = - psi() * nn().dsdr() / (psi().dsdr())  ;

  // We  have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  tmp = (tmp + rho * nn)/(1 + rho) ;

  tmp = tmp - 1 ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur nn_bound (mp.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;  // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}

// Component x of boundary value of beta
//--------------------------------------

const Valeur Single_hor:: boundary_beta_x(double om_orb, 
					  double om_loc)const {

    // Les alignemenents pour le signe des CL.
    double orientation = mp.get_rot_phi() ;
    assert ((orientation == 0) || (orientation == M_PI)) ;
    int aligne = (orientation == 0) ? 1 : -1 ;
   
    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    Vector tmp_vect = nn * get_gam().radial_vect() ;
    tmp_vect.change_triad(mp.get_bvect_cart() ) ;

    //Isol_hor boundary conditions
  
    Valeur lim_x (mp.get_mg()->get_angu()) ;
    
    lim_x = 1 ;  // Juste pour affecter dans espace des configs ;
  
    Mtbl ya_mtbl (mp.get_mg()) ;
    ya_mtbl.set_etat_qcq() ;
    ya_mtbl = mp.ya ;

    Mtbl yy_mtbl (mp.get_mg()) ;
    yy_mtbl.set_etat_qcq() ;
    yy_mtbl = mp.y ;

    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	  lim_x.set(0, k, j, 0) = aligne * om_orb * ya_mtbl(1, k, j, 0)
	    + (om_loc-om_orb)* yy_mtbl(1, k, j, 0)
	    + tmp_vect(1).val_grid_point(1, k, j, 0) ;
    
  lim_x.set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;

  return  lim_x ;
}


// Component y of boundary value of beta 
//--------------------------------------

const Valeur Single_hor:: boundary_beta_y(double om_orb, 
					  double om_loc)const {

    // Les alignemenents pour le signe des CL.
    double orientation = mp.get_rot_phi() ;
    assert ((orientation == 0) || (orientation == M_PI)) ;
    int aligne = (orientation == 0) ? 1 : -1 ;
    

    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    Vector tmp_vect = nn * get_gam().radial_vect() ;
    tmp_vect.change_triad(mp.get_bvect_cart() ) ;

    //Isol_hor boundary conditions
  
    Valeur lim_y (mp.get_mg()->get_angu()) ;
    
    lim_y = 1 ;  // Juste pour affecter dans espace des configs ;
  
    Mtbl xa_mtbl (mp.get_mg()) ;
    xa_mtbl.set_etat_qcq() ;
    xa_mtbl = mp.xa ;

    Mtbl xx_mtbl (mp.get_mg()) ;
    xx_mtbl.set_etat_qcq() ;
    xx_mtbl = mp.x ;
  
    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	  lim_y.set(0, k, j, 0) = - aligne *om_orb * xa_mtbl(1, k, j, 0) -
	    (om_loc-om_orb)*xx_mtbl(1, k, j, 0)
	    + tmp_vect(2).val_grid_point(1, k, j, 0) ;
    
  lim_y.set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;

  return  lim_y ;
}

// Component z of boundary value of beta 
//--------------------------------------

const Valeur Single_hor:: boundary_beta_z()const {
    
    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    Vector tmp_vect = nn * get_gam().radial_vect() ;
    tmp_vect.change_triad(mp.get_bvect_cart() ) ;

    //Isol_hor boundary conditions
  
    Valeur lim_z (mp.get_mg()->get_angu()) ;
    
    lim_z = 1 ;  // Juste pour affecter dans espace des configs ;
  
    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	    lim_z.set(0, k, j, 0) = tmp_vect(3).val_grid_point(1, k, j, 0) ;
    
  lim_z.set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;

  return  lim_z ;
}

}
