/*
 *  Method of class Hor_isol to compute boundary conditions
 *
 *    (see file hor_isol.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
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

char bound_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2004/10/29 15:42:14  jl_jaramillo
 * Static shift boundary conbdition
 *
 * Revision 1.6  2004/10/01 16:46:51  f_limousin
 * Added a pure Dirichlet boundary condition
 *
 * Revision 1.5  2004/09/28 16:06:41  f_limousin
 * Correction of an error when taking the bases of the boundary
 * condition for the shift.
 *
 * Revision 1.4  2004/09/17 13:36:23  f_limousin
 * Add some new boundary conditions
 *
 * Revision 1.2  2004/09/09 16:53:49  f_limousin
 * Add the two lines $Id$Log: for CVS.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "graphique.h"
#include "utilitaires.h"




// Dirichlet boundary condition for Psi 
//-------------------------------------
//(This is implemented as a Dirichlet condition. It could be better to implemented as a Dirichlet one or
// as a mixed. It is analogous to the case for the lapse)
// ONE HAS TO GUARANTEE THAT BETA IS NOT ZERO, BUT IT IS PROPORTIONAL TO THE RADIAL VECTOR

Valeur Isol_hor::boundary_psi_Dir_evol(){

  const Map& map = ff.get_mp() ; 


  Scalar tmp = - 6 * contract(beta(), 0, psi().derive_cov(ff), 0) ;
  tmp = tmp / (contract(beta().derive_cov(ff), 0, 1) - nn() * trk() ) - 1 ;

  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  Valeur psi_bound (map.get_mg()->get_angu() )  ;
  
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;
			
  psi_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  psi_bound.std_base_scal() ;

  return psi_bound ;

}


// Neumann boundary condition for Psi 
//-------------------------------------
// ONE HAS TO GUARANTEE THAT BETA IS NOT ZERO, BUT IT IS PROPORTIONAL TO THE RADIAL VECTOR

Valeur Isol_hor::boundary_psi_Neu_evol(){

  const Map& map = ff.get_mp() ; 
  
  // Introduce 2-trace gamma tilde dot 
  Scalar tmp = - 1./ 6. * psi() * (beta().divergence(ff) - nn() * trk() ) 
    - beta()(2)* psi().derive_cov(ff)(2) - beta()(3)* psi().derive_cov(ff)(3) ;

  tmp = tmp / beta()(1) ;

  // in this case you don't have to substract any value
 
  Valeur psi_bound (map.get_mg()->get_angu() )  ;
  
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;
			
  psi_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  psi_bound.std_base_scal() ;
  
  return psi_bound ;

}


Valeur Isol_hor::boundary_psi_Dir_spat(){

  const Map& map = ff.get_mp() ; 
  
  Scalar tmp = psi() * psi() * psi() * trk() 
      - contract(k_dd(), 0, 1, tradial_vect_hor() * tradial_vect_hor(), 0, 1) 
      / psi()
      - 4.* contract(tradial_vect_hor(), 0, psi().derive_cov(ff), 0) ;

  tmp = tmp / (tradial_vect_hor().divergence(ff)) - 1. ;

  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  
 
  Valeur psi_bound (map.get_mg()->get_angu() )  ;
  
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;
			
  psi_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0)  ;

  psi_bound.std_base_scal() ;
  
  return psi_bound ;

}

Valeur Isol_hor::boundary_psi_Neu_spat(){

  const Map& map = ff.get_mp() ; 
  
  Scalar tmp = psi() * psi() * psi() * trk() 
      - contract(k_dd(), 0, 1, tradial_vect_hor() * tradial_vect_hor(), 0, 1) 
      / psi()
      - psi() * tradial_vect_hor().divergence(ff) 
      - 4 * ( tradial_vect_hor()(2) * psi().derive_cov(ff)(2) 
	      + tradial_vect_hor()(3) * psi().derive_cov(ff)(3) ) ;

  tmp = tmp / (4 * tradial_vect_hor()(1)) ;

  // in this case you don't have to substract any value
 
  Valeur psi_bound (map.get_mg()->get_angu() )  ;
  
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;
			
  psi_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  psi_bound.std_base_scal() ;
  
  return psi_bound ;

}



// Dirichlet boundary condition on nn using the extrinsic curvature
// (No time evolution taken into account! Make this)
//--------------------------------------------------------------------------
Valeur Isol_hor::boundary_nn_Dir_kk(){

  const Map& map = ff.get_mp() ;

  Scalar tmp(map) ;

  Scalar kk_rr = contract( radial_vect_hor() * radial_vect_hor(), 0, 1
			   , k_dd(), 0, 1 ) ;

  Scalar k_kerr (map) ;
  k_kerr = kappa_hor() ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;

  tmp = k_kerr - contract(radial_vect_hor(), 0, nn().derive_cov(ff), 0) ;

  tmp = - tmp / kk_rr - 1;

  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur nn_bound (map.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}



// Neumann boundary condition on nn using the extrinsic curvature
// (No time evolutuon taken into account! Make this)
//--------------------------------------------------------------------------
Valeur Isol_hor::boundary_nn_Neu_kk() {
  
  const Map& map = ff.get_mp() ;
  
  const Vector& dnn= nn().derive_cov(ff) ;

  Scalar kk_rr = contract( radial_vect_hor() * radial_vect_hor(), 0, 1
			   , k_dd(), 0, 1 ) ; 

  Scalar k_kerr (map) ;
  k_kerr = kappa_hor() ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;

  Scalar tmp = k_kerr + nn() * kk_rr 
             - radial_vect_hor()(2) * dnn(2) - radial_vect_hor()(3) * dnn(3)  ;

  tmp = tmp / radial_vect_hor()(1) ;

  // in this case you don't have to substract any value
 
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur nn_bound (map.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}



Valeur Isol_hor::boundary_nn_Dir_eff(double aa){

  const Map& map = ff.get_mp() ;

  Scalar tmp(map) ;

  tmp = - aa * nn().derive_cov(ff)(1) ;
  tmp.dec_dzpuis(2) ;
  tmp = tmp - 1. ;
  
  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur nn_bound (map.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}



Valeur Isol_hor::boundary_nn_Neu_eff(double aa) {
  
  const Map& map = ff.get_mp() ;
  
  Scalar tmp = - aa * nn() ;

  // in this case you don't have to substract any value
 
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur nn_bound (map.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}


Valeur Isol_hor::boundary_nn_Dir(double aa){

  const Map& map = ff.get_mp() ;

  Scalar tmp(map) ;
  tmp = aa - 1 ;
  
  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur nn_bound (map.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}



// Component r of boundary value of beta (using expression in terms of radial vector)
//--------------------------------------
Valeur Isol_hor:: boundary_beta_r(){

  const Map& map = ff.get_mp() ; 

  Scalar tmp (map) ;

  tmp = nn() * radial_vect_hor()(1) ;
 
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur bnd_beta_r (map.get_mg()->get_angu()) ;
    
  bnd_beta_r = 1 ;   // Why is it necessary this and what it is actually doing?
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_r.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  bnd_beta_r.std_base_scal() ;
  
  return  bnd_beta_r ;


}



// Component theta of boundary value of beta (using expression in terms of radial vector)
//------------------------------------------
Valeur Isol_hor::boundary_beta_theta(){
  
  const Map& map = ff.get_mp() ; 

  Scalar tmp(map) ;  
  
  tmp = nn() * radial_vect_hor()(2) ;

  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur bnd_beta_theta (map.get_mg()->get_angu()) ;
    
  bnd_beta_theta = 1 ;   // Why is it necessary this and what it is actually doing?
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_theta.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  bnd_beta_theta.std_base_scal() ;
  
  return  bnd_beta_theta ;


}
 
// Component phi of boundary value of beta (using expression in terms of radial vector) 
//-------------------------------------------------------------------------------------
Valeur Isol_hor::boundary_beta_phi(){

  const Map& map = ff.get_mp() ; 

  Scalar tmp (map) ;

  Scalar vel_ang(map) ;
  vel_ang = omega_hor() ;
  vel_ang.std_spectral_base() ;
  vel_ang.mult_rsint() ;

  tmp = nn() * radial_vect_hor()(3)  -  vel_ang ;

  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur bnd_beta_phi (map.get_mg()->get_angu()) ;
    
  bnd_beta_phi = 1 ; // Why is it necessary this and what it is actually doing?
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_phi.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  bnd_beta_phi.std_base_scal() ;
  
  return  bnd_beta_phi ;

}




// Component x of boundary value of beta (using expression in terms of radial vector)
//--------------------------------------
Valeur Isol_hor:: boundary_beta_x(double velang){

  const Map& map = ff.get_mp() ; 
  
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  //Isol_hor boundary conditions
  
  Valeur lim_x (map.get_mg()->get_angu()) ;
    
  lim_x = 1 ;   // Why is it necessary this and what it is actually doing?
  
  Scalar beta_x = beta_bound_cart(velang)(1) ;

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_x.set(0, k, j, 0) = beta_x.val_grid_point(1, k, j, 0) ;
  
  lim_x.set_base(beta_x.get_spectral_va().get_base()) ;

  return  lim_x ;


}


// Component y of boundary value of beta (using expression in terms of radial vector)
//--------------------------------------
Valeur Isol_hor:: boundary_beta_y(double velang){
  
  const Map& map = ff.get_mp() ;
 
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
  
  Valeur lim_y (map.get_mg()->get_angu()) ;
    
  lim_y = 1 ;   // Why is it necessary this and what it is actually doing?
 
  Scalar beta_y = beta_bound_cart(velang)(2) ;

  for (int k=0 ; k<nnp ; k++)
      for (int j=0 ; j<nnt ; j++)
	  lim_y.set(0, k, j, 0) = beta_y.val_grid_point(1, k, j, 0) ;

  lim_y.set_base(beta_y.get_spectral_va().get_base()) ;

  return  lim_y ;
}


// Component z of boundary value of beta (using expression in terms of radial vector)
//--------------------------------------
Valeur Isol_hor:: boundary_beta_z(double velang){

  const Map& map = ff.get_mp() ; 

  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
 
  Valeur lim_z (map.get_mg()->get_angu()) ;
    
  lim_z = 1 ;   // Why is it necessary this and what it is actually doing?
  
  Scalar beta_z = beta_bound_cart(velang)(3) ;

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_z.set(0, k, j, 0) = beta_z.val_grid_point(1, k, j, 0) ;
 
  lim_z.set_base(beta_z.get_spectral_va().get_base()) ;

  return  lim_z ;
}


Vector Isol_hor::beta_bound_cart(double velang) {

  const Map& map = ff.get_mp() ; 

  Vector tmp_vect = nn() * radial_vect_hor() ;

  Scalar vel_ang (map) ;  
  vel_ang = velang ;
  vel_ang.std_spectral_base() ;
  vel_ang.mult_rsint() ;

  tmp_vect.set(3) = tmp_vect(3) - vel_ang ;

  tmp_vect.change_triad(map.get_bvect_cart() ) ;
  
  return tmp_vect ;

}
