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
 * Revision 1.12  2004/12/31 15:34:37  f_limousin
 * Modifications to avoid warnings
 *
 * Revision 1.11  2004/12/22 18:15:16  f_limousin
 * Many different changes.
 *
 * Revision 1.10  2004/11/24 19:30:58  jl_jaramillo
 * Berlin boundary conditions  vv_bound_cart
 *
 * Revision 1.9  2004/11/18 09:49:44  jl_jaramillo
 * Some new conditions for the shift (Neumann + Dirichlet)
 *
 * Revision 1.8  2004/11/05 10:52:26  f_limousin
 * Replace double aa by double cc in argument of boundary_beta_x
 * boundary_beta_y and boundary_beta_z to avoid warnings.
 *
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

    Scalar tmp = - 6 * contract(beta(), 0, psi().derive_cov(ff), 0) ;
  tmp = tmp / (contract(beta().derive_cov(ff), 0, 1) - nn() * trk() ) - 1 ;

  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

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


// Neumann boundary condition for Psi 
//-------------------------------------
// ONE HAS TO GUARANTEE THAT BETA IS NOT ZERO, BUT IT IS PROPORTIONAL TO THE RADIAL VECTOR

Valeur Isol_hor::boundary_psi_Neu_evol(){

  // Introduce 2-trace gamma tilde dot 
  Scalar tmp = - 1./ 6. * psi() * (beta().divergence(ff) - nn() * trk() ) 
    - beta()(2)* psi().derive_cov(ff)(2) - beta()(3)* psi().derive_cov(ff)(3) ;

  tmp = tmp / beta()(1) ;

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


Valeur Isol_hor::boundary_psi_Dir_spat(){

    Scalar tmp = psi() * psi() * psi() * trk() 
      - contract(k_dd(), 0, 1, tradial_vect_hor() * tradial_vect_hor(), 0, 1) 
      / psi()
      - 4.* contract(tradial_vect_hor(), 0, psi().derive_cov(ff), 0) ;

  tmp = tmp / (tradial_vect_hor().divergence(ff)) - 1. ;

  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  
 
  Valeur psi_bound (mp.get_mg()->get_angu() )  ;
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;
			
  psi_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0)  ;

  psi_bound.std_base_scal() ;
  
  return psi_bound ;

}

Valeur Isol_hor::boundary_psi_Neu_spat(){

    Scalar tmp = psi() * psi() * psi() * trk() 
      - contract(k_dd(), 0, 1, tradial_vect_hor() * tradial_vect_hor(), 0, 1) 
      / psi()
      - psi() * tradial_vect_hor().divergence(ff) 
      - 4 * ( tradial_vect_hor()(2) * psi().derive_cov(ff)(2) 
	      + tradial_vect_hor()(3) * psi().derive_cov(ff)(3) ) ;

  tmp = tmp / (4 * tradial_vect_hor()(1)) ;

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



// Dirichlet boundary condition on nn using the extrinsic curvature
// (No time evolution taken into account! Make this)
//--------------------------------------------------------------------------
Valeur Isol_hor::boundary_nn_Dir_kk(){

  Scalar tmp(mp) ;

  Scalar kk_rr = contract( radial_vect_hor() * radial_vect_hor(), 0, 1
			   , k_dd(), 0, 1 ) ;

  Scalar k_kerr (mp) ;
  k_kerr = kappa_hor() ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;

  tmp = k_kerr - contract(radial_vect_hor(), 0, nn().derive_cov(ff), 0) ;

  tmp = - tmp / kk_rr - 1;

  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur nn_bound (mp.get_mg()->get_angu()) ;
    
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
  
  const Vector& dnnn = nn().derive_cov(ff) ;

  Scalar kk_rr = contract( radial_vect_hor() * radial_vect_hor(), 0, 1
			   , k_dd(), 0, 1 ) ; 

  Scalar k_kerr (mp) ;
  k_kerr = kappa_hor() ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;

  Scalar tmp = k_kerr + nn() * kk_rr 
             - radial_vect_hor()(2) * dnnn(2) - radial_vect_hor()(3) * dnnn(3)  ;

  tmp = tmp / radial_vect_hor()(1) ;

  // in this case you don't have to substract any value
 
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur nn_bound (mp.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}



Valeur Isol_hor::boundary_nn_Dir_eff(double cc){

  Scalar tmp(mp) ;

  tmp = - cc * nn().derive_cov(ff)(1) ;
  tmp.dec_dzpuis(2) ;
  tmp = tmp - 1. ;
  
  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur nn_bound (mp.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}



Valeur Isol_hor::boundary_nn_Neu_eff(double cc) {
  
  Scalar tmp = - cc * nn() ;

  // in this case you don't have to substract any value
 
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur nn_bound (mp.get_mg()->get_angu()) ;
    
  nn_bound = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}


Valeur Isol_hor::boundary_nn_Dir(double cc){

  Scalar tmp(mp) ;
  tmp = cc - 1 ;
  
  // We  have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur nn_bound (mp.get_mg()->get_angu()) ;
    
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

  Scalar tmp (mp) ;

  tmp = nn() * radial_vect_hor()(1) ;
 
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_r (mp.get_mg()->get_angu()) ;
    
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
  
  Scalar tmp(mp) ;  
  
  tmp = nn() * radial_vect_hor()(2) ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_theta (mp.get_mg()->get_angu()) ;
    
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

  Scalar tmp (mp) ;

  Scalar vel_ang(mp) ;
  vel_ang = omega_hor() ;
  vel_ang.std_spectral_base() ;
  vel_ang.mult_rsint() ;

  tmp = nn() * radial_vect_hor()(3)  -  vel_ang ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_phi (mp.get_mg()->get_angu()) ;
    
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
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  //Isol_hor boundary conditions
  
  Valeur lim_x (mp.get_mg()->get_angu()) ;
    
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
    
    int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
  
  Valeur lim_y (mp.get_mg()->get_angu()) ;
    
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

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
 
  Valeur lim_z (mp.get_mg()->get_angu()) ;
    
  lim_z = 1 ;   // Why is it necessary this and what it is actually doing?
  
  Scalar beta_z = beta_bound_cart(velang)(3) ;

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_z.set(0, k, j, 0) = beta_z.val_grid_point(1, k, j, 0) ;
 
  lim_z.set_base(beta_z.get_spectral_va().get_base()) ;

  return  lim_z ;
}


Vector Isol_hor::beta_bound_cart(double velang) {

  Vector tmp_vect = nn() * radial_vect_hor() ;

 
  Scalar vel_ang (mp) ;  
  vel_ang = velang ;
  vel_ang.std_spectral_base() ;
  vel_ang.mult_rsint() ;

  tmp_vect.set(3) = tmp_vect(3) - vel_ang ;

  tmp_vect.change_triad(mp.get_bvect_cart() ) ;
  
  return tmp_vect ;

}


// Neumann boundary condition for b_tilde 
//---------------------------------------
// ONE HAS TO GUARANTEE THAT BETA IS NOT ZERO, BUT IT IS PROPORTIONAL TO THE RADIAL VECTOR

Valeur Isol_hor::boundary_b_tilde_Neu(){
  
  // Introduce 2-trace gamma tilde dot

  Vector s_tilde = met_gamt.radial_vect() ;

  Scalar hh_tilde = contract(s_tilde.derive_cov(met_gamt), 0, 1) ;

  //des_profile(hh_tilde, 1.00001, 10, M_PI/2., 0., "H_tilde") ;

  Scalar tmp (mp) ;

  tmp = + b_tilde() * hh_tilde - 2 * ( s_tilde(2) * b_tilde().derive_cov(ff)(2)
				     + s_tilde(3) * b_tilde().derive_cov(ff)(3) ) ;
  
  Scalar constant (mp) ;
  constant = 0. ;
  constant.std_spectral_base() ;
  constant.inc_dzpuis(2) ;

  tmp += constant ;
    
    

  tmp = tmp / (2 *  s_tilde(1) ) ;


  // in this case you don't have to substract any value
 
  Valeur b_tilde_bound (mp.get_mg()->get_angu() )  ;
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;
			
  b_tilde_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      b_tilde_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  b_tilde_bound.std_base_scal() ;
  
  return b_tilde_bound ;

}


Valeur Isol_hor::boundary_b_tilde_Dir(){

  Vector s_tilde = met_gamt.radial_vect() ;

  Scalar hh_tilde = contract(s_tilde.derive_cov(met_gamt), 0, 1) ;

  Scalar tmp = 2 * contract (s_tilde, 0, b_tilde().derive_cov(ff) , 0) ; 

  Scalar constant (mp) ;
  constant = -1. ;
  constant.std_spectral_base() ;
  constant.inc_dzpuis(2) ;

  tmp -= constant ;
    
  tmp = tmp / hh_tilde   ;





  //  des_profile(tmp, 1.00001, 10, M_PI/2., 0., "tmp") ;


  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  Valeur b_tilde_bound (mp.get_mg()->get_angu() )  ;
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;
			
  b_tilde_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      b_tilde_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  b_tilde_bound.std_base_scal() ;

  return b_tilde_bound ;

}


Vector Isol_hor::vv_bound_cart(double velang) {

  // Preliminaries
  //--------------

  Vector s_tilde =  tradial_vect_hor() ;

  //cout<<"s_tilde"<<endl;
  //  cout<<max(s_tilde)<<endl;

  
  Scalar hh_tilde = contract(s_tilde.derive_cov(met_gamt), 0, 1) ;
  hh_tilde.dec_dzpuis(2) ;

  Vector tmp_vect = b_tilde() * s_tilde ;
  

  Scalar vel_ang (mp) ;  
  vel_ang = velang ;
  vel_ang.std_spectral_base() ;
  vel_ang.mult_rsint() ;

  
  Scalar bc_source (mp) ;
  bc_source = 0. ;
  bc_source.std_spectral_base() ;
  bc_source.inc_dzpuis(2) ;

  // beta^r component
  //-----------------
  double rho = 10. ; // rho>1 ; rho=1 "pure Dirichlet" version
  Scalar beta_r_corr =  (rho - 1) * b_tilde() * hh_tilde;
  beta_r_corr.inc_dzpuis(2) ;

  Scalar beta_r (mp) ;
  beta_r = 2 * contract(s_tilde, 0, b_tilde().derive_cov(ff), 0) + beta_r_corr - bc_source ;
  beta_r = beta_r / (hh_tilde * rho) ;
    
  beta_r.dec_dzpuis(2) ;
  
  beta_r = beta_r - beta()(2)*s_tilde(2) -  beta()(3)*s_tilde(3) ;
  Vector tmp = s_tilde.down(0, met_gamt) ;
  beta_r = beta_r/tmp(1) ;

  //  cout<<"beta_r"<<beta_r<<endl ;
  //  des_profile(beta_r, 1.00001, 10, M_PI/2., 0., "beta_r") ;
  
  
 
  //  Scalar beta_r (mp) ;
  // beta_r = 0.5 ; 
  beta_r.set_spectral_va().set_base(s_tilde(1).get_spectral_va().get_base()) ;
  

  //Vector tmp_vect = nn() * radial_vect_hor() ;
  //  Vector tmp_vect (mp, CON, mp.get_bvect_spher() ) ;
  //  tmp_vect.set_etat_zero() ; 

  tmp_vect.set(1) = beta_r ;
  //  tmp_vect.set(2) = 0. ;
  tmp_vect.set(3) += - vel_ang ;

  //  tmp_vect.std_spectral_base() ;

  tmp_vect.change_triad(mp.get_bvect_cart() ) ;
  

  /*
  Scalar tmp(mp) ;
  
  tmp = tmp_vect(1) ;
  des_profile(tmp, 1.00001, 10, M_PI/2.,  M_PI/3. , "Source(1)") ;
  
  tmp = tmp_vect(2) ;
  des_profile(tmp, 1.00001, 10, M_PI/2.,  M_PI/3., "Source(2)") ;
  
  tmp = tmp_vect(3) ;
  des_profile(tmp, 1.00001, 10, M_PI/2.,  M_PI/3., "Source(3)") ;
  */  

return tmp_vect ;

}



// Component x of boundary value of V^i 
//-------------------------------------
Valeur Isol_hor:: boundary_vv_x(double velang){
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  //Isol_hor boundary conditions
  
  Valeur lim_x (mp.get_mg()->get_angu()) ;
    
  lim_x = 1 ;   // Why is it necessary this and what it is actually doing?
  
  Scalar vv_x = vv_bound_cart(velang)(1) ;

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_x.set(0, k, j, 0) = vv_x.val_grid_point(1, k, j, 0) ;
  
  lim_x.set_base(vv_x.get_spectral_va().get_base()) ;

  return  lim_x ;


}


// Component y of boundary value of V^i
//--------------------------------------
Valeur Isol_hor:: boundary_vv_y(double velang){
 
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
  
  Valeur lim_y (mp.get_mg()->get_angu()) ;
    
  lim_y = 1 ;   // Why is it necessary this and what it is actually doing?
 
  Scalar vv_y = vv_bound_cart(velang)(2) ;

  for (int k=0 ; k<nnp ; k++)
      for (int j=0 ; j<nnt ; j++)
	  lim_y.set(0, k, j, 0) = vv_y.val_grid_point(1, k, j, 0) ;

  lim_y.set_base(vv_y.get_spectral_va().get_base()) ;

  return  lim_y ;
}


// Component z of boundary value of V^i
//-------------------------------------
Valeur Isol_hor:: boundary_vv_z(double velang){

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
 
  Valeur lim_z (mp.get_mg()->get_angu()) ;
    
  lim_z = 1 ;   // Why is it necessary this and what it is actually doing?
  
  Scalar vv_z = vv_bound_cart(velang)(3) ;
  Scalar vv_x = vv_bound_cart(velang)(1) ; // This will be non zero, and we need a 
                                           // non zero quantity in order to obtain a basis
                                           // for the valeur lim_z

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_z.set(0, k, j, 0) = vv_z.val_grid_point(1, k, j, 0) ;
 
  lim_z.set_base(vv_z.get_spectral_va().get_base()) ;

  return  lim_z ;
}

