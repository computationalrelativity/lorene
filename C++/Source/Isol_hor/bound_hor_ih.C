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

Valeur Isol_hor::boundary_psi_Dir_ih(){

  const Map& map = ff.get_mp() ; 


  Scalar tmp= - 6 * contract(beta(), 0, psi().derive_cov(ff), 0) ;
  tmp / (contract(beta().derive_cov(ff), 0, 1) - nn() * trk() ) - 1 ;

  //Se ha restado 1 porque se resuelve con condicion de cero en el infinito y despues se suma 1 a la solucion   

  Valeur psi_bound (map.get_mg()->get_angu() )  ;
  
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;
			
  psi_bound = 1 ;
			
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  psi_bound.std_base_scal() ;


  /*  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      gamt_rr_bound.set(0, k, j, 0) = gamt.con()(1,1).val_grid_point(1, k, j, 0) ;

  gamt_rr_bound.std_base_scal() ;


  // Boundary condition for Psi^2

  Valeur bc_psi2 (map.get_mg()->get_angu() )  ;

  bc_psi2 = pow( (f_tt * f_pp -  f_tp * f_tp) / gamt_rr_bound, 1./4. )  ;
  */

  
  return psi_bound ;

}


// Neumann boundary condition for Psi 
//-------------------------------------
// ONE HAS TO GUARANTEE THAT BETA IS NOT ZERO, BUT IT IS PROPORTIONAL TO THE RADIAL VECTOR

Valeur Isol_hor::boundary_psi_Neu_ih(){

  const Map& map = ff.get_mp() ; 

  Scalar tmp = - 1./ 6. * psi() * (contract(beta().derive_cov(ff), 0, 1) - nn() * trk() ) 
    - beta()(2)* psi().derive_cov(ff)(2) - beta()(3)* psi().derive_cov(ff)(3) ;

  tmp = tmp / beta()(1) ;

  //En este caso no hay que restar nada a ninguna variable   

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
Valeur Isol_hor::boundary_nn_Dir_kk_ih(){

  const Map& map = ff.get_mp() ;

  Scalar tmp(map) ;

  Scalar kk_rr = contract( radial_vect_hor_ih(), 0, contract(radial_vect_hor_ih(), 0, k_dd(), 0), 0) ;


  //  des_profile(qq, 2.01, 10., M_PI/7., 0., "qq ", "radial distance") ;

  //  kk_rr.set_domain(0) = 1. ; 


  Scalar log_nn = log(nn()) ;
  log_nn.std_spectral_base() ;
  
  Scalar k_kerr (map) ;
  k_kerr = kappa_hor_ih() ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;

  tmp = k_kerr - contract(radial_vect_hor_ih(), 0, nn().derive_cov(ff), 0)
    - omega_hor_ih() * log_nn.derive_cov(ff)(3) ;

  tmp = - tmp / kk_rr - 1;
  
  //  tmp = -1. ;  
  //  tmp.std_spectral_base() ;
  
    
  //Se ha restado 1 porque se resuelve con condicion de cero en el infinito y despues se suma 1 a la solucion   


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
Valeur Isol_hor::boundary_nn_Neu_kk_ih() {
  
  const Map& map = ff.get_mp() ;
  
  Vector dnn= nn().derive_cov(ff) ;

  Scalar kk_rr = contract( radial_vect_hor_ih(), 0, contract(radial_vect_hor_ih(), 0, k_dd(), 0), 0) ;

  Scalar log_nn = log(nn()) ;
  log_nn.std_spectral_base() ;
  
  Scalar k_kerr (map) ;
  k_kerr = kappa_hor_ih() ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;

  Scalar tmp = k_kerr + nn() * kk_rr - omega_hor_ih() * log_nn.derive_cov(ff)(3) 
             - radial_vect_hor_ih()(2) * dnn(2) - radial_vect_hor_ih()(3) * dnn(3)  ;

  tmp = tmp / radial_vect_hor_ih()(1) ;

  //  tmp = 1. ;
  //  tmp.std_spectral_base() ;

  //En este caso no hay que restar nada a ninguna variable   

  tmp.set_etat_zero() ; 


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
Valeur Isol_hor:: boundary_beta_r_ih(){

  const Map& map = ff.get_mp() ; 

  Scalar tmp (map) ;

  Scalar beta_r = nn() * radial_vect_hor_ih()(1) ;
 
  //  tmp = 0.;
  
  tmp = beta_r ;

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
Valeur Isol_hor::boundary_beta_theta_ih(){
  
  const Map& map = ff.get_mp() ; 

  Scalar tmp(map) ;  
  
  Scalar beta_theta = nn() * radial_vect_hor_ih()(2) ;

  tmp = beta_theta ;  
  //  tmp = 0. ;

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
Valeur Isol_hor::boundary_beta_phi_ih(){

  const Map& map = ff.get_mp() ; 

  Scalar tmp (map) ;

  Scalar vel_ang(map) ;

  vel_ang = omega_hor_ih() ;

  vel_ang.std_spectral_base() ;

  vel_ang.mult_rsint() ;

  cout<<"dzpuis: "<< vel_ang.get_dzpuis() <<endl ;

  Scalar beta_phi = nn() * radial_vect_hor_ih()(3)  -  vel_ang ;

  tmp = beta_phi ;
  //  tmp = 0. ;

  //  des_profile(beta_phi, 2., 10., 1.57, 0., "beta_phi", "radial distance") ;


  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Valeur bnd_beta_phi (map.get_mg()->get_angu()) ;
    
  bnd_beta_phi = 1 ;   // Why is it necessary this and what it is actually doing?
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_phi.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  bnd_beta_phi.std_base_scal() ;
  
  return  bnd_beta_phi ;

}




// Component x of boundary value of beta (using expression in terms of radial vector)
//--------------------------------------
Valeur Isol_hor:: boundary_beta_x_ih(){

  const Map& map = ff.get_mp() ; 
  
 
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;


  Mtbl x_mtbl (map.get_mg()) ;
  x_mtbl.set_etat_qcq() ;
  Mtbl y_mtbl (map.get_mg()) ;
  y_mtbl.set_etat_qcq() ;
  x_mtbl = map.x ;
  y_mtbl = map.y ;


  
  // Bases for limit conditions  (WHY IS THIS NECESSARY?)
  Base_val** bases = map.get_mg()->std_base_vect_cart() ;
 
  Valeur lim_x (map.get_mg()->get_angu()) ;
  lim_x = 1 ;
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_x.set(0, k, j, 0) = omega_hor_ih()*y_mtbl(1, k, j, 0) ;
      //      lim_x.set(0, k, j, 0) = y_mtbl(1, k, j, 0) ;

  /*
  //Isol_hor boundary conditions
  
  Valeur lim_x (map.get_mg()->get_angu()) ;
    
  lim_x = 1 ;   // Why is it necessary this and what it is actually doing?
  
  Scalar beta_x = beta_bound_cart()(1) ;


  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_x.set(0, k, j, 0) = beta_x.val_grid_point(1, k, j, 0) ;
  */

  lim_x.base = *bases[0] ;
  
     
  // We do not need it any more
  for (int i=0 ; i<3 ; i++)
    delete bases[i] ;
  delete [] bases ;
    
  return  lim_x ;


}


// Component y of boundary value of beta (using expression in terms of radial vector)
//--------------------------------------
Valeur Isol_hor:: boundary_beta_y_ih(){
  
  const Map& map = ff.get_mp() ; 
  
 
  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;


  Mtbl x_mtbl (map.get_mg()) ;
  x_mtbl.set_etat_qcq() ;
  Mtbl y_mtbl (map.get_mg()) ;
  y_mtbl.set_etat_qcq() ;
  x_mtbl = map.x ;
  y_mtbl = map.y ;


  
  // Bases for limit conditions  (WHY IS THIS NECESSARY?)
  Base_val** bases = map.get_mg()->std_base_vect_cart() ;
      
  Valeur lim_y (map.get_mg()->get_angu()) ;
    lim_y = 1 ;
    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	  lim_y.set(0, k, j, 0) = - omega_hor_ih()*x_mtbl(1, k, j, 0) ;
    //lim_y.set(0, k, j, 0) = - x_mtbl(1, k, j, 0) ;
    
    /*
    // Isol_hor boundary conditions
    
    Valeur lim_y (map.get_mg()->get_angu()) ;
    
    lim_y = 1 ;   // Why is it necessary this and what it is actually doing?
 
    Scalar beta_y = beta_bound_cart()(2) ;
   

    for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
    lim_y.set(0, k, j, 0) = beta_y.val_grid_point(1, k, j, 0) ;
    */

  
    lim_y.base = *bases[1] ;
   


  // We do not need it any more
  for (int i=0 ; i<3 ; i++)
    delete bases[i] ;
  delete [] bases ;

  return  lim_y ;


}


// Component z of boundary value of beta (using expression in terms of radial vector)
//--------------------------------------
Valeur Isol_hor:: boundary_beta_z_ih(){

  const Map& map = ff.get_mp() ; 

  

  int nnp = map.get_mg()->get_np(1) ;
  int nnt = map.get_mg()->get_nt(1) ;

  Base_val** bases = map.get_mg()->std_base_vect_cart() ;
  
  Valeur lim_z (map.get_mg()->get_angu()) ;
  lim_z = 1 ;
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_z.set(0, k, j, 0) = 0. ;
  
    /* 
  // Isol_hor boundary conditions
 
  Valeur lim_z (map.get_mg()->get_angu()) ;
    
  lim_z = 1 ;   // Why is it necessary this and what it is actually doing?
  
  Scalar beta_z = beta_bound_cart()(3) ;


  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_z.set(0, k, j, 0) = beta_z.val_grid_point(1, k, j, 0) ;
    */ 

  lim_z.base = *bases[2] ;
   
  // On n'en a plus besoin
  for (int i=0 ; i<3 ; i++)
    delete bases[i] ;
  delete [] bases ;


 
  return  lim_z ;


}
