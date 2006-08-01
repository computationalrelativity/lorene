/*
 *  Method of class Isol_hor to compute boundary conditions
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

char bound_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.31  2006/08/01 14:35:48  f_limousin
 * Many small modifs
 *
 * Revision 1.30  2006/02/22 17:02:04  f_limousin
 * Removal of warnings
 *
 * Revision 1.29  2006/02/22 16:29:55  jl_jaramillo
 * corrections on the relaxation and boundary conditions
 *
 * Revision 1.27  2005/10/24 16:44:40  jl_jaramillo
 * Cook boundary condition ans minot bound of kss
 *
 * Revision 1.26  2005/10/21 16:20:55  jl_jaramillo
 * Version for the paper JaramL05
 *
 * Revision 1.25  2005/09/13 18:33:17  f_limousin
 * New function vv_bound_cart_bin(double) for computing binaries with
 * berlin condition for the shift vector.
 * Suppress all the symy and asymy in the importations.
 *
 * Revision 1.24  2005/09/12 12:33:54  f_limousin
 * Compilation Warning - Change of convention for the angular velocity
 * Add Berlin boundary condition in the case of binary horizons.
 *
 * Revision 1.23  2005/07/08 13:15:23  f_limousin
 * Improvements of boundary_vv_cart(), boundary_nn_lapl().
 * Add a fonction to compute the departure of axisymmetry.
 *
 * Revision 1.22  2005/06/09 08:05:32  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * Revision 1.21  2005/05/12 14:48:07  f_limousin
 * New boundary condition for the lapse : boundary_nn_lapl().
 *
 * Revision 1.20  2005/04/29 14:04:17  f_limousin
 * Implementation of boundary_vv_x (y,z) for binary black holes.
 *
 * Revision 1.19  2005/04/19 16:40:51  jl_jaramillo
 * change of sign of ang_vel in vv_bound_cart. Convention of Phys. Rep.
 *
 * Revision 1.18  2005/04/08 12:16:52  f_limousin
 * Function set_psi(). And dependance in phi.
 *
 * Revision 1.17  2005/04/02 15:49:21  f_limousin
 * New choice (Lichnerowicz) for aaquad. New member data nz.
 *
 * Revision 1.16  2005/03/22 13:25:36  f_limousin
 * Small changes. The angular velocity and A^{ij} are computed
 * with a differnet sign.
 *
 * Revision 1.15  2005/03/10 10:19:42  f_limousin
 * Add the regularisation of the shift in the case of a single black hole
 * and lapse zero on the horizon.
 *
 * Revision 1.14  2005/03/03 10:00:55  f_limousin
 * The funtions beta_boost_x() and beta_boost_z() have been added.
 *
 * Revision 1.13  2005/02/24 17:21:04  f_limousin
 * Suppression of the function beta_bound_cart() and direct computation
 * of boundary_beta_x, y and z.
 *
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
#include "param.h"


// Dirichlet boundary condition for Psi 
//-------------------------------------
// ONE HAS TO GUARANTEE THAT BETA IS NOT ZERO, BUT IT IS PROPORTIONAL TO THE RADIAL VECTOR

const Valeur Isol_hor::boundary_psi_Dir_evol() const{

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

const Valeur Isol_hor::boundary_psi_Neu_evol()const {

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


const Valeur Isol_hor::boundary_psi_Dir_spat()const {

    Scalar tmp = psi() * psi() * psi() * trk() 
      - contract(k_dd(), 0, 1, tradial_vect_hor() * tradial_vect_hor(), 0, 1) 
      / psi()
      - 4.* contract(tradial_vect_hor(), 0, psi().derive_cov(ff), 0) ;

    // rho = 1 is the standart case.
    double rho = 1. ;
    tmp += (tradial_vect_hor().divergence(ff)) * (rho - 1.)*psi() ;

    tmp = tmp / (rho * tradial_vect_hor().divergence(ff)) - 1. ;

    tmp.std_spectral_base() ;
    tmp.inc_dzpuis(2) ;

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

const Valeur Isol_hor::boundary_psi_Neu_spat()const {

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

const Valeur Isol_hor::boundary_psi_app_hor()const {

    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    Valeur psi_bound (mp.get_mg()->get_angu()) ;

    psi_bound = 1 ; // Juste pour affecter dans espace des configs ;

//    if (psi_comp_evol.is_known(jtime)) {
//    for (int k=0 ; k<nnp ; k++)
//	for (int j=0 ; j<nnt ; j++)
//	    psi_bound.set(0, k, j, 0) = - 0.5/radius*(psi_auto().val_grid_point(1, k, j, 0) + psi_comp().val_grid_point(1, k, j, 0)) ;
//    }
//    else {
    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	    psi_bound.set(0, k, j, 0) = - 0.5/radius*psi().val_grid_point(1, k, j, 0) ;
//    }


    psi_bound.std_base_scal() ;

    return psi_bound ;

}

const Valeur Isol_hor::boundary_psi_Dir()const {

  Scalar rho (mp) ;
  rho = 5. ;
  rho.std_spectral_base() ;


  Scalar tmp(mp) ;
  //  tmp = nn()+1. - 1 ;

  
  //Scalar b_tilde_tmp = b_tilde() ;
  //b_tilde_tmp.set_domain(0) = 1. ;
  //tmp = pow(nn()/b_tilde_tmp, 0.5) ;
  

  tmp = 1.5 ;
  tmp.std_spectral_base() ;

  //tmp = 1/ (2* nn()) ;

  tmp = (tmp + rho * psi())/(1 + rho) ; 

  tmp = tmp - 1. ;



  Valeur psi_bound (mp.get_mg()->get_angu()) ;
    
  psi_bound = 1 ;  // Juste pour affecter dans espace des configs ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  psi_bound.std_base_scal() ;
  
  return  psi_bound ;

}

// Dirichlet boundary condition on nn using the extrinsic curvature
// (No time evolution taken into account! Make this)
//--------------------------------------------------------------------------
const Valeur Isol_hor::boundary_nn_Dir_kk()const {

  Scalar tmp(mp) ;
  double rho = 0. ;

  //  Scalar kk_rr = 0.8*psi().derive_cov(mp.flat_met_spher())(1) / psi() ;
  Scalar kk_rr = contract( gam().radial_vect() * gam().radial_vect(), 0, 1
			   , k_dd(), 0, 1 ) ;

  Scalar k_kerr (mp) ;
  k_kerr = 0.1 ;//1.*kappa_hor() ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;

  Scalar temp (rho*nn()) ;
  temp.inc_dzpuis(2) ;

  tmp = contract(gam().radial_vect(), 0, nn().derive_cov(ff), 0) + temp
      - k_kerr ;

  tmp = tmp / (kk_rr + rho) - 1;
  
  Scalar diN (contract(gam().radial_vect(), 0, nn().derive_cov(ff), 0)) ;
  cout << "k_kerr = " << k_kerr.val_grid_point(1, 0, 0, 0) << endl ;
  cout << "D^i N  = " << diN.val_grid_point(1, 0, 0, 0) << endl ;
  cout << "kss = " << kk_rr.val_grid_point(1, 0, 0, 0) << endl ;

  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

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


const Valeur Isol_hor::boundary_nn_Neu_kk(int step) const {
  
  const Vector& dnnn = nn().derive_cov(ff) ;
  double rho = 5. ;
  
  step = 100 ;  // Truc bidon pour eviter warning

  Scalar kk_rr = contract( gam().radial_vect() * gam().radial_vect(), 0, 1
			   , k_dd(), 0, 1 ) ; 
 
  Scalar k_kerr (mp) ;
  k_kerr = kappa_hor() ; 
  //k_kerr = 0.06 ;
  k_kerr.std_spectral_base() ;
  k_kerr.inc_dzpuis(2) ;
  
  Scalar k_hor (mp) ;
  k_hor =  kappa_hor() ; 
  k_hor.std_spectral_base() ;

  //  Vector beta_tilde_d (beta().down(0, met_gamt)) ;
  //Scalar naass = 1./2. * contract( met_gamt.radial_vect() * met_gamt.radial_vect(), 0, 1, beta_tilde_d.ope_killing_conf(met_gamt) , 0, 1 ) ;
  
  Scalar naass = contract( met_gamt.radial_vect() * met_gamt.radial_vect(), 
			   0, 1, aa_auto().up_down(met_gamt), 0, 1) ;

  Scalar traceK = 1./3. * nn() * trk() * 
               contract( met_gamt.radial_vect() * met_gamt.radial_vect(), 0, 1
			   , met_gamt.cov() , 0, 1 ) ;

  Scalar sdN (contract(gam().radial_vect(), 0, nn().derive_cov(ff), 0)) ;



  Scalar tmp =  psi() * psi() * ( k_kerr +  naass + 1.* traceK) 
    -  met_gamt.radial_vect()(2) * dnnn(2) 
    - met_gamt.radial_vect()(3) * dnnn(3) 
    + ( rho - 1 ) * met_gamt.radial_vect()(1)  * dnnn(1)  ;


  tmp = tmp / (rho * met_gamt.radial_vect()(1))  ;




  Scalar rhs ( sdN - nn() * kk_rr) ;
  cout << "kappa_pres = " << k_kerr.val_grid_point(1, 0, 0, 0) << endl ;
  cout << "kappa_hor = " << k_hor.val_grid_point(1, 0, 0, 0) << endl ;
  cout << "sDN  = " << sdN.val_grid_point(1, 0, 0, 0) << endl ;

  // in this case you don't have to substract any value
 
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

const Valeur Isol_hor::boundary_nn_Neu_Cook() const {
  
  const Vector& dnnn = nn().derive_cov(ff) ;
  double rho = 1. ;


  Scalar kk_rr = contract( gam().radial_vect() * gam().radial_vect(), 0, 1
			   , k_dd(), 0, 1 ) ; 

  Sym_tensor qq_uu = gam_uu() - gam().radial_vect() * gam().radial_vect() ;

  // Free function L_theta = h
  //--------------------------
  Scalar L_theta (mp) ;
  L_theta = .0;
  L_theta.std_spectral_base() ;
  L_theta.inc_dzpuis(4) ;

  //Divergence of Omega
  //-------------------
  Vector hom1 = nn().derive_cov(met_gamt)/nn()  ;
  hom1 = contract(qq_uu.down(1, gam()), 0, hom1, 0) ;

  Vector hom2 = -contract( k_dd(), 1, gam().radial_vect() , 0) ;
  hom2 = contract(qq_uu.down(1, gam()), 0, hom2, 0) ;

  Vector hom = hom1 + hom2 ;
  
  Scalar div_omega = 1.*contract( qq_uu, 0, 1,  (1.*hom1 + 1.*hom2).derive_cov(gam()), 0, 1) ;
  div_omega.inc_dzpuis() ;

  //---------------------

  
  // Two-dimensional Ricci scalar
  //-----------------------------

  Scalar rr (mp) ;
  rr = mp.r ;

  Scalar Ricci_conf(mp) ;
  Ricci_conf = 2 / rr / rr ;
  Ricci_conf.std_spectral_base() ;

  Scalar Ricci(mp) ;
  Scalar log_psi (log(psi())) ;
  log_psi.std_spectral_base() ;
  Ricci = pow(psi(), -4) * (Ricci_conf - 4*log_psi.lapang())/rr /rr ;
  Ricci.std_spectral_base() ;
  Ricci.inc_dzpuis(4) ;
  //-------------------------------

 Scalar theta_k = -1/(2*nn()) * (gam().radial_vect().divergence(gam()) -
				  kk_rr + trk() ) ;
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;
  /*  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++){
       cout << theta_k.val_grid_point(1, k, j, 0) << endl ;
       cout << nn().val_grid_point(1, k, j, 0) << endl ;
    }
  */



  //Source
  //------
  Scalar  source = div_omega + contract( qq_uu, 0, 1,  hom * hom , 0, 1) - 0.5 * Ricci - L_theta;
  source = source / theta_k ;

  Scalar tmp = ( source + nn() * kk_rr + rho * contract(gam().radial_vect(), 0,
							nn().derive_cov(gam()), 0))/(1+rho) 
    - gam().radial_vect()(2) * dnnn(2) - gam().radial_vect()(3) * dnnn(3)  ;
  
  tmp = tmp / gam().radial_vect()(1) ;

  // in this case you don't have to substract any value
 
  Valeur nn_bound (mp.get_mg()->get_angu()) ;
  
  nn_bound = 1 ;  // Juste pour affecter dans espace des configs ; 
  
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  nn_bound.std_base_scal() ;
  
  return  nn_bound ;
  
}



const Valeur Isol_hor::boundary_nn_Dir_eff(double cc)const {

  Scalar tmp(mp) ;

  tmp = - nn().derive_cov(ff)(1) ;

  // rho = 1 is the standart case
  double rho = 1. ;
  tmp.dec_dzpuis(2) ;
  tmp += cc * (rho -1)*nn() ;
  tmp = tmp / (rho*cc) ;

  tmp = tmp - 1. ;
  
  // We have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

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



const Valeur Isol_hor::boundary_nn_Neu_eff(double cc)const  {
  
  Scalar tmp = - cc * nn() ;
  //  Scalar tmp = - nn()/psi()*psi().dsdr() ;

  // in this case you don't have to substract any value
 
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


const Valeur Isol_hor::boundary_nn_Dir(double cc)const {

  Scalar rho(mp);
  rho = 0. ;  // 0 is the standard case

  Scalar tmp(mp) ;
  tmp = cc ; 
  
  /*
  if (b_tilde().val_grid_point(1, 0, 0, 0) < 0.08)
    tmp = 0.25   ;
  else {
    cout << "OK" << endl ;
    //   des_profile(nn(), 0, 20, M_PI/2, M_PI) ;
    rho = 5. ;
    tmp =   b_tilde()*psi()*psi() ;
  }
  */
  
  //tmp = 1./(2*psi()) ;
  //  tmp = - psi() * nn().dsdr() / (psi().dsdr())  ;

  // We  have substracted 1, since we solve for zero condition at infinity 
  //and then we add 1 to the solution  

  tmp = (tmp + rho * nn())/(1 + rho) ;

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

const Valeur Isol_hor::boundary_nn_Dir_lapl(int mer) const {

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;
   

  //Preliminaries
  //-------------
  
  //Metrics on S^2
  //--------------

  Sym_tensor q_uu = gam_uu() - gam().radial_vect() * gam().radial_vect() ;
  Sym_tensor q_dd = q_uu.up_down(gam()) ;

  Scalar det_q = q_dd(2,2) * q_dd(3,3) - q_dd(2,3) * q_dd(3,2) ;
  Scalar square_q (pow(det_q, 0.5)) ;
  square_q.std_spectral_base() ;

  Sym_tensor qhat_uu = square_q * q_uu ;

  Sym_tensor ffr_uu = ff.con() - ff.radial_vect() * ff.radial_vect() ;
  Sym_tensor ffr_dd = ff.cov() - ff.radial_vect().down(0, ff) * ff.radial_vect().down(0,ff) ;

  Sym_tensor h_uu = qhat_uu - ffr_uu ;

  //------------------


  // Source of the "round" laplacian: 
  //---------------------------------

  //Source 1: "Divergence" term
  //---------------------------
  Vector kqs = contract(k_dd(), 1, gam().radial_vect(), 0 ) ;
  kqs = contract( q_uu.down(0, gam()), 1, kqs, 0) ;
  Scalar div_kqs = contract( q_uu, 0, 1, kqs.derive_cov(gam()), 0, 1) ;

  Scalar source1 = div_kqs ;
  source1 *= square_q ;


  // Source 2: correction term
  //--------------------------
    
  Vector corr = - contract(h_uu, 1, nn().derive_cov(ff), 0) / nn() ;
  Scalar source2 = contract( ffr_dd, 0, 1, corr.derive_con(ff), 0, 1 ) ;


  // Source 3: (physical) divergence of Omega
  //-----------------------------------------

  Scalar div_omega(mp) ;
  
  //Source from (L_l\theta_k=0)
  
  Scalar L_theta (mp) ;
  L_theta = 0. ;
  L_theta.std_spectral_base() ;
  
  Scalar kk_rr = contract( gam().radial_vect() * gam().radial_vect(), 0, 1
			   , k_dd(), 0, 1 ) ; 
     
 //   Scalar kappa (mp) ;
 //   kappa = kappa_hor() ;
 //   kappa.std_spectral_base() ;
 //   kappa.set_dzpuis(2) ;
  
 
  Scalar kappa = contract(gam().radial_vect(), 0, nn().derive_cov(gam()), 0) - nn() * kk_rr ;
  Scalar theta_k = -1/(2*nn()) * (gam().radial_vect().divergence(gam()) -
				  kk_rr  + trk() ) ;
  
  Sym_tensor qqq = gam_uu() - gam().radial_vect() * gam().radial_vect() ;
  
  Vector hom = nn().derive_cov(met_gamt) - contract( k_dd(), 1, gam().radial_vect() , 0) ;
  hom = contract(qqq.down(1, gam()), 0, hom, 0) ;
  
  Scalar rr(mp) ;
  rr = mp.r ;
  
  Scalar Ricci_conf = 2 / rr /rr ;
  Ricci_conf.std_spectral_base() ;
  
  Scalar log_psi (log(psi())) ;
  log_psi.std_spectral_base() ;
  Scalar Ricci = pow(psi(), -4) * (Ricci_conf - 4*log_psi.lapang())/rr /rr ;
  Ricci.std_spectral_base() ;
  Ricci.inc_dzpuis(4) ;
  
  
  div_omega =  L_theta - contract(qqq, 0, 1, hom * hom, 0, 1) + 0.5 * Ricci
    + theta_k * kappa ;
  div_omega.dec_dzpuis() ;
    
  //End of Source from (L_l\theta_k=0)
  //----------------------------------

  div_omega = 0. ;
  //  div_omega = 0.01*log(1/(2*psi())) ;
  div_omega.std_spectral_base() ;
  div_omega.inc_dzpuis(3) ;

 
  //Construction of source3 (correction of div_omega by the square root of the determinant)
  Scalar source3 =  div_omega ;
  source3 *= square_q ;


  //"Correction" to the div_omega term (no Y_00 component must be present)
  
  double corr_const = mp.integrale_surface(source3, radius  + 1e-15)/(4. * M_PI) ;
  cout << "Constant part of div_omega = " << corr_const <<endl ;

  Scalar corr_div_omega (mp) ;
  corr_div_omega = corr_const ;
  corr_div_omega.std_spectral_base() ;
  corr_div_omega.set_dzpuis(3) ;

  source3 -=  corr_div_omega ;

 

  //Source
  //------
  
  Scalar source =  source1 + source2 + source3 ;



  //Resolution of round angular laplacian
  //-------------------------------------
  
  Scalar source_tmp = source ;

  Scalar logn (mp) ;
  logn = source.poisson_angu() ;

  double cc = 0.2 ; // Integration constant

  Scalar lapse (mp) ;
  lapse = (exp(logn)*cc) ;
  lapse.std_spectral_base() ;
  

  //Test of Divergence of Omega
  //---------------------------
  
  ofstream err ;
  err.open ("err_laplBC.d", ofstream::app ) ;
  
  hom = nn().derive_cov(gam())/nn() 
    - contract( k_dd(), 1, gam().radial_vect() , 0) ;
  hom = contract(q_uu.down(1, gam()), 0, hom, 0) ;
  
  Scalar div_omega_test = contract( q_uu, 0, 1,  hom.derive_cov(gam()), 0, 1) ;
  
  
  Scalar err_lapl = div_omega_test - div_omega ;
  
  double max_err = err_lapl.val_grid_point(1, 0, 0, 0) ;
  double min_err = err_lapl.val_grid_point(1, 0, 0, 0) ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++){
      if (err_lapl.val_grid_point(1, k, j, 0) > max_err)
	max_err = err_lapl.val_grid_point(1, k, j, 0) ;
      if (err_lapl.val_grid_point(1, k, j, 0) < min_err)
	min_err = err_lapl.val_grid_point(1, k, j, 0) ;
    }
  
   err << mer << " " << max_err << " " << min_err << endl ;

   err.close() ;



  //Construction of the Valeur
  //--------------------------

  lapse = lapse - 1. ;

  //  int nnp = mp.get_mg()->get_np(1) ;
  //  int nnt = mp.get_mg()->get_nt(1) ;
  
  Valeur nn_bound (mp.get_mg()->get_angu()) ;
  
  nn_bound = 1 ;  // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = lapse.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
    return  nn_bound ;

}


// Component r of boundary value of beta (using expression in terms 
// of radial vector)
// --------------------------------------

const Valeur Isol_hor:: boundary_beta_r()const {

  Scalar tmp (mp) ;

  tmp = nn() * radial_vect_hor()(1) ;

  Base_val base_tmp (radial_vect_hor()(1).get_spectral_va().base) ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_r (mp.get_mg()->get_angu()) ;
    
  bnd_beta_r = 1. ;  // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_r.set(1, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;

  bnd_beta_r.set_base_r(0, base_tmp.get_base_r(0)) ;
  for (int l=0 ; l<(*mp.get_mg()).get_nzone()-1 ; l++)
      bnd_beta_r.set_base_r(l, base_tmp.get_base_r(l)) ;
  bnd_beta_r.set_base_r((*mp.get_mg()).get_nzone()-1, base_tmp.get_base_r((*mp.get_mg()).get_nzone()-1)) ;
  bnd_beta_r.set_base_t(tmp.get_spectral_va().get_base().get_base_t(1)) ;
  bnd_beta_r.set_base_p(tmp.get_spectral_va().get_base().get_base_p(1)) ;

//  bnd_beta_r.set_base(*(mp.get_mg()->std_base_vect_spher()[0])) ;

  cout << "norme de lim_vr" << endl << norme(bnd_beta_r) << endl ;
  cout << "bases" << endl << bnd_beta_r.base << endl ;
  
  return  bnd_beta_r ;


}


// Component theta of boundary value of beta (using expression in terms 
// of radial vector)
// ------------------------------------------

const Valeur Isol_hor::boundary_beta_theta()const {
  
  Scalar tmp(mp) ;  
  
  tmp = nn() * radial_vect_hor()(2) ;
  Base_val base_tmp (radial_vect_hor()(2).get_spectral_va().base) ;


  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_theta (mp.get_mg()->get_angu()) ;
    
  bnd_beta_theta = 1. ;   // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_theta.set(1, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
//  bnd_beta_theta.set_base(*(mp.get_mg()->std_base_vect_spher()[1])) ;

  bnd_beta_theta.set_base_r(0, base_tmp.get_base_r(0)) ;
  for (int l=0 ; l<(*mp.get_mg()).get_nzone()-1 ; l++)
      bnd_beta_theta.set_base_r(l, base_tmp.get_base_r(l)) ;
  bnd_beta_theta.set_base_r((*mp.get_mg()).get_nzone()-1, base_tmp.get_base_r((*mp.get_mg()).get_nzone()-1)) ;
  bnd_beta_theta.set_base_t(base_tmp.get_base_t(1)) ;
  bnd_beta_theta.set_base_p(base_tmp.get_base_p(1)) ;

  cout << "bases" << endl << bnd_beta_theta.base << endl ;


  return  bnd_beta_theta ;


}
 
// Component phi of boundary value of beta (using expression in terms 
// of radial vector) 
// -----------------------------------------------------------

const Valeur Isol_hor::boundary_beta_phi(double om)const {

  Scalar tmp (mp) ;

  Scalar ang_vel(mp) ;
  ang_vel = om ;
  ang_vel.std_spectral_base() ;
  ang_vel.mult_rsint() ;

  tmp = nn() * radial_vect_hor()(3)  -  ang_vel ;
  Base_val base_tmp (ang_vel.get_spectral_va().base) ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_phi (mp.get_mg()->get_angu()) ;
    
  bnd_beta_phi = 1. ; // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_phi.set(1, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
    
  for (int l=0 ; l<(*mp.get_mg()).get_nzone()-1 ; l++)

//  bnd_beta_phi.set_base(*(mp.get_mg()->std_base_vect_spher()[2])) ;

  bnd_beta_phi.set_base_r(0, base_tmp.get_base_r(0)) ;
  for (int l=0 ; l<(*mp.get_mg()).get_nzone()-1 ; l++)
      bnd_beta_phi.set_base_r(l, base_tmp.get_base_r(l)) ;
  bnd_beta_phi.set_base_r((*mp.get_mg()).get_nzone()-1, base_tmp.get_base_r((*mp.get_mg()).get_nzone()-1)) ;
  bnd_beta_phi.set_base_t(base_tmp.get_base_t(1)) ;
  bnd_beta_phi.set_base_p(base_tmp.get_base_p(1)) ;

//  cout << "bound beta_phi" << endl << bnd_beta_phi << endl ;

  return  bnd_beta_phi ;

}

// Component x of boundary value of beta
//--------------------------------------

const Valeur Isol_hor:: boundary_beta_x(double om)const {

    // Les alignemenents pour le signe des CL.
    double orientation = mp.get_rot_phi() ;
    assert ((orientation == 0) || (orientation == M_PI)) ;
    int aligne = (orientation == 0) ? 1 : -1 ;
   
    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    Vector tmp_vect = nn() * gam().radial_vect() ;
    tmp_vect.change_triad(mp.get_bvect_cart() ) ;

    //Isol_hor boundary conditions
  
    Valeur lim_x (mp.get_mg()->get_angu()) ;
    
    lim_x = 1 ;  // Juste pour affecter dans espace des configs ;
  
    Mtbl ya_mtbl (mp.get_mg()) ;
    ya_mtbl.set_etat_qcq() ;
    ya_mtbl = mp.ya ;

    Scalar cosp (mp) ;
    cosp = mp.cosp ;
    Scalar cost (mp) ;
    cost = mp.cost ;
    Scalar sinp (mp) ;
    sinp = mp.sinp ;
    Scalar sint (mp) ;
    sint = mp.sint ;

    Scalar dep_phi (mp) ;
    dep_phi = 0.0*sint*cosp ;

    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	  lim_x.set(0, k, j, 0) = aligne * om * ya_mtbl(1, k, j, 0) * (1 + 
				  dep_phi.val_grid_point(1, k, j, 0))
	    + tmp_vect(1).val_grid_point(1, k, j, 0) ;
    
  lim_x.set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;

  return  lim_x ;
}


// Component y of boundary value of beta 
//--------------------------------------

const Valeur Isol_hor:: boundary_beta_y(double om)const {

    // Les alignemenents pour le signe des CL.
    double orientation = mp.get_rot_phi() ;
    assert ((orientation == 0) || (orientation == M_PI)) ;
    int aligne = (orientation == 0) ? 1 : -1 ;
    

    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    Vector tmp_vect = nn() * gam().radial_vect() ;
    tmp_vect.change_triad(mp.get_bvect_cart() ) ;

    //Isol_hor boundary conditions
  
    Valeur lim_y (mp.get_mg()->get_angu()) ;
    
    lim_y = 1 ;  // Juste pour affecter dans espace des configs ;
  
    Mtbl xa_mtbl (mp.get_mg()) ;
    xa_mtbl.set_etat_qcq() ;
    xa_mtbl = mp.xa ;

    Scalar cosp (mp) ;
    cosp = mp.cosp ;
    Scalar cost (mp) ;
    cost = mp.cost ;
    Scalar sinp (mp) ;
    sinp = mp.sinp ;
    Scalar sint (mp) ;
    sint = mp.sint ;

    Scalar dep_phi (mp) ;
    dep_phi = 0.0*sint*cosp ;
    
    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	  lim_y.set(0, k, j, 0) = - aligne * om * xa_mtbl(1, k, j, 0) *(1 +
					dep_phi.val_grid_point(1, k, j, 0))
	    + tmp_vect(2).val_grid_point(1, k, j, 0) ;
    
  lim_y.set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;

  return  lim_y ;
}

// Component z of boundary value of beta 
//--------------------------------------

const Valeur Isol_hor:: boundary_beta_z()const {
    
    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    Vector tmp_vect = nn() * gam().radial_vect() ;
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

const Valeur Isol_hor::beta_boost_x() const {
 
    Valeur lim_x (mp.get_mg()->get_angu()) ;
    lim_x = boost_x ;  // Juste pour affecter dans espace des configs ;
  
    lim_x.set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;

    return  lim_x ;
}
   

const Valeur Isol_hor::beta_boost_z() const {
 
    Valeur lim_z (mp.get_mg()->get_angu()) ;
    lim_z = boost_z ;  // Juste pour affecter dans espace des configs ;
      
    lim_z.set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;

  return  lim_z ;
}
      

// Neumann boundary condition for b_tilde 
//---------------------------------------

const Valeur Isol_hor::boundary_b_tilde_Neu()const {
  
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


const Valeur Isol_hor::boundary_b_tilde_Dir()const {

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

const Vector Isol_hor::vv_bound_cart(double om) const{

    // Preliminaries
    //--------------

    // HH_tilde
    Vector s_tilde =  met_gamt.radial_vect() ;

    Scalar hh_tilde = contract(s_tilde.derive_cov(met_gamt), 0, 1) ;
    hh_tilde.dec_dzpuis(2) ;

    // Tangential part of the shift
    Vector tmp_vect = b_tilde() * s_tilde ;

    Vector VV =  b_tilde() * s_tilde - beta() ;

    // "Acceleration" term V^a \tilde{D}_a ln M
    Scalar accel = 2 * contract( VV, 0, contract( s_tilde, 0, s_tilde.down(0,
									   met_gamt).derive_cov(met_gamt), 1), 0) ;


    // Divergence of V^a
    Sym_tensor qq_spher = met_gamt.con() - s_tilde * s_tilde ;
    Scalar div_VV = contract( qq_spher.down(0, met_gamt), 0, 1, VV.derive_cov(met_gamt), 0, 1) ;


    Scalar cosp (mp) ;
    cosp = mp.cosp ;
    Scalar cost (mp) ;
    cost = mp.cost ;
    Scalar sinp (mp) ;
    sinp = mp.sinp ;
    Scalar sint (mp) ;
    sint = mp.sint ;

    Scalar dep_phi (mp) ;
    dep_phi = 0.0*sint*cosp ;

    // Les alignemenents pour le signe des CL.
    double orientation = mp.get_rot_phi() ;
    assert ((orientation == 0) || (orientation == M_PI)) ;
    int aligne = (orientation == 0) ? 1 : -1 ;

    Vector angular (mp, CON, mp.get_bvect_cart()) ;
    Scalar yya (mp) ;
    yya = mp.ya ;
    Scalar xxa (mp) ;
    xxa = mp.xa ;

    angular.set(1) = - aligne * om * yya * (1 + dep_phi) ;
    angular.set(2) = aligne * om * xxa * (1 + dep_phi) ;
    angular.set(3).annule_hard() ;


    angular.set(1).set_spectral_va()
	.set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;
    angular.set(2).set_spectral_va()
	.set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;
    angular.set(3).set_spectral_va()
	.set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;

    angular.change_triad(mp.get_bvect_spher()) ;

/*
  Scalar ang_vel (mp) ;
  ang_vel = om * (1 + dep_phi) ;
  ang_vel.std_spectral_base() ;
  ang_vel.mult_rsint() ;
*/

    Scalar kss (mp) ;
    kss = - 0.2 ;
    //    kss.std_spectral_base() ;
    //    kss.inc_dzpuis(2) ;
  
    
    //kss from from L_lN=0 condition
    //------------------------------
    kss = - 0.15 / nn() ;
    kss.inc_dzpuis(2) ;
    kss += contract(gam().radial_vect(), 0, nn().derive_cov(ff), 0) / nn() ;
    

    cout << "kappa_hor = " << kappa_hor() << endl ;

    /*
    // Apparent horizon condition
    //---------------------------
    kss = trk() ;
    kss -= (4* contract(psi().derive_cov(met_gamt), 0, met_gamt.radial_vect(), 
			0)/psi() +
	    contract(met_gamt.radial_vect().derive_cov(met_gamt), 0, 1)) /
      (psi() * psi()) ;
    */
  

    // beta^r component
    //-----------------
    double rho = 5. ; // rho>1 ; rho=1 "pure Dirichlet" version
    Scalar beta_r_corr = (rho - 1) * b_tilde() * hh_tilde;
    beta_r_corr.inc_dzpuis(2) ;
    Scalar nn_KK = nn() * trk() ;
    nn_KK.inc_dzpuis(2) ;
    
    beta_r_corr.set_dzpuis(2) ;
    nn_KK.set_dzpuis(2) ;
    accel.set_dzpuis(2) ;
    div_VV.set_dzpuis(2) ;
    
    Scalar b_tilde_new (mp) ;
    b_tilde_new = 2 * contract(s_tilde, 0, b_tilde().derive_cov(ff), 0)
	+ beta_r_corr
	- 3 * nn() * kss + nn_KK + accel + div_VV  ;
    
    b_tilde_new = b_tilde_new / (hh_tilde * rho) ;
    
    b_tilde_new.dec_dzpuis(2) ;
    
    tmp_vect.set(1) = b_tilde_new * s_tilde(1) + angular(1) ;
    tmp_vect.set(2) = b_tilde_new * s_tilde(2) + angular(2) ;
    tmp_vect.set(3) = b_tilde_new * s_tilde(3) + angular(3) ;
    
    tmp_vect.change_triad(mp.get_bvect_cart() ) ;
    
    return tmp_vect ;    
}


const Vector Isol_hor::vv_bound_cart_bin(double om, int jj) const{

  // Les alignemenents pour le signe des CL.
  double orientation = mp.get_rot_phi() ;
  assert ((orientation == 0) || (orientation == M_PI)) ;
  int aligne = (orientation == 0) ? 1 : -1 ;
  
  Vector angular (mp, CON, mp.get_bvect_cart()) ;
  Scalar yya (mp) ;
  yya = mp.ya ;
  Scalar xxa (mp) ;
  xxa = mp.xa ;
  
  angular.set(1) = aligne * om * yya ;
  angular.set(2) = - aligne * om * xxa ;
  angular.set(3).annule_hard() ;
  
  angular.set(1).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;
  angular.set(2).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;
  angular.set(3).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;
  
  angular.change_triad(mp.get_bvect_spher()) ;
  
  // HH_tilde
  Vector s_tilde =  met_gamt.radial_vect() ;
  
  Scalar hh_tilde = contract(s_tilde.derive_cov(met_gamt), 0, 1) ;
  hh_tilde.dec_dzpuis(2) ;
  
  Scalar btilde = b_tilde() - contract(angular, 0, s_tilde.up_down(met_gamt), 0) ;
  
  // Tangential part of the shift
  Vector tmp_vect = btilde * s_tilde ;
  
  // Value of kss
  // --------------
  
  Scalar kss (mp) ;
    kss = - 0.26 ;
  kss.std_spectral_base() ;
  kss.inc_dzpuis(2) ;
  
  // Apparent horizon boundary condition
  // -----------------------------------

  
  // kss frome a fich

  // Construction of the binary
  // --------------------------
  
  int depth_s = 3 ;
  char* name_fich = "/home/francois/resu/bin_hor/Test/b11_9x9x8/bin.dat" ;
  
  FILE* fich_s = fopen(name_fich, "r") ;
  Mg3d grid_s (fich_s) ;
  Map_af map_un_s (grid_s, fich_s) ;
  Map_af map_deux_s (grid_s, fich_s) ;
  Bin_hor bin (map_un_s, map_deux_s, fich_s, true, depth_s) ;
  fclose(fich_s) ;
  
  // Inititialisation of fields :
  // ---------------------------- 
  
  bin.set(1).n_comp (bin(2)) ;
  bin.set(1).psi_comp (bin(2)) ;
  bin.set(2).n_comp (bin(1)) ;
  bin.set(2).psi_comp (bin(1)) ;
  bin.decouple() ;
  bin.extrinsic_curvature() ;
  
  kss = contract(bin(jj).k_dd(), 0, 1, bin(jj).gam().radial_vect()*
		 bin(jj).gam().radial_vect(), 0, 1) ;


  /*
  kss = (contract(gam().radial_vect(), 0, nn().derive_cov(ff), 0)
	 - 0.1) / nn() ;
  kss.inc_dzpuis(2) ;
  */

  cout << "kappa_hor = " << kappa_hor() << endl ;
  
  // beta^r component
  //-----------------
  double rho = 5. ; // rho=0 "pure Dirichlet" version
  Scalar beta_r_corr = rho * btilde * hh_tilde;
  beta_r_corr.inc_dzpuis(2) ;
  Scalar nn_KK = nn() * trk() ;
  nn_KK.inc_dzpuis(2) ;
  
  beta_r_corr.set_dzpuis(2) ;
  nn_KK.set_dzpuis(2) ;
  
  Scalar b_tilde_new (mp) ;
  b_tilde_new = 2 * contract(s_tilde, 0, btilde.derive_cov(ff), 0)
    + beta_r_corr - 3 * nn() * kss + nn_KK ;
  
  b_tilde_new = b_tilde_new / (hh_tilde * (rho+1)) ;
  
  b_tilde_new.dec_dzpuis(2) ;
  
  tmp_vect.set(1) = b_tilde_new * s_tilde(1) + angular(1) ;
  tmp_vect.set(2) = b_tilde_new * s_tilde(2) + angular(2) ;
  tmp_vect.set(3) = b_tilde_new * s_tilde(3) + angular(3) ; 
  
  tmp_vect.change_triad(mp.get_bvect_cart() ) ;
  
  return tmp_vect ;
  
}


// Component x of boundary value of V^i 
//-------------------------------------

const Valeur Isol_hor:: boundary_vv_x(double om)const {
  
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  //Isol_hor boundary conditions
  
  Valeur lim_x (mp.get_mg()->get_angu()) ;
    
  lim_x = 1 ;  // Juste pour affecter dans espace des configs ;
 
  Scalar vv_x = vv_bound_cart(om)(1) ;

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_x.set(0, k, j, 0) = vv_x.val_grid_point(1, k, j, 0) ;
  
  lim_x.set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;

  return  lim_x ;


}


// Component y of boundary value of V^i
//--------------------------------------

const Valeur Isol_hor:: boundary_vv_y(double om)const {
 
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
  
  Valeur lim_y (mp.get_mg()->get_angu()) ;
    
  lim_y = 1 ;  // Juste pour affecter dans espace des configs ;
 
  Scalar vv_y = vv_bound_cart(om)(2) ;

  for (int k=0 ; k<nnp ; k++)
      for (int j=0 ; j<nnt ; j++)
	  lim_y.set(0, k, j, 0) = vv_y.val_grid_point(1, k, j, 0) ;

  lim_y.set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;

  return  lim_y ;
}


// Component z of boundary value of V^i
//-------------------------------------

const Valeur Isol_hor:: boundary_vv_z(double om)const {

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  // Isol_hor boundary conditions
 
  Valeur lim_z (mp.get_mg()->get_angu()) ;
    
  lim_z = 1 ;   // Juste pour affecter dans espace des configs ;
  
  Scalar vv_z = vv_bound_cart(om)(3) ;

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      lim_z.set(0, k, j, 0) = vv_z.val_grid_point(1, k, j, 0) ;
 
  lim_z.set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;

  return  lim_z ;
}

// Component x of boundary value of V^i
//-------------------------------------

const Valeur Isol_hor:: boundary_vv_x_bin(double om, int jj)const {

    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    //Isol_hor boundary conditions

    Valeur lim_x (mp.get_mg()->get_angu()) ;

    lim_x = 1 ;  // Juste pour affecter dans espace des configs ;

    Scalar vv_x = vv_bound_cart_bin(om, jj)(1) ;

    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	    lim_x.set(0, k, j, 0) = vv_x.val_grid_point(1, k, j, 0) ;

    lim_x.set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;

    return  lim_x ;


}

// Component y of boundary value of V^i
//--------------------------------------

const Valeur Isol_hor:: boundary_vv_y_bin(double om, int jj)const {

    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    // Isol_hor boundary conditions

    Valeur lim_y (mp.get_mg()->get_angu()) ;

    lim_y = 1 ;  // Juste pour affecter dans espace des configs ;

    Scalar vv_y = vv_bound_cart_bin(om, jj)(2) ;

    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	    lim_y.set(0, k, j, 0) = vv_y.val_grid_point(1, k, j, 0) ;

    lim_y.set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;

    return  lim_y ;
}


// Component z of boundary value of V^i
//-------------------------------------

const Valeur Isol_hor:: boundary_vv_z_bin(double om, int jj)const {

    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;

    // Isol_hor boundary conditions

    Valeur lim_z (mp.get_mg()->get_angu()) ;

    lim_z = 1 ;   // Juste pour affecter dans espace des configs ;

    Scalar vv_z = vv_bound_cart_bin(om, jj)(3) ;

    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++)
	    lim_z.set(0, k, j, 0) = vv_z.val_grid_point(1, k, j, 0) ;

    lim_z.set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;

    return  lim_z ;
}
