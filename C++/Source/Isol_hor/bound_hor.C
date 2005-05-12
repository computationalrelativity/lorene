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


// Dirichlet boundary condition on nn using the extrinsic curvature
// (No time evolution taken into account! Make this)
//--------------------------------------------------------------------------
const Valeur Isol_hor::boundary_nn_Dir_kk()const {

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
    
  nn_bound = 1 ;  // Juste pour affecter dans espace des configs ;
  

  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      nn_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  nn_bound.std_base_scal() ;
  
  return  nn_bound ;

}


const Valeur Isol_hor::boundary_nn_Neu_kk() const {
  
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
  double rho = 25. ;
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

  Scalar tmp(mp) ;
  tmp = cc - 1 ;
  
  // We  have substracted 1, since we solve for zero condition at infinity 
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


const Valeur Isol_hor::boundary_nn_Dir_lapl() const {

    Scalar det_q = gam_dd()(2,2) * gam_dd()(3,3) - 
	                    gam_dd()(2,3) * gam_dd()(2,3) ;
    Scalar square_q (pow(det_q, 0.5)) ;
    square_q.std_spectral_base() ;

    Sym_tensor qq_uu (mp, CON, mp.get_bvect_spher()) ;
    for (int i=1 ; i<=3 ; i++)
	for (int j=i ; j<=3 ; j++)
	    qq_uu.set(i,j).annule_hard() ;
    
    qq_uu.set(2,2) = gam_dd()(2,2) / det_q ;
    qq_uu.set(3,2) = - gam_dd()(3,2) / det_q ;
    qq_uu.set(3,3) = gam_dd()(3,3) / det_q ;
    qq_uu.std_spectral_base() ;

    Metric met_q (qq_uu) ;

    Sym_tensor hh_uu (mp, CON, mp.get_bvect_spher()) ;
    hh_uu = qq_uu * square_q ;
    hh_uu.set(2,2) = hh_uu(2,2) - 1. ;
    hh_uu.set(3,3) = hh_uu(3,3) - 1. ;

    Vector temp_vect (contract(hh_uu, 1, nn().derive_cov(ff) / nn(), 0)) ;
    temp_vect.set(2).mult_sint() ; 
    
    Scalar temp_scal (mp) ;
    temp_scal = temp_vect(2).dsdt() ;
    temp_scal.div_sint() ;
    temp_scal += temp_vect(3).stdsdp() ;
    temp_scal.div_r() ;    


    Vector temp_vect2 (mp, CON, mp.get_bvect_spher()) ;
    temp_vect2 = contract(k_dd(), 1, gam().radial_vect(), 0).up(0, gam())  -
	contract(k_dd(), 0, 1, gam().radial_vect()*gam().radial_vect(), 0, 1) *
	gam().radial_vect() ;
    temp_vect2.set(2).mult_sint() ;
    temp_vect2 = temp_vect2 * square_q ;

    Scalar temp_scal2 (mp) ;
    temp_scal2 = temp_vect2(2).dsdt() ;
    temp_scal2.div_sint() ;
    temp_scal2 += temp_vect2(3).stdsdp() ;
    temp_scal2.div_r() ;
    
    Scalar lew_pal (mp) ;
    lew_pal = 0. ;
    lew_pal *= square_q ;
    
    Scalar source (temp_scal2 - temp_scal + lew_pal) ;
 
    Scalar logn (mp) ;
    Param bidon ;
    mp.poisson_angu(source, bidon, logn) ;

    double cc = 1. ; // Integration constant
    Scalar lapse (exp(logn)*cc) ;
    lapse.std_spectral_base() ;
    
    lapse = lapse - 1. ;

    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;
    
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
 
  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_r (mp.get_mg()->get_angu()) ;
    
  bnd_beta_r = 1 ;  // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_r.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  bnd_beta_r.std_base_scal() ;
  
  return  bnd_beta_r ;


}


// Component theta of boundary value of beta (using expression in terms 
// of radial vector)
// ------------------------------------------

const Valeur Isol_hor::boundary_beta_theta()const {
  
  Scalar tmp(mp) ;  
  
  tmp = nn() * radial_vect_hor()(2) ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_theta (mp.get_mg()->get_angu()) ;
    
  bnd_beta_theta = 1 ;   // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_theta.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  bnd_beta_theta.std_base_scal() ;
  
  return  bnd_beta_theta ;


}
 
// Component phi of boundary value of beta (using expression in terms 
// of radial vector) 
// -----------------------------------------------------------

const Valeur Isol_hor::boundary_beta_phi()const {

  Scalar tmp (mp) ;

  Scalar ang_vel(mp) ;
  ang_vel = omega_hor() ;
  ang_vel.std_spectral_base() ;
  ang_vel.mult_rsint() ;

  tmp = nn() * radial_vect_hor()(3)  -  ang_vel ;

  int nnp = mp.get_mg()->get_np(1) ;
  int nnt = mp.get_mg()->get_nt(1) ;

  Valeur bnd_beta_phi (mp.get_mg()->get_angu()) ;
    
  bnd_beta_phi = 1 ; // Juste pour affecter dans espace des configs ;
  
  for (int k=0 ; k<nnp ; k++)
    for (int j=0 ; j<nnt ; j++)
      bnd_beta_phi.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
  
  bnd_beta_phi.std_base_scal() ;
  
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

    Vector tmp_vect = nn() * radial_vect_hor() ;
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
	  lim_x.set(0, k, j, 0) =  - aligne * om * ya_mtbl(1, k, j, 0) * (1 + 
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

    Vector tmp_vect = nn() * radial_vect_hor() ;
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
	  lim_y.set(0, k, j, 0) =  aligne * om * xa_mtbl(1, k, j, 0) *(1 +
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

    Vector tmp_vect = nn() * radial_vect_hor() ;
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

  Vector s_tilde =  tradial_vect_hor() ;
  
  Scalar hh_tilde = contract(s_tilde.derive_cov(met_gamt), 0, 1) ;
  hh_tilde.dec_dzpuis(2) ;

  Vector tmp_vect = b_tilde() * s_tilde ;  

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

  Scalar bc_source (mp) ;
  bc_source = 0. ;
  bc_source.std_spectral_base() ;
  bc_source.inc_dzpuis(2) ;

  // beta^r component
  //-----------------
  double rho = 5. ; // rho>1 ; rho=1 "pure Dirichlet" version
  Scalar beta_r_corr = (rho - 1) * b_tilde() * hh_tilde;
  beta_r_corr.inc_dzpuis(2) ;

  Scalar beta_r (mp) ;
  beta_r = 2 * contract(s_tilde, 0, b_tilde().derive_cov(ff), 0) 
      + beta_r_corr - bc_source ;
  beta_r = beta_r / (hh_tilde * rho) ;
    
  beta_r.dec_dzpuis(2) ;
  
  beta_r = beta_r - beta()(2)*s_tilde(2) -  beta()(3)*s_tilde(3) ;
  Vector tmp = s_tilde.down(0, met_gamt) ;
  beta_r = beta_r/tmp(1) ;

  beta_r.set_spectral_va().set_base(beta()(1).get_spectral_va().get_base()) ;
  
  tmp_vect.set(1) = beta_r ;
  tmp_vect.set(3) += angular(3) ; //ang_vel ;


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

