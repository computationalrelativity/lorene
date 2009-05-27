/*
 *  Definition of methods for the class Spheroid and its subclass App_hor
 *
 */

/*
 *   Copyright (c) 2009  Jose-Luis Jaramillo & Jerome Novak & Nicolas Vasset
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 
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

char excision_hor_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "metric.h"
#include "spheroid.h"
#include "utilitaires.h"
#include "param.h"
#include "itbl.h"
#include "map.h"
#include <assert.h>
#include "nbr_spx.h"
#include "math.h"
#include "param_elliptic.h"
#include "tensor.h"
#include "sym_tensor.h"
#include "diff.h"
#include "proto.h"
#include "param.h" 
#include "excision_hor.h"

//---------------//
//  Constructors //
//--------------// 
 

Excision_hor::Excision_hor(const Scalar& h_in, const Metric& gij, const Sym_tensor& Kij2, const Scalar& ppsi, const Scalar& nn, const Vector& beta, const Sym_tensor& Tij2, double timestep, int int_nos):
  sph(h_in, gij, Kij2),
  conf_fact(ppsi),
  lapse(nn),
  shift(beta),
  gamij (gij),
  Kij(Kij2),
  delta_t(timestep),
  no_of_steps(int_nos),
  Tij(Tij2)
{
  
 set_der_0x0() ;

}





//Copy constructor//

Excision_hor::Excision_hor(const Excision_hor &exc_in) :sph(exc_in.sph),
							    conf_fact(exc_in.conf_fact),
							    lapse(exc_in.lapse),                          				
							    shift(exc_in.shift),
							    gamij (exc_in.gamij),
							    Kij (exc_in.Kij),
							    delta_t(exc_in.delta_t),
							no_of_steps(exc_in.no_of_steps),
							Tij(exc_in.Tij)
							    
{
  set_der_0x0() ; 
  
}
//------------//
//Destructor //
//-----------//

Excision_hor::~Excision_hor()
{
  del_deriv() ;
}

// -----------------//
// Memory management//
//------------------//
void Excision_hor::del_deriv() const {
  
 if (p_get_BC_conf_fact != 0x0) delete p_get_BC_conf_fact ;
 if (p_get_BC_bmN != 0x0) delete p_get_BC_bmN ;
 if (p_get_BC_bpN != 0x0) delete p_get_BC_bpN ;
 if (p_get_BC_shift != 0x0) delete p_get_BC_shift ;
  set_der_0x0() ;
}

void Excision_hor::set_der_0x0() const {
  p_get_BC_conf_fact = 0x0 ;
  p_get_BC_bmN = 0x0 ;
  p_get_BC_bpN = 0x0 ;
  p_get_BC_shift = 0x0 ;

} 


 
//---------//
//Accessors//
//---------//



// Source for the Neumann BC on the conformal factor
const Scalar& Excision_hor::get_BC_conf_fact() const{
  if (p_get_BC_conf_fact == 0x0){
    Sym_tensor gamconfcov = gamij.cov()/pow(conf_fact, 4); 
    gamconfcov.std_spectral_base();
    Metric gamconf(gamconfcov);
    Vector tilde_s = gamconf.radial_vect();
    Scalar bound_psi =   -((1./conf_fact)*contract((contract(Kij,1,tilde_s,0)),0, tilde_s,0));    
    bound_psi += -conf_fact*tilde_s.divergence(gamconf);
    bound_psi = 0.25*bound_psi;
    bound_psi += -contract(conf_fact.derive_cov(gamconf),0,tilde_s,0) + conf_fact.dsdr();
    bound_psi.std_spectral_base();
    bound_psi.set_spectral_va().ylm();    
    p_get_BC_conf_fact = new Scalar(bound_psi);

}
  return *p_get_BC_conf_fact ;
}



// Case 0: Source of Dirichlet BC for (b-N), based on an entropy prescription.
// WARNING: the argument value has to be carefully fixed w.r.t initial data for (attempted) continuity.

//  Case 1: Source of Dirichlet BC for (b-N), from a component of projected Einstein Equations.
//  Requires a 2d poisson solver for a non-flat metric.


const Scalar& Excision_hor::get_BC_bmN(int choice_bmN, double value) const{
  if (p_get_BC_bmN == 0x0){

    switch(choice_bmN){

    case 0 : {

  Scalar thetaminus = sph.theta_minus();
      Scalar theta_minus3 (lapse.get_mp()); 
 
  theta_minus3.allocate_all();
  theta_minus3.std_spectral_base();

  int nz = (*lapse.get_mp().get_mg()).get_nzone();
  int nr = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);


  for (int f= 0; f<nz; f++)
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	for (int l=0; l<nr; l++) {
		
	  theta_minus3.set_grid_point(f,k,j,l) = thetaminus.val_grid_point(0,k,j,0);
	
	}
      }
  if (nz >2){
     theta_minus3.annule_domain(0);
    theta_minus3.annule_domain(nz - 1);
   }


    Scalar bound_bmN(lapse.get_mp()); 
    bound_bmN = - value*theta_minus3; bound_bmN.std_spectral_base();
    bound_bmN.set_spectral_va().ylm();
    p_get_BC_bmN = new Scalar(bound_bmN); 
    }

    case 1 : {

    Scalar bound_bmN(lapse.get_mp());
    bound_bmN.allocate_all();
    bound_bmN.std_spectral_base();

// Radial vector for the full 3-metric.
     Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
      Scalar bb = contract (shift,0, sss_down,0);
      Scalar bmN3 = bb - lapse; bmN3.set_spectral_va().ylm();
      Scalar bpN3 = bb + lapse; bpN3.set_spectral_va().ylm();  
      
  int nz = (*lapse.get_mp().get_mg()).get_nzone();
  int nr = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);
    
    Scalar bmN(sph.get_hsurf().get_mp());
    bmN.allocate_all();
    bmN.std_spectral_base();
    bmN.set_spectral_va().ylm();
     Scalar bpN(sph.get_hsurf().get_mp());
    bpN.allocate_all();
    bpN.std_spectral_base();
    bpN.set_spectral_va().ylm();

    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {

       	  bmN.set_grid_point(0,k,j,0) = bmN3.val_grid_point(1,k,j,0);
	  bpN.set_grid_point(0,k,j,0) = bpN3.val_grid_point(1,k,j,0);
	}

    Scalar bmN_new(bmN.get_mp());
    bmN_new.allocate_all();
    bmN_new.std_spectral_base();
     
    double diff_ent = 1.; 
    double precis = 1.e-9;
    int mer_max = 200;
    double relax = 0.;
    for(int mer=0 ;(diff_ent > precis) && (mer<mer_max) ; mer++) {    
  
      // Calculation of some source terms.
     
      Scalar hsurf = sph.get_hsurf();
      hsurf.set_spectral_va().ylm();
      const Metric_flat& fmets = hsurf.get_mp().flat_met_spher() ;
 
    Scalar shear_up = sph.shear(); shear_up.up_down(sph.get_qab());
 
    Scalar B_source = 0.5*contract(contract(sph.shear(),0, shear_up, 0),0,1) + 4.*M_PI*Tij.trace(sph.get_qab()); // Redo the matter terms.
    Scalar A_source = 0.5*sph.get_ricci() - contract(sph.derive_cov2d(sph.get_ll()), 0, 1) - contract(sph.get_ll(),0, sph.get_ll().up_down(sph.get_qab()),0) - 8.*M_PI*Tij.trace(sph.get_qab()); // Redo the matter terms.
   
    Scalar op_bmN_add = - 2.*contract(sph.derive_cov2d(bmN),0, sph.get_ll(),0) + A_source*bmN;

    Scalar source_bmN = B_source*bpN - op_bmN_add;  
    source_bmN.set_spectral_va().ylm();

    Scalar sqrtqh2 = sph.sqrt_q()*hsurf*hsurf;
    sqrtqh2.set_spectral_va().ylm();
   
    source_bmN = sqrtqh2*source_bmN;
     
     // Conformal decomposition of the 2-metric
        
    Sym_tensor qab_con = sph.get_qab().con();
    qab_con = qab_con/(hsurf*hsurf); // Renormalization due to the triad still not built-in spheroid class
    //This is provisory work.

    // h^ab as q^ab = (f^ab + h^ab) / sqrt_q
    Sym_tensor hab =(qab_con*sqrtqh2 - fmets.con()) / (hsurf*hsurf) ;
    // for the sake of clarity
    hab.set(1,1) = 1. ;
    hab.set(1,2) = 0. ;
    hab.set(1,3) = 0. ;
    hab.std_spectral_base() ;
    //end
    // Complete source for the angular laplacian.
    Scalar d_bmN = sph.derive_cov2dflat(bmN);
    d_bmN.set_spectral_va().ylm();
      Scalar d2_bmN = sph.derive_cov2dflat(d_bmN);
    d2_bmN.set_spectral_va().ylm();
 
    Scalar source_add = - hsurf*hsurf*contract(hab, 0,1, d2_bmN, 0,1) + sqrtqh2*contract(contract(qab_con,0,1,sph.delta(),1,2),0,d_bmN,0) ;   
    source_add.set_spectral_va().ylm();
    source_bmN = source_bmN + source_add;
    //
    
    // System inversion
    bmN_new = source_bmN.poisson_angu(0.);
    
    // Actualisation of the principal variable, convergence parameter.
    diff_ent = max(maxabs(bmN - bmN_new));
    
    bmN = relax*bmN + (1. - relax)*bmN_new;   
     
    }
      bound_bmN = bmN;
    bound_bmN.set_spectral_va().ylm();
          p_get_BC_bmN = new Scalar(bound_bmN); 
    }

}
  }
  return *p_get_BC_bmN ;

}


// Case 0:  Arbitrary Dirichlet BC for (b+N), fixed by a parabolic driver towards a constant value.
//  Case 1: Source of Dirichlet BC for (b+N), from a component of projected Einstein Equations.
const Scalar& Excision_hor::get_BC_bpN(int choice_bpN, double c_bpn_lap, double c_bpn_fin, Scalar *bpN_fin) const{
  if (p_get_BC_bpN == 0x0){

    switch(choice_bpN) {
 
    case 0 : {

      Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
      Scalar bb = contract (shift,0, sss_down,0);
  Scalar bpN = bb + lapse;
  Scalar ff = lapse*(c_bpn_lap*bpN.lapang() + c_bpn_fin*(bpN- *bpN_fin));
  ff.std_spectral_base();
 

  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of b-N, for Runge-Kutta 2nd order scheme
  Scalar bpN_int = bpN + k_1; bpN_int.std_spectral_base();

  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse*(c_bpn_lap*bpN_int.lapang() + c_bpn_fin*bpN_int);
 
  // Definition of k_2
  Scalar k_2 = delta_t*ff_int; 
  k_2.std_spectral_base();
 
  // Result of RK2 evolution
  Scalar bound_bpN =  bpN + k_2;
  bound_bpN.std_spectral_base();
  bound_bpN.set_spectral_va().ylm();

        p_get_BC_bpN = new Scalar(bound_bpN); 
    

    }

    case 1 : {

    
  int nz = (*lapse.get_mp().get_mg()).get_nzone();
  int nr = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);

    
    
  Scalar bmN3 = get_BC_bmN(0, 1.); // change the argument.
    
    Scalar bmN(sph.get_hsurf().get_mp());
    bmN.allocate_all();
    bmN.std_spectral_base();
    bmN.set_spectral_va().ylm();
 
    

    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {

		
	  bmN.set_grid_point(0,k,j,0) = bmN3.val_grid_point(1,k,j,0);
	
	}


    Scalar bound_bpN(lapse.get_mp());
    bound_bpN.allocate_all();
    bound_bpN.std_spectral_base();
 
    // Definition of source terms in relation (6) of Jaramillo et al. 2007.
    // All is done on the spheroid of radius r=1.
 
    Scalar shear_up = sph.shear(); shear_up.up_down(sph.get_qab());
 
    Scalar B_source = 0.5*contract(contract(sph.shear(),0, shear_up, 0),0,1) + 4.*M_PI*Tij.trace(sph.get_qab()); // Redo the matter terms.
    Scalar A_source = 0.5*sph.get_ricci() - contract(sph.derive_cov2d(sph.get_ll()), 0, 1) - contract(sph.get_ll(),0, sph.get_ll().up_down(sph.get_qab()),0) - 8.*M_PI*Tij.trace(sph.get_qab()); // Redo the matter terms.
  
    // Curved 2d Laplacian of (b -N).
   
    Sym_tensor interlap = sph.derive_cov2d(sph.derive_cov2d(bmN)); 
    interlap.up(0,sph.get_qab()); 
    Sym_tensor lap_bmN = contract(interlap,0,1);
 
    Scalar op_bmN = lap_bmN - 2.*contract(sph.derive_cov2d(bmN),0, sph.get_ll(),0) + A_source*bmN;
          
   Scalar  bound_bpN2 = op_bmN/B_source;
    bound_bpN2.std_spectral_base();
    bound_bpN2.set_spectral_va().ylm();


  for (int f= 0; f<nz; f++)
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	for (int l=0; l<nr; l++) {
		
	  bound_bpN.set_grid_point(f,k,j,l) = bound_bpN2.val_grid_point(0,k,j,0);
	
	}
      }
  if (nz >2){
     bound_bpN.annule_domain(0);
    bound_bpN.annule_domain(nz - 1);
   }
 
 

        p_get_BC_bpN = new Scalar(bound_bpN); 
    
 
}

}
}
  return *p_get_BC_bpN ;
}




// Source for the Dirichlet BC on the shift
// The tangential shift is fixed using a parabolic driver based on the conformal Killing equation in the dynamical case.

const Vector& Excision_hor::get_BC_shift( double c_V_lap ) const{
  if (p_get_BC_shift == 0x0){

    // Radial vector for the full 3-metric.
     Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
     
//     // Boundary value for the radial part of the shift: parabolic driver for (b-N)
     //  Scalar bound = lapse ; 
     Scalar bb = 0.5*(*p_get_BC_bpN + *p_get_BC_bmN) ; // TO CHANGE: additional function?-> put choice-bb
  

  // Tangent part of the shift, with parabolic driver
  
  
  Vector V_par = shift - bb*sss;
  Sym_tensor q_upup = gamij.con() - sss*sss;

  
  // Calculation of the conformal 2d laplacian of V
  Tensor q_updown = q_upup.down(1, gamij); 
  Tensor dd_V = V_par.derive_con(gamij);
  dd_V = contract(q_updown, 1, contract(q_updown,1 ,dd_V, 0), 1);
  Vector lap_V = contract(q_updown, 1, contract(dd_V.derive_cov(gamij),1,2), 0);
  
  // 3d interpolation of the Ricci scalar on the surface.
  
  Scalar ricci2 = sph.get_ricci();
  
     // Start Mapping interpolation
 
      Scalar ricci3 (lapse.get_mp()); 
 
  ricci3.allocate_all();
  ricci3.std_spectral_base();

  int nz = (*lapse.get_mp().get_mg()).get_nzone();
  int nr = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);


  for (int f= 0; f<nz; f++)
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	for (int l=0; l<nr; l++) {
		
	  ricci3.set_grid_point(f,k,j,l) = ricci2.val_grid_point(0,k,j,0);

	
	}
      }
  if (nz >2){
     ricci3.annule_domain(0);
    ricci3.annule_domain(nz - 1);

  }
    // End Mapping interpolation
  
    // Construction of the Ricci COV tensor on the sphere
 
  Sym_tensor ricci_t = gamij.cov() - sss_down*sss_down;
  ricci_t = 0.5*ricci3*ricci_t;
  ricci_t.std_spectral_base();
 
  Tensor ricci_t_updown = contract(q_upup,0, ricci_t,0); 
  
  // Calculation of ff 

  Vector ffV = c_V_lap*lapse*(lap_V + contract(ricci_t_updown,1, V_par,0));
  ffV.std_spectral_base();


  // Definition of k_1
  Vector k_1V =delta_t*ffV;
   
  // Intermediate value of Npsi, for Runge-Kutta 2nd order scheme
  if (nz >2){
    k_1V.annule_domain(nz-1);
  }                             // Patch to avoid dzpuis problems if existent.
  Vector V_par_int = V_par + k_1V;// V_par_int.std_spectral_base();

  // Recalculation of ff with intermediate values.

  Sym_tensor dd_V_int = V_par_int.derive_con(gamij);
  dd_V_int = contract(q_updown, 1, contract(q_updown,1 ,dd_V_int, 0), 1);
  Vector lap_V_int = contract(q_updown, 1, contract(dd_V_int.derive_cov(gamij),1,2), 0);
 
  Vector ffV_int =  c_V_lap*lapse*(lap_V_int + contract(ricci_t_updown,1, V_par_int,0));
 
  // Definition of k_2
  Vector k_2V = delta_t*ffV_int; 
  //  k_2.std_spectral_base();
 
  // Result of RK2 evolution
  if (nz >2){
    k_2V.annule_domain(nz-1);
  }
  Vector bound_V = V_par + k_2V;
  //  bound_V.std_spectral_base();

       // Construction of the total shift boundary condition
       Vector bound_shift = bb*sss + bound_V;
       bound_shift.std_spectral_base();
       p_get_BC_shift = new Vector(bound_shift);
}
  return *p_get_BC_shift ;
}
 


 void Excision_hor::sauve(FILE* ) const {

   cout << "c'est pas fait!" << endl ;
   return ; 

 }
