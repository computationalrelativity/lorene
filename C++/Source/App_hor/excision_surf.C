/*
 *  Definition of methods for the class Spheroid and its subclass App_hor
 *
 */

/*
 *   Copyright (c) 2008  Jose-Luis Jaramillo & Jerome Novak & Nicolas Vasset
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

char excision_surf_C[] = "$Header$" ;

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
#include "excision_surf.h"

//---------------//
//  Constructors //
//--------------// 
 

Excision_surf::Excision_surf(const Scalar& h_in, const Metric& gij, const Sym_tensor& Kij2, const Scalar& ppsi, const Scalar& nn, const Vector& beta, double timestep, int int_nos = 1):
  sph(h_in, gij, Kij2),
  conf_fact(ppsi),
  lapse(nn),
  shift(beta),
  gamij (gij),
  Kij(Kij2),
  delta_t(timestep),
  no_of_steps(int_nos)

{
  
 set_der_0x0() ;

}





//Copy constructor//

Excision_surf::Excision_surf (const Excision_surf &exc_in) :sph(exc_in.sph),
							    conf_fact(exc_in.conf_fact),
							    lapse(exc_in.lapse),                          				
							    shift(exc_in.shift),
							    gamij (exc_in.gamij),
							    Kij (exc_in.Kij),
							    delta_t(exc_in.delta_t),
							    no_of_steps(exc_in.no_of_steps)
							    
{
  set_der_0x0() ; 
  
}
//------------//
//Destructor //
//-----------//

Excision_surf::~Excision_surf()
{
  del_deriv() ;
}

// -----------------//
// Memory management//
//------------------//
void Excision_surf::del_deriv() const {
  if (p_get_BC_conf_fact_1 != 0x0) delete p_get_BC_conf_fact_1 ;
  if (p_get_BC_lapse_1 != 0x0) delete p_get_BC_lapse_1 ;
  if (p_get_BC_shift_1 != 0x0) delete p_get_BC_shift_1 ;
  if (p_get_BC_Npsi_1 != 0x0) delete p_get_BC_Npsi_1 ;
  if (p_get_BC_conf_fact_2 != 0x0) delete p_get_BC_conf_fact_2 ;
  if (p_get_BC_conf_fact_3 != 0x0) delete p_get_BC_conf_fact_3 ;
  if (p_get_BC_conf_fact_4 != 0x0) delete p_get_BC_conf_fact_4 ;
  if (p_get_BC_lapse_2 != 0x0) delete p_get_BC_lapse_2 ;
  if (p_get_BC_shift_2 != 0x0) delete p_get_BC_shift_2 ;
  if (p_get_BC_Npsi_2 != 0x0) delete p_get_BC_Npsi_2 ;
  set_der_0x0() ;
}

void Excision_surf::set_der_0x0() const {
  p_get_BC_conf_fact_1 = 0x0 ;
  p_get_BC_lapse_1 = 0x0 ;
  p_get_BC_shift_1 = 0x0 ;
  p_get_BC_Npsi_1 = 0x0 ;
  p_get_BC_conf_fact_2 = 0x0 ;
  p_get_BC_conf_fact_3 = 0x0 ;
  p_get_BC_conf_fact_4 = 0x0 ;
  p_get_BC_lapse_2 = 0x0 ;
  p_get_BC_shift_2 = 0x0 ;
  p_get_BC_Npsi_2 = 0x0 ;

} 


 
//---------//
//Accessors//
//---------//


// Source for the Neumann BC on the conformal factor
const Scalar& Excision_surf::get_BC_conf_fact_1() const{
  if (p_get_BC_conf_fact_1 == 0x0){
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
    p_get_BC_conf_fact_1 = new Scalar(bound_psi);
    // WARNING: Implement Spheroid.expansion instead, only multiplying by \psi^{3} (check again)
    //Generalize to non-zero expansions
}
  return *p_get_BC_conf_fact_1 ;
}


// Source for the Dirichlet BC on the conformal factor, based on a parabolic driver
const Scalar& Excision_surf::get_BC_conf_fact_2(double c_psi_lap, double c_psi_fin, Scalar& expa_fin) const{
  if (p_get_BC_conf_fact_2 == 0x0){

    // Definition of ff
    // ================

    // Start Mapping interpolation
    Scalar thetaplus = sph.theta_plus();
      Scalar theta_plus3 (lapse.get_mp()); 
 
  theta_plus3.allocate_all();
  theta_plus3.std_spectral_base();

      Scalar expa_fin3 (lapse.get_mp()); 
 
  expa_fin3.allocate_all();
  expa_fin3.std_spectral_base();

  int nz = (*lapse.get_mp().get_mg()).get_nzone();
  int nr = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);


  for (int f= 0; f<nz; f++)
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	for (int l=0; l<nr; l++) {
		
	  theta_plus3.set_grid_point(f,k,j,l) = thetaplus.val_grid_point(0,k,j,0);
	  expa_fin3.set_grid_point(f,k,j,l) = expa_fin.val_grid_point(0,k,j,0);
	
	}
      }
  if (nz >2){
     theta_plus3.annule_domain(0);
    theta_plus3.annule_domain(nz - 1);
     expa_fin3.annule_domain(0);
    expa_fin3.annule_domain(nz - 1);
  }
    // End Mapping interpolation
  

  Scalar ff = lapse*(c_psi_lap*theta_plus3.lapang() + c_psi_fin*(theta_plus3 - expa_fin3));
  ff.std_spectral_base();

  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of the expansion, for Runge-Kutta 2nd order scheme
  Scalar psi_int = conf_fact + k_1; psi_int.std_spectral_base();
  Sym_tensor gamconfcov = gamij.cov()/pow(psi_int, 4); // think about the consistency of redifining the conformal metric
                                                       // since in this manner the unimodular conditions is lost
  gamconfcov.std_spectral_base();
  Metric gamconf(gamconfcov);
  Vector tilde_s = gamconf.radial_vect();
  Scalar theta_int =   ((1./psi_int)*contract((contract(Kij,1,tilde_s,0)),0, tilde_s,0));    
  theta_int += psi_int*tilde_s.divergence(gamconf);
  theta_int += 4.*contract(psi_int.derive_cov(gamconf),0,tilde_s,0);
  theta_int = theta_int/pow(psi_int,3);
  theta_int.std_spectral_base();
 
  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse*(c_psi_lap*theta_int.lapang() + c_psi_fin*(theta_int - expa_fin3));
 
  // Definition of k_2
  Scalar k_2 = delta_t*ff_int; 
  k_2.std_spectral_base();
 
  // Result of RK2 evolution
  if (nz >2){
    k_2.annule_domain(nz-1);}
  Scalar bound_psi = conf_fact + k_2;
  bound_psi.std_spectral_base();
  bound_psi.set_spectral_va().ylm();
  
  // Assignment of output
  p_get_BC_conf_fact_2 = new Scalar(bound_psi);

}
  return *p_get_BC_conf_fact_2 ;
}




// Source for the Neumann BC on the conformal factor, based on a parabolic driver for the expansion
const Scalar& Excision_surf::get_BC_conf_fact_3(double c_theta_lap, double c_theta_fin, Scalar& expa_fin) const{
  if (p_get_BC_conf_fact_3 == 0x0){

    // Definition of ff
    // ================

    // Start Mapping interpolation
    Scalar thetaplus = sph.theta_plus();
      Scalar theta_plus3 (lapse.get_mp()); 
 
  theta_plus3.allocate_all();
  theta_plus3.std_spectral_base();

      Scalar expa_fin3 (lapse.get_mp()); 
 
  expa_fin3.allocate_all();
  expa_fin3.std_spectral_base();

  int nz = (*lapse.get_mp().get_mg()).get_nzone();
  int nr = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);


  for (int f= 0; f<nz; f++)
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	for (int l=0; l<nr; l++) {
		
	  theta_plus3.set_grid_point(f,k,j,l) = thetaplus.val_grid_point(0,k,j,0);
	  expa_fin3.set_grid_point(f,k,j,l) = expa_fin.val_grid_point(0,k,j,0);
	
	}
      }
  if (nz >2){
     theta_plus3.annule_domain(0);
    theta_plus3.annule_domain(nz - 1);
     expa_fin3.annule_domain(0);
    expa_fin3.annule_domain(nz - 1);
  }
    // End Mapping interpolation
  

  Scalar ff = lapse*(c_theta_lap*theta_plus3.lapang() + c_theta_fin*(theta_plus3 - expa_fin3));
  ff.std_spectral_base();

  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of the expansion, for Runge-Kutta 2nd order scheme
  Scalar theta_int = theta_plus3 + k_1; theta_int.std_spectral_base();

  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse*(c_theta_lap*theta_int.lapang() + c_theta_fin*(theta_int - expa_fin3));
 
  // Definition of k_2
  Scalar k_2 = delta_t*ff_int; 
  k_2.std_spectral_base();
 
  // Result of RK2 evolution
  Scalar bound_theta = theta_plus3 + k_2;
  bound_theta.std_spectral_base();
  
  // Calculation of Neumann BC for Psi
  Scalar bound_psi = get_BC_conf_fact_1() + bound_theta*pow(conf_fact,3)/4.;
  bound_psi.std_spectral_base();
  bound_psi.set_spectral_va().ylm();

  // Assignment of output
  p_get_BC_conf_fact_3 = new Scalar(bound_psi);

}
  return *p_get_BC_conf_fact_3 ;
}

// Source for the Dirchlet BC on the conformal factor, based on the consistency condition derived from the trace
// of the kinematical relation
const Scalar& Excision_surf::get_BC_conf_fact_4() const{
  if (p_get_BC_conf_fact_3 == 0x0){

    // Definition of ff
    // ================
    const Metric_flat& flat = lapse.get_mp().flat_met_spher() ; 
    

    Scalar ff = contract(shift, 0, conf_fact.derive_cov(flat), 0) 
             + 1./6. * conf_fact * (shift.divergence(flat)) ; // Add he N K term
    // Divergence with respect to the conformal metric coincides with divergence withh respect to the
    // flat metric (from the unimodular condition on the conformal metric)
    // In this way, we do not need to recalculate a conformal metric in the intermediate RK step that would violated
    // the unimodular condition

    ff.std_spectral_base() ;

    // Definition of k_1
    Scalar k_1 =delta_t*ff;
    
    // Intermediate value of the expansion, for Runge-Kutta 2nd order scheme
    Scalar psi_int = conf_fact + k_1; psi_int.std_spectral_base();
    
    // Recalculation of ff with intermediate values. 
    Scalar ff_int =  contract(shift, 0, psi_int.derive_cov(flat), 0) 
                    + 1./6. * conf_fact * (shift.divergence(flat)) ;  // Add he N K term

 
    // Definition of k_2
    Scalar k_2 = delta_t*ff_int; 
    k_2.std_spectral_base();

    // Result of RK2 evolution
    Scalar bound_psi = conf_fact + k_2;
    bound_psi.std_spectral_base();
    bound_psi.set_spectral_va().ylm();

    // Assignment of output
    p_get_BC_conf_fact_4 = new Scalar(bound_psi);

}
  return *p_get_BC_conf_fact_4 ;
}




// Source for the Dirichlet BC on the lapse
const Scalar& Excision_surf::get_BC_lapse_1(double value) const{
  if (p_get_BC_lapse_1 == 0x0){

    Scalar bound_lapse(lapse.get_mp()); 
    bound_lapse = value; bound_lapse.std_spectral_base();
    bound_lapse.set_spectral_va().ylm();
    p_get_BC_lapse_1 = new Scalar(bound_lapse); 
    
}
  return *p_get_BC_lapse_1 ;
}

// Source for Dirichlet BC on the lapse, based on a parabolic driver
const Scalar& Excision_surf::get_BC_lapse_2(double lapse_fin, double c_lapse_lap, double c_lapse_fin) const{
  if (p_get_BC_lapse_2 == 0x0){

  
  
  Scalar ff = lapse*(c_lapse_lap*lapse.lapang() + c_lapse_fin*lapse);
  ff.std_spectral_base();
  ff += -lapse*c_lapse_fin*lapse_fin;
  ff.std_spectral_base();


  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of Npsi, for Runge-Kutta 2nd order scheme
  Scalar lapse_int = lapse + k_1; lapse_int.std_spectral_base();

  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse*(c_lapse_lap*lapse_int.lapang() + c_lapse_fin*lapse_int);
  ff_int.std_spectral_base();
  ff_int += -lapse*c_lapse_fin*lapse_fin;
  ff_int.std_spectral_base();
  
 
  // Definition of k_2
  Scalar k_2 = delta_t*ff_int; 
  k_2.std_spectral_base();
 
  // Result of RK2 evolution
  Scalar bound_lapse = lapse + k_2;
  bound_lapse.std_spectral_base();
  bound_lapse.set_spectral_va().ylm();

  
        p_get_BC_lapse_2 = new Scalar(bound_lapse); 
    
}
  return *p_get_BC_lapse_2 ;
}

// Source for the Dirichlet BC on the shift
const Vector& Excision_surf::get_BC_shift_1(double Omega) const{
  if (p_get_BC_shift_1 == 0x0){

    // Radial vector for the full 3-metric.
    Vector sss = gamij.radial_vect();

    // Boundary value for the radial part of the shift
    Scalar bound = lapse ; 
 

      // Tangent part of the shift
      // (For now, only axisymmetric configurations are envisaged)
  
      Vector V_par = shift;
      Scalar V_phi = lapse; V_phi.annule_hard(); V_phi = 1.; // Rotation parameter for spacetime
      V_phi.std_spectral_base() ; V_phi.mult_rsint();
      V_par.set(1).annule_hard();
      V_par.set(2).annule_hard();
      V_par.set(3) = V_phi;

      V_par = V_par*Omega;

  
      // Construction of the total shift boundary condition
      Vector bound_shift = bound*sss + V_par;
      bound_shift.std_spectral_base(); // Boundary is fixed by value of 3 components of a vector (rather than value of potentials) 
      p_get_BC_shift_1 = new Vector(bound_shift);
}
  return *p_get_BC_shift_1 ;
}


// Source for the Dirichlet BC on the shift, based on a Parabolic driver
const Vector& Excision_surf::get_BC_shift_2(double c_bb_lap, double c_bb_fin, double c_V_lap ) const{
  if (p_get_BC_shift_2 == 0x0){

    // Radial vector for the full 3-metric.
     Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
     
//     // Boundary value for the radial part of the shift: parabolic driver for (b-N)
     //  Scalar bound = lapse ; 
     Scalar bb = contract (shift,0, sss_down,0); 
   
  Scalar b_min_N = bb - lapse;
  Scalar ff = lapse*(c_bb_lap*b_min_N.lapang() + c_bb_fin*b_min_N);
  ff.std_spectral_base();
 

  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of b-N, for Runge-Kutta 2nd order scheme
  Scalar b_min_N_int = b_min_N + k_1; b_min_N_int.std_spectral_base();

  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse*(c_bb_lap*b_min_N_int.lapang() + c_bb_fin*b_min_N_int);
 
  // Definition of k_2
  Scalar k_2 = delta_t*ff_int; 
  k_2.std_spectral_base();
 
  // Result of RK2 evolution
  Scalar bound_b_min_N =  b_min_N + k_2;
  bound_b_min_N.std_spectral_base();
  bound_b_min_N.set_spectral_va().ylm();

  Scalar bb2 = bound_b_min_N + lapse; // Look out for additional term for time variation of lapse and sss
  

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
       Vector bound_shift = bb2*sss + bound_V;
       bound_shift.std_spectral_base();
       p_get_BC_shift_2 = new Vector(bound_shift);
}
  return *p_get_BC_shift_2 ;
}

// Source for the Dirichlet BC on (N*Psi1)
const Scalar& Excision_surf::get_BC_Npsi_1(double value) const{
  if (p_get_BC_Npsi_1 == 0x0){

    Scalar bound_Npsi = value*conf_fact;
    bound_Npsi.set_spectral_va().ylm();
    p_get_BC_Npsi_1 = new Scalar(bound_Npsi);
    
}
  return *p_get_BC_Npsi_1 ;
}

// Source for the Dirichlet BC on (N*Psi1), based on a parabolic driver.
const Scalar& Excision_surf::get_BC_Npsi_2(double npsi_fin, double c_npsi_lap, double c_npsi_fin) const{
  if (p_get_BC_Npsi_2 == 0x0){


    Scalar npsi = lapse*conf_fact; npsi.std_spectral_base();
  Scalar ff = lapse*(c_npsi_lap*npsi.lapang() + c_npsi_fin*npsi);
  ff.std_spectral_base();
  ff += -lapse*c_npsi_fin*npsi_fin;
  ff.std_spectral_base();


  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of Npsi, for Runge-Kutta 2nd order scheme
  Scalar npsi_int = npsi + k_1; npsi_int.std_spectral_base();

  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse*(c_npsi_lap*npsi_int.lapang() + c_npsi_fin*(npsi_int - npsi_fin));
 
  // Definition of k_2
  Scalar k_2 = delta_t*ff_int; 
  k_2.std_spectral_base();
 
  // Result of RK2 evolution
  Scalar bound_npsi = npsi + k_2;
  bound_npsi.std_spectral_base();
  bound_npsi.set_spectral_va().ylm();



//     Scalar bound_Npsi = value*conf_fact;
//     bound_Npsi.set_spectral_va().ylm();
     p_get_BC_Npsi_2 = new Scalar(bound_npsi);
    
}
  return *p_get_BC_Npsi_2 ;
}


 void Excision_surf::sauve(FILE* ) const {

   cout << "c'est pas fait!" << endl ;
   return ; 

 }
