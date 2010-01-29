/*
 *  Definition of methods for the class Excision_surf and its subclass App_hor
 *
 */
////////////////////////////////////////////////////////////////////////////////
//                               LOCAL VERSION
//  Several additional functions that may be temporary (on the testing pad...)
///////////////////////////////////////////////////////////////////////////////

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
  no_of_steps(int_nos),
  expa(sph.theta_plus()),
  dt_expa(sph.theta_plus())

{
  
  dt_expa.set_etat_zero();
 
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
							    no_of_steps(exc_in.no_of_steps),
							    expa(exc_in.expa),
							    dt_expa(exc_in.dt_expa)

{
  set_der_0x0() ; 
  
}




// Assignment to another Excision_surf
void Excision_surf::operator=(const Excision_surf& surf_in)
{

  sph = surf_in.sph ;
  conf_fact = surf_in.conf_fact ;
  lapse = surf_in.lapse ;
  shift = surf_in.shift ;
  gamij  = surf_in.gamij ;
  Kij = surf_in.Kij ;
  delta_t = surf_in.delta_t ;
  no_of_steps = surf_in.no_of_steps ;
  expa = surf_in.expa;
  dt_expa = surf_in.dt_expa;
  del_deriv() ;  // Deletes all derived quantities

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
  if (p_get_BC_lapse_3 != 0x0) delete p_get_BC_lapse_3 ;
  if (p_get_BC_lapse_4 != 0x0) delete p_get_BC_lapse_4 ;
  if (p_derive_t_expa != 0x0) delete p_derive_t_expa ;
  if (p_get_BC_shift_2 != 0x0) delete p_get_BC_shift_2 ;
  if (p_get_BC_shift_3 != 0x0) delete p_get_BC_shift_3 ;
  if (p_get_BC_shift_4 != 0x0) delete p_get_BC_shift_4 ;
  if (p_get_BC_Npsi_2 != 0x0) delete p_get_BC_Npsi_2 ;
  if (p_get_BC_Npsi_3 != 0x0) delete p_get_BC_Npsi_3 ;
  if (p_get_BC_Npsi_4 != 0x0) delete p_get_BC_Npsi_4 ;
  if (p_get_BC_Npsi_5 != 0x0) delete p_get_BC_Npsi_5 ;
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
  p_get_BC_lapse_3 = 0x0 ;
  p_get_BC_lapse_4 = 0x0 ;
  p_derive_t_expa = 0x0 ;
  p_get_BC_shift_2 = 0x0 ;
  p_get_BC_shift_3 = 0x0 ;
  p_get_BC_shift_4 = 0x0 ;
  p_get_BC_Npsi_2 = 0x0 ;
  p_get_BC_Npsi_3 = 0x0 ;
  p_get_BC_Npsi_4 = 0x0 ;
  p_get_BC_Npsi_5 = 0x0 ;

} 



			    //--------------//
			    //	  Outputs   //
			    //--------------//



// Provides the three parameters to use for hyperbolic evolution of the expansion,
// Using expansion and its time derivative on initial data.
// WARNING: to be called once and for all on initial data. Re-calling it does not make sense.

void Excision_surf::get_evol_params_from_ID(double alpha, double beta, double gamma, Scalar& Ee, Vector& Jj, Sym_tensor& Ss){
  cout << "===========================================================================================" << endl;
  cout << "Starting the routine that computes parameters for the hyperbolic evolution of the expansion" << endl;
  cout << "Those parameters should be set once and for all: do not call that routine twice" <<endl;
  cout << "===========================================================================================" << endl;

  Scalar expa_0= sph.theta_plus();
  expa_0.set_spectral_va().ylm();

  set_expa() = expa_0;

  if ((max(expa_0))(0)<=1.e-5)
    
    {  
      cout << "=============================================================================" << endl;
      cout << " WARNING: Routine failure..." << endl;
      cout << "the horizon expansion is already pretty close to zero..."<<  endl;
       cout << "We advise to switch from hyperbolic evolution to parabolic evolution" << endl;
       cout << "This routine will not do anything relevant on parameters alpha, beta, gamma" << endl;
     cout << "=============================================================================" << endl;
     double is_expansion_adapted = 1.;
     assert (is_expansion_adapted ==2.);
      
	return;
    }
  else{

  Scalar fbN = derive_t_expa(Ee, Jj, Ss);       // Basic "right" time derivative of the expansion from ID
  set_dt_expa()=fbN;  

  Scalar tf_zero = 2.*fbN/expa_0; tf_zero.std_spectral_base(); 
  tf_zero.set_spectral_va().ylm(); // Relevant quantity for parameter characterization (see notes).
  
  Mtbl_cf* tfz = tf_zero.set_spectral_va().c_cf;// Mtbl_cf storing the spectral coefficients of tf_zero.
  
  // Choice for \beta (see notes for all choices)
  beta = 2.*(*tfz).val_in_bound_jk(0,0,0);

  cout << "beta:   " << beta << endl;
  
  // Choice for \gamma
  gamma = -beta*beta* (1./8.); // Any choice between zero and -beta*beta*(1./4.) is acceptable.

    cout << "gamma:    " << gamma << endl; 
    
  // Choice for \alpha, with an additional assertion.
    alpha = 0.;
    double alpha_aux=0;  
  double nb_l = (*tfz).set(0).get_dim(1);
    double nb_m = (*tfz).set(0).get_dim(2);
      
    cout << "nb_l:   " << nb_l << endl;
    cout << "nb_m:   " << nb_m << endl;
      
    if (nb_m==1){
    alpha =  2.*abs(beta - (*tfz).val_in_bound_jk(0,1,0));
    }
    else{
     for (int ii=0; ii <nb_m; ii++){
     alpha_aux = 2.*abs(beta - (*tfz).val_in_bound_jk(0,1,ii));
     if (alpha_aux >=alpha){
       alpha = alpha_aux;
       }
    }
    }
    cout << "alpha:   " << alpha << endl;
    // Test for higher harmonics, supposedly of lower amplitude, but... 
    for (int ii=2; ii <nb_l; ii++)
      for (int jj=0; jj <nb_m; jj++){
	alpha_aux =  2.*(beta - (*tfz).val_in_bound_jk(0,ii,jj))/(double(ii*(ii+1)));
	   if (alpha_aux >= alpha){
	     cout << "higher Ylm harmonics are predominant in the expansion variation on the horizon!" << endl;
	     cout << "changing the values of the parameters accordingly..." << endl;
	     alpha = alpha_aux;}
      }

  return;                                                                         
}
}

//---------//
//Accessors//
//---------//

// Source for the Neumann BC on the conformal factor, with arbitrary expansion
const Scalar& Excision_surf::get_BC_conf_fact_1(bool isMOTS) const{
 
 
 if (p_get_BC_conf_fact_1 == 0x0){

    //Initialization of expa in the trivial case
    Scalar exppa(sph.theta_plus().get_mp());
    if ( isMOTS == true)
      {
	exppa.set_etat_zero();
	
      } 
    else {
      assert (expa.get_spectral_va().get_etat() !=0);
	exppa = expa;
   
 
// 	expa.spectral_display();
// 	int tgh; cin >> tgh;
    }

   // 3d interpolation for expa
     Scalar expa_3 (lapse.get_mp()); 
 
  expa_3.allocate_all();
  expa_3.std_spectral_base();

  int nz = (*lapse.get_mp().get_mg()).get_nzone();
  int nr = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);


  for (int f= 0; f<nz; f++)
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	for (int l=0; l<nr; l++) {

	  expa_3.set_grid_point(f,k,j,l) = exppa.val_grid_point(0,k,j,0);
	
	}
      }

  if (nz >2){
  
     expa_3.annule_domain(0);
    expa_3.annule_domain(nz - 1);
  }
  expa_3.std_spectral_base();
  expa_3.set_spectral_va().ylm();
 // End Mapping interpolation


    Sym_tensor gamconfcov = gamij.cov()/pow(conf_fact, 4); 
    gamconfcov.std_spectral_base();
    Metric gamconf(gamconfcov);
    Vector tilde_s = gamconf.radial_vect();
    Scalar bound_psi =   -((1./conf_fact)*contract((contract(Kij,1,tilde_s,0)),0, tilde_s,0));    
    bound_psi.annule_domain(nz-1);
    bound_psi += -conf_fact*tilde_s.divergence(gamconf);
    bound_psi.annule_domain(nz-1);
    bound_psi += (conf_fact*conf_fact*conf_fact)*Kij.trace(gamij);
    bound_psi = 0.25*bound_psi;
    bound_psi.annule_domain(nz-1);
    bound_psi = (bound_psi -contract(conf_fact.derive_cov(gamconf),0,tilde_s,0)) + conf_fact.dsdr();
    Scalar bound_psi2 = expa_3*((conf_fact*conf_fact*conf_fact))/4.;
    // Remark: the used expansion term is actually \frac{\theta^{+}}{N}, as it is not normalized by the lapse in the Spheroid class.
    bound_psi2.annule_domain(nz -1);
    bound_psi = bound_psi + bound_psi2;
    bound_psi.std_spectral_base();
    bound_psi.set_spectral_va().ylm();    

  p_get_BC_conf_fact_1 = new Scalar(bound_psi);

  
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
  theta_int += -Kij.trace(gamij);
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
//   cout << "convergence?" << endl;
//   cout << expa_fin3.val_grid_point(1,0,0,0) << endl;
//   cout << theta_plus3.val_grid_point(1,0,0,0) << endl;
  

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
  Scalar bound_psi = get_BC_conf_fact_1(true) + bound_theta*pow(conf_fact,3)/4.;
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
  if (p_get_BC_conf_fact_4 == 0x0){

    // Definition of ff
    // ================
    const Metric_flat& flat = lapse.get_mp().flat_met_spher() ; 
    
 int nz = (*lapse.get_mp().get_mg()).get_nzone();
    Scalar ff = contract(shift, 0, conf_fact.derive_cov(flat), 0) 
      + 1./6. * (conf_fact * (shift.divergence(flat) - lapse*Kij.trace(gamij))) ; // Add he N K term
    // Divergence with respect to the conformal metric coincides with divergence withh respect to the
    // flat metric (from the unimodular condition on the conformal metric)
    // In this way, we do not need to recalculate a conformal metric in the intermediate RK step that would violated
    // the unimodular condition

    ff.std_spectral_base() ;

    // Definition of k_1
    Scalar k_1 =delta_t*ff;
    k_1.annule_domain(nz-1);

    
    // Intermediate value of the expansion, for Runge-Kutta 2nd order scheme
    Scalar psi_int = conf_fact + k_1; psi_int.std_spectral_base();
    psi_int.annule_domain(nz-1);

    // Recalculation of ff with intermediate values. 
    Scalar ff_int =  contract(shift, 0, psi_int.derive_cov(flat), 0) 
                    + 1./6. * psi_int * (shift.divergence(flat) -  lapse*Kij.trace(gamij)) ;  // Add he N K term

 
    // Definition of k_2
    Scalar k_2 = delta_t*ff_int; 
    
    //  k_2 = k_2*1000; //TO REMOVE

    k_2.std_spectral_base();

    k_2.annule_domain(nz-1);
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
  ff = ff -lapse*c_lapse_fin*lapse_fin;
  ff.std_spectral_base();


  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of lapse, for Runge-Kutta 2nd order scheme
  Scalar lapse_int = lapse + k_1; lapse_int.std_spectral_base();

  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse*(c_lapse_lap*lapse_int.lapang() + c_lapse_fin*lapse_int);
  ff_int.std_spectral_base();
  ff_int = ff_int -lapse*c_lapse_fin*lapse_fin;
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




// Source for Dirichlet BC on the lapse, based on einstein equations!! dttheta is 2d scalar, and matter quantities are 3d.
const Scalar& Excision_surf::get_BC_lapse_3(Scalar& dttheta, Scalar& Ee, Vector& Jj, Sym_tensor& Sij, bool sph_sym) const{
  if (p_get_BC_lapse_3 == 0x0){
  int nz2 = (*lapse.get_mp().get_mg()).get_nzone();
    // Radial vector for the full 3-metric.
     Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
     Scalar bb = contract (shift,0, sss_down,0); 

     Scalar bound_N3(lapse.get_mp()); bound_N3.allocate_all(); bound_N3.std_spectral_base();// Result to be stored there.
    
    Scalar Kappa  = bb*contract(contract(Kij,0,sss,0),0,sss,0) - contract(sss,0, lapse.derive_cov(gamij),0);
    Scalar Matter = 8.*M_PI*(bb*Ee);
    Matter.annule_domain(nz2-1);
    Matter += - 8.*M_PI*(bb + lapse)*contract(Jj,0,sss_down,0);
    Matter.annule_domain(nz2-1);
    Matter +=  8.*M_PI*lapse*contract(contract(Sij,0,sss_down,0),0,sss_down,0);
  
  
  // 2d interpolation of the Kappa constant and the matter terms.
  
  Scalar Kappa2 = sph.get_ricci();
  Kappa2.annule_hard();
  Scalar bb2 = Kappa2;
  Scalar Matter2 = Kappa2;
  Scalar nn2 = Kappa2;

  const Metric_flat& flat2 = Kappa2.get_mp().flat_met_spher();

  int nt = (*Kappa2.get_mp().get_mg()).get_nt(0);
  int np = (*Kappa2.get_mp().get_mg()).get_np(0);
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	  Kappa2.set_grid_point(0,k,j,0) = Kappa.val_grid_point(1,k,j,0);
	  bb2.set_grid_point(0,k,j,0) = bb.val_grid_point(1,k,j,0);
	  Matter2.set_grid_point(0,k,j,0) = Matter.val_grid_point(1,k,j,0);
	  nn2.set_grid_point(0,k,j,0) = lapse.val_grid_point(1,k,j,0);
      }
      
    Kappa2.std_spectral_base();
    bb2.std_spectral_base();
    Matter2.std_spectral_base();

    Scalar Tplus = sph.theta_plus();
    Scalar Tminus = sph.theta_minus();
    Scalar Ricci = sph.get_ricci();
    Vector Ll = sph.get_ll();
    Sym_tensor Shear = sph.shear();
    Metric qab = sph.get_qab();

    Scalar source = 2.*contract(Ll.up_down(qab),0,sph.derive_cov2d(bb2),0) + bb2*(contract(sph.derive_cov2d(Ll.up_down(qab)),0,1)+ contract(Ll,0,Ll.up_down(qab),0)- 0.5*(Ricci + Tplus*Tminus) + 0.25*Tplus*Tplus + 0.5*contract(Shear.up_down(qab),0,1,Shear,0,1)) + Matter2 + Tplus*Kappa2 + dttheta;
 
    Scalar lapb = contract(sph.derive_cov2d(sph.derive_cov2d(bb2)).up(0,qab),0,1);
    
    source = source + lapb;

  
    Scalar source_add =  2.*contract(Ll.up_down(qab),0,sph.derive_cov2d(nn2),0); 

    source = source - source_add; 
    
    Scalar sqrt_q_h2 = sph.sqrt_q() * sph.get_hsurf() * sph.get_hsurf() ;
    Tensor Delta2 = sph.delta();
    Sym_tensor qab_con = sph.get_qab().con() / (sph.get_hsurf() * sph.get_hsurf()) ;
  qab_con.set(1,1) = 1. ;
  qab_con.set(1,2) = 0. ;
  qab_con.set(1,3) = 0. ;
  qab_con.std_spectral_base() ;

  Sym_tensor hab =(qab_con*sqrt_q_h2 - flat2.con()) / (sph.get_hsurf()*sph.get_hsurf()) ;
  hab.set(1,1) = 1. ;
  hab.set(1,2) = 0. ;
  hab.set(1,3) = 0. ;
  hab.std_spectral_base() ;

  Vector dN = sph.derive_cov2dflat(nn2);
  Tensor ddN = sph.derive_cov2dflat(dN);

  Scalar lap_rem = contract(qab_con, 0,1, contract(Delta2,0,dN,0),0,1)*sqrt_q_h2 - sph.get_hsurf()*sph.get_hsurf()*contract(hab,0,1,ddN,0,1);// What remains of the laplacian

  source = sqrt_q_h2*source + lap_rem;

  Scalar bound_N=nn2;

    if (sph_sym == true){
      Scalar lapang_s_par = (contract(sph.derive_cov2d(Ll.up_down(qab)),0,1)+ contract(Ll,0,Ll.up_down(qab),0)- 0.5*(Ricci + Tplus*Tminus) - 0.25*Tplus*Tplus - 0.5*contract(Shear.up_down(qab),0,1,Shear,0,1))*sqrt_q_h2;
      double lapang_par = lapang_s_par.val_grid_point(0,0,0,0);
     
      bound_N = source.poisson_angu(lapang_par);
      bound_N.set_spectral_va().coef_i();
      bound_N.set_spectral_va().ylm();
    }

    else{
      cout << "non_spherical case has not been treated yet!" << endl;}

    // Interpolation to get a 3d value (as poisson solvers take this for boundary...)
  
  
 

  int nr2 = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt2 = (*lapse.get_mp().get_mg()).get_nt(1);
  int np2 = (*lapse.get_mp().get_mg()).get_np(1);


  for (int f= 0; f<nz2; f++)
    for (int k=0; k<np2; k++)
      for (int j=0; j<nt2; j++) {
	for (int l=0; l<nr2; l++) {
		
	  bound_N3.set_grid_point(f,k,j,l) = bound_N.val_grid_point(0,k,j,0);

	
	}
      }
  if (nz2 >2){
    bound_N3.annule_domain(nz2 - 1);

  }
    
  bound_N3.std_spectral_base();
  bound_N3.set_spectral_va().ylm();

  // End interpolation

    p_get_BC_lapse_3 = new Scalar(bound_N3); 
    
  }
  return *p_get_BC_lapse_3 ;
}





// Used to retrieve d (theta)/dt over only initial data. Probably obsolete in a practical use on CoCoNuT. Uses the same formalism as get_BC_lapse_3().
const Scalar& Excision_surf::derive_t_expa( Scalar& Ee, Vector& Jj, Sym_tensor& Sij)const{
if (p_derive_t_expa == 0x0){

  int nz2 = (*lapse.get_mp().get_mg()).get_nzone();
    // Radial vector for the full 3-metric.
     Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
     Scalar bb = contract (shift,0, sss_down,0); 
    
    Scalar Kappa  = bb*contract(contract(Kij,0,sss,0),0,sss,0) - contract(sss,0, lapse.derive_cov(gamij),0);
    Scalar Matter = 8.*M_PI*(bb*Ee);
    Matter.annule_domain(nz2-1);
    Matter += - 8.*M_PI*(bb + lapse)*contract(Jj,0,sss_down,0);
    Matter.annule_domain(nz2-1);
    Matter +=  8.*M_PI*lapse*contract(contract(Sij,0,sss_down,0),0,sss_down,0);
  
  
  // 2d interpolation of the Kappa constant and the matter terms.
  
  Scalar Kappa2 = sph.get_ricci();
  Kappa2.annule_hard();
  Scalar bb2 = Kappa2;
  Scalar Matter2 = Kappa2;
  Scalar nn2 = Kappa2;

  int nt = (*Kappa2.get_mp().get_mg()).get_nt(0);
  int np = (*Kappa2.get_mp().get_mg()).get_np(0);
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	  Kappa2.set_grid_point(0,k,j,0) = Kappa.val_grid_point(1,k,j,0);
	  bb2.set_grid_point(0,k,j,0) = bb.val_grid_point(1,k,j,0);
	  Matter2.set_grid_point(0,k,j,0) = Matter.val_grid_point(1,k,j,0);
	  nn2.set_grid_point(0,k,j,0) = lapse.val_grid_point(1,k,j,0);
      }
      
    Kappa2.std_spectral_base();
    bb2.std_spectral_base();
    Matter2.std_spectral_base();

    Scalar Tplus = sph.theta_plus();
    Scalar Tminus = sph.theta_minus();
    Scalar Ricci2 = sph.get_ricci();
    Vector Ll = sph.get_ll();
    Sym_tensor Shear = sph.shear();
    Metric qab = sph.get_qab();

    Scalar source = 2.*contract(Ll.up_down(qab),0,sph.derive_cov2d(nn2 -bb2),0) + (nn2- bb2)*(contract(sph.derive_cov2d(Ll.up_down(qab)),0,1)+ contract(Ll,0,Ll.up_down(qab),0)- 0.5*(Ricci2 + Tplus*Tminus)) -(bb2 +nn2)*(0.25*Tplus*Tplus + 0.5*contract(Shear.up_down(qab),0,1,Shear,0,1)) - Matter2 - Tplus*Kappa2 ;
 
    Scalar lapNmb = contract(sph.derive_cov2d(sph.derive_cov2d(nn2-bb2)).up(0,qab),0,1);
    
    source = source + lapNmb;
    source.set_spectral_va().ylm();

    Scalar difftheta = source*delta_t;
    cout << "mean difference between old and new expansion" << endl;
    cout << difftheta.val_grid_point(0,0,0,0) << endl;

 
  p_derive_t_expa = new Scalar (source);
 }

 return *p_derive_t_expa ;
}








// Source for the Dirichlet BC on the shift
const Vector& Excision_surf::get_BC_shift_1(double Omega) const{
  if (p_get_BC_shift_1 == 0x0){
    
    // Radial vector for the full 3-metric.
    Vector sss = gamij.radial_vect();
    
    // Boundary value for the radial part of the shift
    Scalar bound = lapse ;
    //   bound = bound + 0.05; // REMOVE real quick

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
const Vector& Excision_surf::get_BC_shift_2(double c_bb_lap, double c_bb_fin, double c_V_lap, double epsilon) const{
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
  Sym_tensor q_downdown = gamij.cov() - sss_down*sss_down;
 Tensor q_updown = q_upup.down(1, gamij); 

 

 

  // Calculation of the conformal 2d laplacian of V
 Tensor d_V = contract(q_updown, 1, contract(q_updown,0 , V_par.derive_cov(gamij),1 ),1);
 Tensor d_V_con = contract(d_V,1,q_upup,1);
  Tensor dd_V = d_V_con.derive_cov(gamij);
  // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
  dd_V = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V, 2), 2), 2);
  Tensor dd_Vdown = contract (dd_V,1,q_downdown,1);
  Vector Ricci1 = contract(dd_Vdown,0,1) - contract(dd_Vdown,0,2);    
  Vector Ricci2 = contract (q_upup,1,Ricci1,0);
  Vector lap_V = contract(dd_V, 1,2);



//   Tensor dd_V = V_par.derive_con(gamij);
//   dd_V = contract(q_updown, 1, contract(q_updown,1 ,dd_V, 0), 1);
//   Vector lap_V = contract(q_updown, 1, contract(dd_V.derive_cov(gamij),1,2), 0);


  
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
    //     ricci3.annule_domain(0);
    ricci3.annule_domain(nz - 1);

  }
    // End Mapping interpolation
  
    // Construction of the Ricci COV tensor on the sphere
 
  Sym_tensor ricci_t = gamij.cov() - sss_down*sss_down;
  ricci_t = 0.5*ricci3*ricci_t;
  ricci_t.std_spectral_base();
 
  Tensor ricci_t_updown = contract(q_upup,0, ricci_t,0); 
  
  // Calculation of ff 

  // Vector ffV = c_V_lap*lapse*(lap_V + contract(ricci_t_updown,1, V_par,0));
   Vector ffV = c_V_lap*lapse*(lap_V + Ricci2);

  ffV.annule_domain(nz-1);


//   cout << "verification of vanishing shear" << endl;
//   cout << "points on the surface" << endl;
  
//   Tbl val_ffV(3, nt, np);
//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();
 
 

   ffV = ffV + epsilon*lapse*contract(V_par,0,V_par.derive_cov(gamij),1); // Add of dragging term for time transport.
  ffV.annule_domain(nz-1);

  ffV.std_spectral_base();

//   cout << "verification of vanishing shear with dragging" << endl;
//  cout << "points on the surface" << endl;
//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();
 
 


  // Definition of k_1
  Vector k_1V =delta_t*ffV;
   
  // Intermediate value of Npsi, for Runge-Kutta 2nd order scheme
  if (nz >2){
    k_1V.annule_domain(nz-1);
  }                             // Patch to avoid dzpuis problems if existent.
  Vector V_par_int = V_par + k_1V;// V_par_int.std_spectral_base();

  // Calculation of the conformal 2d laplacian of V
 Tensor d_V_int = contract(q_updown, 1, contract(q_updown,0 , V_par_int.derive_cov(gamij),1 ),1);
 Tensor d_V_con_int = contract(d_V_int,1,q_upup,1);
  Tensor dd_V_int = d_V_con_int.derive_cov(gamij);
  // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
  dd_V_int = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V_int, 2), 2), 2);
  Tensor dd_Vdown_int = contract (dd_V_int,1,q_downdown,1);
  Vector Ricci1_int = contract(dd_Vdown_int,0,1) - contract(dd_Vdown_int,0,2);    
  Vector Ricci2_int = contract (q_upup,1,Ricci1_int,0);
  Vector lap_V_int = contract(dd_V_int, 1,2);


  // Recalculation of ff with intermediate values.
//   Tensor dd_V_int = V_par_int.derive_con(gamij).derive_cov(gamij);
//   // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
//   dd_V_int = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V_int, 2), 2), 2);
//   Vector lap_V_int = contract(dd_V_int, 1,2);
 
  Vector ffV_int =  c_V_lap*lapse*(lap_V_int + contract(ricci_t_updown,1, V_par_int,0));
  ffV_int.annule_domain(nz-1);
   ffV_int = ffV_int + epsilon*lapse*contract(V_par_int,0,V_par_int.derive_cov(gamij),1); // Add of dragging term for time transport.

   //  ffV_int = -ffV_int; // Only a test..


  ffV_int.annule_domain(nz-1);
 

//   cout << "verification of vanishing shear for intermediate RK value" << endl;
//   cout << "points on the surface" << endl;


//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();



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





// Source for the Dirichlet BC on the shift, based on kinematical relation for the radial part, and a parabolic evolution for the tangent part.
// It takes as a fixed argument the time derivative of the psi factor, that can be deduced from the function.get_BC_conf_fact_4(). 
const Vector& Excision_surf::get_BC_shift_3(Scalar& dtpsi, double c_V_lap, double epsilon) const{
  if (p_get_BC_shift_3 == 0x0){

    // Radial vector for the full 3-metric.
     Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
     
    const Metric_flat& flat = lapse.get_mp().flat_met_spher() ; 
    
 int nz = (*lapse.get_mp().get_mg()).get_nzone();

//     // Boundary value for the radial part of the shift: parabolic driver for (b-N)
     //  Scalar bound = lapse ; 
     Scalar bb = contract (shift,0, sss_down,0); 

     Scalar pure_source = contract(sss,0,bb.derive_cov(flat),0)*conf_fact*1./6.;
     pure_source.annule_domain(nz-1); // CAREFUL
     pure_source = dtpsi - pure_source;
     pure_source.annule_domain(nz-1);
     
     Scalar factor = conf_fact*sss.divergence(flat)*1./6.;
     factor.annule_domain(nz-1);
     factor = factor + contract (sss, 0, conf_fact.derive_cov(flat),0);
    
     Scalar source = pure_source/factor;
     


     Scalar bb2 = source; // Look out for additional term for time variation of lapse and sss
  

  // Tangent part of the shift, with parabolic driver
  
  
  Vector V_par = shift - bb*sss;
  Sym_tensor q_upup = gamij.con() - sss*sss;
  Sym_tensor q_downdown = gamij.cov() - sss_down*sss_down;
 Tensor q_updown = q_upup.down(1, gamij); 

 
  // Calculation of the conformal 2d laplacian of V
 Tensor d_V = contract(q_updown, 1, contract(q_updown,0 , V_par.derive_cov(gamij),1 ),1);
 Tensor d_V_con = contract(d_V,1,q_upup,1);
  Tensor dd_V = d_V_con.derive_cov(gamij);
  // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
  dd_V = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V, 2), 2), 2);
  Tensor dd_Vdown = contract (dd_V,1,q_downdown,1);
  Vector Ricci1 = contract(dd_Vdown,0,1) - contract(dd_Vdown,0,2);    
  Vector Ricci2 = contract (q_upup,1,Ricci1,0);
  Vector lap_V = contract(dd_V, 1,2);



//   Tensor dd_V = V_par.derive_con(gamij);
//   dd_V = contract(q_updown, 1, contract(q_updown,1 ,dd_V, 0), 1);
//   Vector lap_V = contract(q_updown, 1, contract(dd_V.derive_cov(gamij),1,2), 0);


  
  // 3d interpolation of the Ricci scalar on the surface.
  
  Scalar ricci2 = sph.get_ricci();
  
     // Start Mapping interpolation
 
      Scalar ricci3 (lapse.get_mp()); 
 
  ricci3.allocate_all();
  ricci3.std_spectral_base();


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
    //     ricci3.annule_domain(0);
    ricci3.annule_domain(nz - 1);

  }
    // End Mapping interpolation
  
    // Construction of the Ricci COV tensor on the sphere
 
  Sym_tensor ricci_t = gamij.cov() - sss_down*sss_down;
  ricci_t = 0.5*ricci3*ricci_t;
  ricci_t.std_spectral_base();
 
  Tensor ricci_t_updown = contract(q_upup,0, ricci_t,0); 
  
  // Calculation of ff 

  // Vector ffV = c_V_lap*lapse*(lap_V + contract(ricci_t_updown,1, V_par,0));
   Vector ffV = c_V_lap*lapse*(lap_V + Ricci2);

  ffV.annule_domain(nz-1);


//   cout << "verification of vanishing shear" << endl;
//   cout << "points on the surface" << endl;
  
//   Tbl val_ffV(3, nt, np);
//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();
 
 

   ffV = ffV + epsilon*lapse*contract(V_par,0,V_par.derive_cov(gamij),1); // Add of dragging term for time transport.
  ffV.annule_domain(nz-1);

  ffV.std_spectral_base();

//   cout << "verification of vanishing shear with dragging" << endl;
//  cout << "points on the surface" << endl;
//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();
 
 


  // Definition of k_1
  Vector k_1V =delta_t*ffV;
   
  // Intermediate value of Npsi, for Runge-Kutta 2nd order scheme
  if (nz >2){
    k_1V.annule_domain(nz-1);
  }                             // Patch to avoid dzpuis problems if existent.
  Vector V_par_int = V_par + k_1V;// V_par_int.std_spectral_base();

  // Calculation of the conformal 2d laplacian of V
 Tensor d_V_int = contract(q_updown, 1, contract(q_updown,0 , V_par_int.derive_cov(gamij),1 ),1);
 Tensor d_V_con_int = contract(d_V_int,1,q_upup,1);
  Tensor dd_V_int = d_V_con_int.derive_cov(gamij);
  // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
  dd_V_int = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V_int, 2), 2), 2);
  Tensor dd_Vdown_int = contract (dd_V_int,1,q_downdown,1);
  Vector Ricci1_int = contract(dd_Vdown_int,0,1) - contract(dd_Vdown_int,0,2);    
  Vector Ricci2_int = contract (q_upup,1,Ricci1_int,0);
  Vector lap_V_int = contract(dd_V_int, 1,2);


  // Recalculation of ff with intermediate values.
//   Tensor dd_V_int = V_par_int.derive_con(gamij).derive_cov(gamij);
//   // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
//   dd_V_int = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V_int, 2), 2), 2);
//   Vector lap_V_int = contract(dd_V_int, 1,2);
 
  Vector ffV_int =  c_V_lap*lapse*(lap_V_int + contract(ricci_t_updown,1, V_par_int,0));
  ffV_int.annule_domain(nz-1);
   ffV_int = ffV_int + epsilon*lapse*contract(V_par_int,0,V_par_int.derive_cov(gamij),1); // Add of dragging term for time transport.

   //  ffV_int = -ffV_int; // Only a test..


  ffV_int.annule_domain(nz-1);
 

//   cout << "verification of vanishing shear for intermediate RK value" << endl;
//   cout << "points on the surface" << endl;


//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();



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
       p_get_BC_shift_3 = new Vector(bound_shift);
}
  return *p_get_BC_shift_3 ;
}










// Source for the Dirichlet BC on the shift, based on projection of Einstein Equations for the radial part, and a parabolic evolution for the tangent part.
// It takes as fixed arguments the time derivative of the expansion, the matter terms and a boolean specifying whether or not we are in spherical symmetry.. 
const Vector& Excision_surf::get_BC_shift_4(Scalar& dttheta, Scalar& Ee, Vector& Jj, Sym_tensor& Sij, double c_V_lap, double epsilon, bool sph_sym) const{
  if (p_get_BC_shift_4 == 0x0){

    // Radial vector for the full 3-metric.
     Vector sss = gamij.radial_vect();
     Vector sss_down = sss.up_down(gamij);
         
     int nz = (*lapse.get_mp().get_mg()).get_nzone();

//     // Boundary value for the radial part of the shift: parabolic driver for (b-N)
     //  Scalar bound = lapse ; 
     Scalar bb = contract (shift,0, sss_down,0); 



     Scalar bound_b3(lapse.get_mp()); bound_b3.allocate_all(); bound_b3.set_spectral_base(bb.get_spectral_base());// Result to be stored there.
    
    Scalar Kappa  = bb*contract(contract(Kij,0,sss,0),0,sss,0) - contract(sss,0, lapse.derive_cov(gamij),0);
    Scalar Matter = 8.*M_PI*(bb*Ee);
    Matter.annule_domain(nz-1);
    Matter += - 8.*M_PI*(bb + lapse)*contract(Jj,0,sss_down,0);
    Matter.annule_domain(nz-1);
    Matter +=  8.*M_PI*lapse*contract(contract(Sij,0,sss_down,0),0,sss_down,0);
  
  
  // 2d interpolation of the Kappa constant and the matter terms.
  
  Scalar Kappa2 = sph.get_ricci();
  Kappa2.annule_hard();
  Scalar bb2 = Kappa2;
  Scalar Matter2 = Kappa2;
  Scalar nn2 = Kappa2;

  const Metric_flat& flat2 = Kappa2.get_mp().flat_met_spher();

  int nt = (*Kappa2.get_mp().get_mg()).get_nt(0);
  int np = (*Kappa2.get_mp().get_mg()).get_np(0);
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	  Kappa2.set_grid_point(0,k,j,0) = Kappa.val_grid_point(1,k,j,0);
	  bb2.set_grid_point(0,k,j,0) = bb.val_grid_point(1,k,j,0);
	  Matter2.set_grid_point(0,k,j,0) = Matter.val_grid_point(1,k,j,0);
	  nn2.set_grid_point(0,k,j,0) = lapse.val_grid_point(1,k,j,0);
      }
      
    Kappa2.std_spectral_base();
    bb2.std_spectral_base();
    Matter2.std_spectral_base();
    nn2.std_spectral_base();


    Scalar Tplus = sph.theta_plus();
    Scalar Tminus = sph.theta_minus();
    Scalar Ricci = sph.get_ricci();
    Vector Ll = sph.get_ll();
    Sym_tensor Shear = sph.shear();
    Metric qab = sph.get_qab();

    Scalar source = 2.*contract(Ll.up_down(qab),0,sph.derive_cov2d(nn2),0) + nn2*(contract(sph.derive_cov2d(Ll.up_down(qab)),0,1)+ contract(Ll,0,Ll.up_down(qab),0)- 0.5*(Ricci + Tplus*Tminus) - 0.25*Tplus*Tplus - 0.5*contract(Shear.up_down(qab),0,1,Shear,0,1)) - Matter2 - Tplus*Kappa2 - dttheta;
 
    Scalar lapN = contract(sph.derive_cov2d(sph.derive_cov2d(nn2)).up(0,qab),0,1);
    
    source = source + lapN;

  
    Scalar source_add =  2.*contract(Ll.up_down(qab),0,sph.derive_cov2d(bb2),0); 

    source = source - source_add; 
    
    Scalar sqrt_q_h2 = sph.sqrt_q() * sph.get_hsurf() * sph.get_hsurf() ;
    Tensor Delta2 = sph.delta();
    Sym_tensor qab_con = sph.get_qab().con() / (sph.get_hsurf() * sph.get_hsurf()) ;
  qab_con.set(1,1) = 1. ;
  qab_con.set(1,2) = 0. ;
  qab_con.set(1,3) = 0. ;
  qab_con.std_spectral_base() ;

  Sym_tensor hab =(qab_con*sqrt_q_h2 - flat2.con()) / (sph.get_hsurf()*sph.get_hsurf()) ;
  hab.set(1,1) = 1. ;
  hab.set(1,2) = 0. ;
  hab.set(1,3) = 0. ;
  hab.std_spectral_base() ;

  Vector db = sph.derive_cov2dflat(bb2);
  Tensor ddb = sph.derive_cov2dflat(db);

  Scalar lap_rem = contract(qab_con, 0,1, contract(Delta2,0,db,0),0,1)*sqrt_q_h2 - sph.get_hsurf()*sph.get_hsurf()*contract(hab,0,1,ddb,0,1);// What remains of the laplacian

  source = sqrt_q_h2*source + lap_rem;

  Scalar bound_b=bb2;
      Scalar lapang_s_par = (contract(sph.derive_cov2d(Ll.up_down(qab)),0,1)+ contract(Ll,0,Ll.up_down(qab),0)- 0.5*(Ricci + Tplus*Tminus) + 0.25*Tplus*Tplus + 0.5*contract(Shear.up_down(qab),0,1,Shear,0,1))*sqrt_q_h2;
      double lapang_par = lapang_s_par.val_grid_point(0,0,0,0);

    if (sph_sym == true){
      bound_b = source/lapang_par;
      bound_b.set_spectral_va().coef_i();
      bound_b.set_spectral_va().ylm();
    }

    else{
      cout << "non spherical case has not been treated yet!" << endl;  
    }

    // Interpolation to get a 3d value (as poisson solvers take this for boundary...)
  
  
 

  int nr2 = (*lapse.get_mp().get_mg()).get_nr(1);
  int nt2 = (*lapse.get_mp().get_mg()).get_nt(1);
  int np2 = (*lapse.get_mp().get_mg()).get_np(1);


  for (int f= 0; f<nz; f++)
    for (int k=0; k<np2; k++)
      for (int j=0; j<nt2; j++) {
	for (int l=0; l<nr2; l++) {
		
	  bound_b3.set_grid_point(f,k,j,l) = bound_b.val_grid_point(0,k,j,0);

	
	}
      }
  if (nz >2){
    bound_b3.annule_domain(nz - 1);

  }
    
  bound_b3.std_spectral_base();
  bound_b3.set_spectral_va().ylm();




     Scalar bbb2 = bound_b3; // Look out for additional term for time variation of lapse and sss
  

  // Tangent part of the shift, with parabolic driver
  
  
  Vector V_par = shift - bb*sss;
  Sym_tensor q_upup = gamij.con() - sss*sss;
  Sym_tensor q_downdown = gamij.cov() - sss_down*sss_down;
 Tensor q_updown = q_upup.down(1, gamij); 

 
  // Calculation of the conformal 2d laplacian of V
 Tensor d_V = contract(q_updown, 1, contract(q_updown,0 , V_par.derive_cov(gamij),1 ),1);
 Tensor d_V_con = contract(d_V,1,q_upup,1);
  Tensor dd_V = d_V_con.derive_cov(gamij);
  // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
  dd_V = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V, 2), 2), 2);
  Tensor dd_Vdown = contract (dd_V,1,q_downdown,1);
  Vector Ricci1 = contract(dd_Vdown,0,1) - contract(dd_Vdown,0,2);    
  Vector Ricci2 = contract (q_upup,1,Ricci1,0);
  Vector lap_V = contract(dd_V, 1,2);



//   Tensor dd_V = V_par.derive_con(gamij);
//   dd_V = contract(q_updown, 1, contract(q_updown,1 ,dd_V, 0), 1);
//   Vector lap_V = contract(q_updown, 1, contract(dd_V.derive_cov(gamij),1,2), 0);


  
  // 3d interpolation of the Ricci scalar on the surface.
  
  Scalar ricci2 = sph.get_ricci();
  
     // Start Mapping interpolation
 
      Scalar ricci3 (lapse.get_mp()); 
 
  ricci3.allocate_all();
  ricci3.std_spectral_base();

  for (int f= 0; f<nz; f++)
    for (int k=0; k<np2; k++)
      for (int j=0; j<nt2; j++) {
	for (int l=0; l<nr2; l++) {
		
	  ricci3.set_grid_point(f,k,j,l) = ricci2.val_grid_point(0,k,j,0);

	
	}
      }
  if (nz >2){
    //     ricci3.annule_domain(0);
    ricci3.annule_domain(nz - 1);

  }
    // End Mapping interpolation
  
    // Construction of the Ricci COV tensor on the sphere
 
  Sym_tensor ricci_t = gamij.cov() - sss_down*sss_down;
  ricci_t = 0.5*ricci3*ricci_t;
  ricci_t.std_spectral_base();
 
  Tensor ricci_t_updown = contract(q_upup,0, ricci_t,0); 
  
  // Calculation of ff 

  // Vector ffV = c_V_lap*lapse*(lap_V + contract(ricci_t_updown,1, V_par,0));
   Vector ffV = c_V_lap*lapse*(lap_V + Ricci2);

  ffV.annule_domain(nz-1);


//   cout << "verification of vanishing shear" << endl;
//   cout << "points on the surface" << endl;
  
//   Tbl val_ffV(3, nt, np);
//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();
 
 

   ffV = ffV + epsilon*lapse*contract(V_par,0,V_par.derive_cov(gamij),1); // Add of dragging term for time transport.
  ffV.annule_domain(nz-1);

  ffV.std_spectral_base();

//   cout << "verification of vanishing shear with dragging" << endl;
//  cout << "points on the surface" << endl;
//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();
 
 


  // Definition of k_1
  Vector k_1V =delta_t*ffV;
   
  // Intermediate value of Npsi, for Runge-Kutta 2nd order scheme
  if (nz >2){
    k_1V.annule_domain(nz-1);
  }                             // Patch to avoid dzpuis problems if existent.
  Vector V_par_int = V_par + k_1V;// V_par_int.std_spectral_base();

  // Calculation of the conformal 2d laplacian of V
 Tensor d_V_int = contract(q_updown, 1, contract(q_updown,0 , V_par_int.derive_cov(gamij),1 ),1);
 Tensor d_V_con_int = contract(d_V_int,1,q_upup,1);
  Tensor dd_V_int = d_V_con_int.derive_cov(gamij);
  // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
  dd_V_int = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V_int, 2), 2), 2);
  Tensor dd_Vdown_int = contract (dd_V_int,1,q_downdown,1);
  Vector Ricci1_int = contract(dd_Vdown_int,0,1) - contract(dd_Vdown_int,0,2);    
  Vector Ricci2_int = contract (q_upup,1,Ricci1_int,0);
  Vector lap_V_int = contract(dd_V_int, 1,2);


  // Recalculation of ff with intermediate values.
//   Tensor dd_V_int = V_par_int.derive_con(gamij).derive_cov(gamij);
//   // Vector lap_V = contract(dd_V.derive_cov(gamij),1,2);
//   dd_V_int = contract(q_updown, 1, contract(q_updown,1 ,contract(q_updown,0,dd_V_int, 2), 2), 2);
//   Vector lap_V_int = contract(dd_V_int, 1,2);
 
  Vector ffV_int =  c_V_lap*lapse*(lap_V_int + contract(ricci_t_updown,1, V_par_int,0));
  ffV_int.annule_domain(nz-1);
   ffV_int = ffV_int + epsilon*lapse*contract(V_par_int,0,V_par_int.derive_cov(gamij),1); // Add of dragging term for time transport.

   //  ffV_int = -ffV_int; // Only a test..


  ffV_int.annule_domain(nz-1);
 

//   cout << "verification of vanishing shear for intermediate RK value" << endl;
//   cout << "points on the surface" << endl;


//   val_ffV.set_etat_qcq();
//   for(int mm=0; mm <np; mm++)
//     for (int pp=0; pp<nt; pp++){
//       val_ffV.set(0,pp,mm) = ffV(1).val_grid_point(1,mm,pp,0);
//       val_ffV.set(1,pp,mm) = ffV(2).val_grid_point(1,mm,pp,0);
//       val_ffV.set(2,pp,mm) = ffV(3).val_grid_point(1,mm,pp,0);
//     }

//   //  cout << val_ffV<< endl;
//   cout << max(val_ffV) << endl;
//   val_ffV.set_etat_nondef();



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
       Vector bound_shift = bbb2*sss + bound_V;
       bound_shift.std_spectral_base();
       p_get_BC_shift_4 = new Vector(bound_shift);
}
  return *p_get_BC_shift_4 ;
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
  ff = ff -lapse*c_npsi_fin*npsi_fin;
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

// Source for the Dirichlet BC on (N*Psi1), with a Kerr_Schild-like form for the lapse.
const Scalar& Excision_surf::get_BC_Npsi_3(double n_0, double beta) const{
  if (p_get_BC_Npsi_3 == 0x0){


    const Coord& ct = lapse.get_mp().cost;
    Scalar boundN(lapse.get_mp()); boundN = sqrt(n_0 + beta*ct*ct); // Kerr_Schild form for the lapse.
    boundN.std_spectral_base();
    Scalar bound_npsi = boundN*conf_fact;
    bound_npsi.set_spectral_va().ylm();

//     Scalar npsi = lapse*conf_fact; npsi.std_spectral_base();
//   Scalar ff = lapse*(c_npsi_lap*npsi.lapang() + c_npsi_fin*npsi);
//   ff.std_spectral_base();
//   ff += -lapse*c_npsi_fin*npsi_fin;
//   ff.std_spectral_base();


//   // Definition of k_1
//   Scalar k_1 =delta_t*ff;
   
//   // Intermediate value of Npsi, for Runge-Kutta 2nd order scheme
//   Scalar npsi_int = npsi + k_1; npsi_int.std_spectral_base();

//   // Recalculation of ff with intermediate values. 
//   Scalar ff_int =  lapse*(c_npsi_lap*npsi_int.lapang() + c_npsi_fin*(npsi_int - npsi_fin));
 
//   // Definition of k_2
//   Scalar k_2 = delta_t*ff_int; 
//   k_2.std_spectral_base();
 
//   // Result of RK2 evolution
//   Scalar bound_npsi = npsi + k_2;
//   bound_npsi.std_spectral_base();
//   bound_npsi.set_spectral_va().ylm();



//     Scalar bound_Npsi = value*conf_fact;
//     bound_Npsi.set_spectral_va().ylm();
     p_get_BC_Npsi_3 = new Scalar(bound_npsi);
    
}
  return *p_get_BC_Npsi_3 ;
}


// Source for a Dirichlet BC on (N*Psi1), fixing the surface gravity as a constant in space and time.
const Scalar& Excision_surf::get_BC_Npsi_4(double Kappa) const{
  if (p_get_BC_Npsi_4 == 0x0){


    const Vector s_i = gamij.radial_vect();
    Scalar bound_npsi = contract(s_i,0, lapse.derive_cov(gamij),0); bound_npsi.set_dzpuis(0);
    bound_npsi = bound_npsi - Kappa;
    bound_npsi.std_spectral_base();
    bound_npsi = bound_npsi*(conf_fact/contract(contract(Kij,0,s_i,0),0,s_i,0));
    bound_npsi.set_dzpuis(0); bound_npsi.std_spectral_base();
    bound_npsi.set_spectral_va().ylm();

     p_get_BC_Npsi_4 = new Scalar(bound_npsi);
    
}
  return *p_get_BC_Npsi_4 ;
}



// Source for a Neumann BC on (N*Psi1), fixing the surface gravity as a constant in space and time. // Check redundancy with the previous one.
const Scalar& Excision_surf::get_BC_Npsi_5(double Kappa) const{
  if (p_get_BC_Npsi_5 == 0x0){


    const Vector s_i = gamij.radial_vect();
  int nz = (*lapse.get_mp().get_mg()).get_nzone();
    Scalar bound_npsi = lapse*contract(s_i,0, conf_fact.derive_cov(gamij),0); bound_npsi.annule_domain(nz -1);
    Scalar bound_npsi2 = pow(conf_fact,3)*lapse*contract(contract(Kij,0,s_i,0),0,s_i,0);
    bound_npsi2.annule_domain(nz-1);
    bound_npsi += bound_npsi2;
    bound_npsi = bound_npsi + pow(conf_fact,3)*(Kappa); bound_npsi.annule_domain(nz -1);
    bound_npsi = bound_npsi -(conf_fact*conf_fact*contract(s_i,0, (conf_fact*lapse).derive_cov(gamij),0) - (conf_fact*lapse).dsdr());
    //    bound_npsi = bound_npsi*(conf_fact/contract(contract(Kij,0,s_i,0),0,s_i,0));
    bound_npsi.annule_domain(nz-1); bound_npsi.std_spectral_base();
    bound_npsi.set_spectral_va().ylm();

     p_get_BC_Npsi_5 = new Scalar(bound_npsi);
    
}
  return *p_get_BC_Npsi_5 ;
}





 void Excision_surf::sauve(FILE* ) const {

   cout << "c'est pas fait!" << endl ;
   return ; 

 }
