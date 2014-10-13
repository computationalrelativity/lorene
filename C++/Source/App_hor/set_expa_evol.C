
// Header Lorene:
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "math.h"
#include "metric.h"
#include "param.h"
#include "param_elliptic.h"
#include "tensor.h"
#include "unites.h"
#include "excision_surf.h"




namespace Lorene {
//gives a new value for expansion rescaled with lapse (and its time derivative) obtained by parabolic evolution.
//All manipulated quantities are 2-dimensional.
void Excision_surf::set_expa_parab(double c_theta_lap, double c_theta_fin, Scalar& expa_fin){

    // Definition of ff
    // ================

    // Start Mapping interpolation
   
  if (expa.get_spectral_va().get_etat() == 0)
    {
      //       Scalar theta = sph.theta_plus();
      Scalar theta (lapse.get_mp()); theta = 3.; theta.std_spectral_base();
 
       set_expa() = theta;}


  
  Scalar thetaplus = expa;
   
  // Interpolation for the lapse (to get a 2d quantity);
    Scalar lapse2(thetaplus.get_mp());
    lapse2.annule_hard();
    lapse2.std_spectral_base();

  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {

		
	  lapse2.set_grid_point(0,k,j,0) = lapse.val_grid_point(1,k,j,0);
	
	}

  lapse2.std_spectral_base();

 
    // End Mapping interpolation
//   cout << "convergence?" << endl;
//   cout << expa_fin3.val_grid_point(1,0,0,0) << endl;
//   cout << theta_plus3.val_grid_point(1,0,0,0) << endl;
  

  Scalar ff = lapse2*(c_theta_lap*thetaplus.lapang() + c_theta_fin*(thetaplus - expa_fin));
 

  // Definition of k_1
  Scalar k_1 =delta_t*ff;
   
  // Intermediate value of the expansion, for Runge-Kutta 2nd order scheme
  Scalar theta_int = thetaplus + k_1; theta_int.std_spectral_base();

  // Recalculation of ff with intermediate values. 
  Scalar ff_int =  lapse2*(c_theta_lap*theta_int.lapang() + c_theta_fin*(theta_int - expa_fin));
  ff_int.set_spectral_va().ylm();
  
  // Definition of k_2
  Scalar k_2 = delta_t*ff_int; 
 
  // Result of RK2 evolution
  Scalar bound_theta = thetaplus + k_2;
 

  // Assigns new values to the members expa() and dt_expa(). 
  set_expa()=bound_theta;
  set_dt_expa()= ff_int;
 
  return;
}



//gives a new value for the expansion (rescaled with lapse) and its time derivative, obtained by hyperbolic evolution.
// Parameters for the hyperbolic driver are determined by the function Excision_surf::get_evol_params_from_ID() 
// so that the expansion stays of regularity $C^{1}$ throughout.
// The scheme used is a RK4
// Final value is here necessarily zero for the expansion
//All manipulated quantities are 2-dimensional.
// Warning: no rescaling of time dimensionality by the lapse yet.

void Excision_surf::set_expa_hyperb(double alpha0, double beta0, double gamma0) {

    // Definition of ff
    // ================

    // Start Mapping interpolation
    Scalar thetaplus = expa;
    thetaplus.set_spectral_va().ylm();
    assert (dt_expa.get_spectral_va().get_etat() != 0);
    Scalar d_thetaplus = dt_expa;
    d_thetaplus.set_spectral_va().ylm();

    //////////////////////////////////////////////////////////////////////:
    // Interpolating for the lapse into the 2-dimensional surface (if necessary...)
    Scalar lapse2(thetaplus.get_mp());
    lapse2.annule_hard();
    lapse2.std_spectral_base();
 
  int nt = (*lapse.get_mp().get_mg()).get_nt(1);
  int np = (*lapse.get_mp().get_mg()).get_np(1);
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	  lapse2.set_grid_point(0,k,j,0) = lapse.val_grid_point(1,k,j,0);
	}

  lapse2.std_spectral_base();

  /////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////:
  // Starting the RK4 1step evolution for the second-order 2d PDE in spherical harmonics.
  //////////////////////////////////////////////////////////////////////////////////


  //Successive derivative estimates for the scheme

  //Step1
  Scalar k_1 = d_thetaplus; 
  Scalar K_1 = beta0*d_thetaplus + alpha0*d_thetaplus.lapang() + gamma0*thetaplus ;

  thetaplus = expa + 0.5*delta_t*k_1;
  d_thetaplus = dt_expa + 0.5*delta_t*K_1;

  //Step2
  Scalar k_2 =  d_thetaplus;
  Scalar K_2 = beta0*d_thetaplus + alpha0*d_thetaplus.lapang() + gamma0*thetaplus;
  
  thetaplus = expa + 0.5*delta_t*k_2;
  d_thetaplus = dt_expa + 0.5*delta_t*K_2;
 
  //Step3
  Scalar k_3 = d_thetaplus;
  Scalar K_3 = beta0*d_thetaplus + alpha0*d_thetaplus.lapang() + gamma0*thetaplus;
 
  thetaplus = expa + delta_t*k_3;
  d_thetaplus = dt_expa + delta_t*K_3;

  //Step4
  Scalar k_4 = d_thetaplus;
  Scalar K_4 = beta0*d_thetaplus + alpha0*d_thetaplus.lapang() + gamma0*thetaplus;


  // Result of RK2 evolution: assignment at evolved time.
  thetaplus = expa + (1./6.)*delta_t*(k_1 + 2.*k_2 + 2.*k_3 + k_4);
  d_thetaplus = dt_expa + (1./6.)*delta_t*(K_1 + 2.*K_2 + 2.*K_3 + K_4);
   

  set_expa() = thetaplus;
  set_dt_expa() = d_thetaplus;
 
  

  return;

}

}
