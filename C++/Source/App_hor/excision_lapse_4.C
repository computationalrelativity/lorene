
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
//Sym_tensor secmembre_kerr ( const Sym_tensor& hij, const Sym_tensor& aa,const Scalar& nn,const Scalar& ppsi,const Vector& bb) { 
const Scalar& Excision_surf::get_BC_lapse_4 (Scalar& old_nn, Vector& beta_point, Sym_tensor& strain_tens) const{ 

  using namespace Unites;
  if (p_get_BC_lapse_4 == 0x0){

    Scalar nn = lapse;
    Scalar ppsi = conf_fact;
    Vector bb = shift;
    Sym_tensor aa = Kij.up_down(gamij)*ppsi*ppsi*ppsi*ppsi;
    
  const int nz = (*(aa.get_mp().get_mg())).get_nzone(); 

  Sym_tensor hij = gamij.con();
   for (int ii=1; ii<=3; ii++)
    for(int jj=1; jj<=3; jj++)
      { hij.set(ii,jj).annule_hard();
      }
  hij.annule_domain(nz-1); 
  hij.std_spectral_base();    // Set non conformally flat part to zero.


  const Vector& beta = bb;
 
  const Scalar& psi4 = ppsi*ppsi*ppsi*ppsi;
  Scalar ln_psi = log(ppsi); ln_psi.std_spectral_base();
  ln_psi.annule_domain(nz-1);

  const Scalar qq = nn*ppsi*ppsi; 
   
  
  const Metric_flat& ff = (hij.get_mp()).flat_met_spher() ;
   
  const Sym_tensor& tgam_uu = ff.con(); 
    
  const Base_vect_spher& otriad = hij.get_mp().get_bvect_spher();

  Scalar nn_point(hij.get_mp());
  nn_point.annule_hard();
  nn_point.annule_domain(nz -1);
 
  nn_point.std_spectral_base(); 
 

			               
  const Sym_tensor& tgam_dd = ff.cov() ;    // {\tilde \gamma}_{ij}
  const Vector& dln_psi = ln_psi.derive_cov(ff) ; // D_i ln(Psi)
  const Vector& tdln_psi_u = ln_psi.derive_con(ff) ; // tD^i ln(Psi)
  const Vector& tdnn_u = nn.derive_con(ff) ;       // tD^i N
  const Vector& dqq = qq.derive_cov(ff) ;         // D_i Q
  const Scalar& div_beta = beta.divergence(ff) ;  // D_k beta^k
  Sym_tensor l_beta = beta.ope_killing_conf(ff) ; // Attention aux headers a inclure

  //==================================
  // Source for hij
  //==================================


  Scalar tmp(hij.get_mp()) ;
  Sym_tensor sym_tmp(hij.get_mp(), CON, otriad) ; 
       
  // Full quadratic part of source for h : S^{ij}
  // --------------------------------------------
        
  Sym_tensor ss(hij.get_mp(), CON, otriad) ;  // Source secondaire 
        
  sym_tmp = nn * (8.* tdln_psi_u * tdln_psi_u)
    + 4.*( tdln_psi_u * tdnn_u + tdnn_u * tdln_psi_u ) 
    - 0.3333333333333333 * 
    ( nn * ( 8.* contract(dln_psi, 0, tdln_psi_u, 0) )
      + 8.* contract(dln_psi, 0, tdnn_u, 0) ) *tgam_uu ;

  ss = sym_tmp / psi4  ;
        
  sym_tmp = contract(tgam_uu, 1, 
		     contract(tgam_uu, 1, dqq.derive_cov(ff), 0), 1) ;
                            
  sym_tmp.inc_dzpuis() ; // dzpuis : 3 --> 4

  tmp = qq.derive_con(ff).divergence(ff) ; 
  tmp.inc_dzpuis() ; // dzpuis : 3 --> 4  

  sym_tmp -= 0.3333333333333333 * tmp *tgam_uu ; 
                    
  ss -= sym_tmp / (psi4*ppsi*ppsi) ;  // Voir dans quel sens sont construits psi et psi4 (eviter les multiplications d'erreurs)

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += tgam_dd(k,l) * aa(i,k) * aa(j,l) ; 
	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
        
  tmp = psi4 * strain_tens.trace(ff) ; // S = S_i^i 

  ss += (2.*nn) * ( sym_tmp );
  Sym_tensor ss2 =2.*nn*( qpig*(psi4*strain_tens - 0.33333333333333 * tmp * tgam_uu));

  ss2.annule_domain(nz-1); //Dont care, i want only data on the horizon...

  ss += -ss2;
  
  // Source for h^{ij} 
  // -----------------
                

 Sym_tensor source_hh = 2.* nn * ss ;
  source_hh += 2. *  nn.derive_lie(beta) * aa  ;  //HERE

  // Term (d/dt - Lie_beta) (L beta)^{ij}--> sym_tmp        
  // ---------------------------------------

  sym_tmp = beta_point.ope_killing_conf(ff);
  sym_tmp.annule_domain(nz-1);
  sym_tmp = sym_tmp - l_beta.derive_lie(beta) ;
  
  sym_tmp.annule_domain(nz-1); 
  

  // Final source:
  // ------------
  source_hh += 0.6666666666666666* div_beta * l_beta - sym_tmp ; 

Scalar dNdt(hij.get_mp());

 dNdt = -source_hh(1,1)/2.*aa(1,1); // Take any component as long as the aa part does not vanish...

 dNdt.annule_domain(nz-1);

Scalar bound_N = old_nn + dNdt*delta_t; bound_N.std_spectral_base();
 bound_N.set_spectral_va().ylm();
 
 p_get_BC_lapse_4 = new Scalar(bound_N); 
    
}
return *p_get_BC_lapse_4 ;
}


}
