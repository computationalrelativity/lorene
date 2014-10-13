
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
#include "isol_hole.h"

//Computes the rhs of hyperbolic equation for conformal metric assuming statioarity; WARNING; up to now, we are only able to handle void spacetimes.


namespace Lorene {
void Isol_hole::secmembre_kerr(Sym_tensor& source_hh){ 

  // Getting: hij; hatA; lapse; conf_fact; shift;
  // hij; aa; nn; ppsi; bb;


  using namespace Unites;

  const int nz = (*(hij.get_mp().get_mg())).get_nzone(); 
 

  const Vector& beta = shift;
  const Sym_tensor& hh  = hij;

 
  const Scalar& psi4 = conf_fact*conf_fact*conf_fact*conf_fact;
  Scalar ln_psi = log(conf_fact); ln_psi.std_spectral_base();
 
  const Scalar qq = lapse*conf_fact*conf_fact; 
   

  Sym_tensor aa = hatA/(psi4*sqrt(psi4));
  aa.std_spectral_base(); //(check...)
 
  
  const Metric_flat& ff = (hij.get_mp()).flat_met_spher() ;
   
  const Sym_tensor& tgam_uu = ff.con() + hh;
   
  const Metric tgam(tgam_uu);
    
  const Base_vect_spher& otriad = hij.get_mp().get_bvect_spher();

  // On met a zero les quantités supposees etre de "matiere"

  Sym_tensor strain_tens = hij; 
  for (int ii=1; ii<=3; ii++)
    for(int jj=1; jj<=3; jj++)
      { strain_tens.set(ii,jj).annule_hard();
      }
  strain_tens.std_spectral_base();
   
  //     On met a zero les quantités derivee temporelle

  Vector beta_point = shift;
 
  for (int ii=1; ii<=3; ii++)
   
    { beta_point.set(ii).annule_hard();
    }
  beta_point.annule_domain(nz-1) ; // Pour faire passer un assert, je ne comprends pas pourquoi...  

  
  
  beta_point.std_spectral_base();
  Scalar lapse_point(hij.get_mp());
  lapse_point.annule_hard();
  lapse_point.annule_domain(nz -1);
 
  lapse_point.std_spectral_base(); 
   
  Sym_tensor hh_point = hij; 
  for (int ii=1; ii<=3; ii++)
    for(int jj=1; jj<=3; jj++)
      { hh_point.set(ii,jj).annule_hard();
      }
  hh_point.annule_domain(nz-1); 
  hh_point.std_spectral_base();  


  // Note: Il sera probablement nécessaire de ne pas mettre a zero hh point;
   

  //Sym_tensor Rrij(map, CON, map.get_bvect_spher());


  // Rrij = 0.5*[ contract (( h_iju, 0,  h_iju.derive_cov(mets).derive_cov(mets), 3), 0,3) - contract((contract(h_iju.derive_cov(mets),1, h_iju.derive_cov(mets),2)), 1,3)] ;
				   
  //==================================
  // Source for hij
  //==================================
        
  const Sym_tensor& tgam_dd = tgam.cov() ;    // {\tilde \gamma}_{ij}
  //  const Sym_tensor& tgam_uu = tgam().con() ;    // {\tilde \gamma}^{ij}
  const Tensor_sym& dtgam = tgam_dd.derive_cov(ff) ;// D_k {\tilde \gamma}_{ij} // ff etant la métrique plate
  const Tensor_sym& dhh = hh.derive_cov(ff) ; // D_k h^{ij}
  const Vector& dln_psi = ln_psi.derive_cov(ff) ; // D_i ln(Psi)
  const Vector& tdln_psi_u = ln_psi.derive_con(tgam) ; // tD^i ln(Psi)
  const Vector& tdnn_u = lapse.derive_con(tgam) ;       // tD^i N
  const Vector& dqq = qq.derive_cov(ff) ;         // D_i Q
  const Scalar& div_beta = beta.divergence(ff) ;  // D_k beta^k

  Scalar tmp(hij.get_mp()) ;
  Sym_tensor sym_tmp(hij.get_mp(), CON, otriad) ; 

  // Quadratic part of the Ricci tensor of gam_tilde 
  // ------------------------------------------------
        
  Sym_tensor ricci_star(hij.get_mp(), CON, otriad) ; 
        
  ricci_star = contract(hh, 0, 1, dhh.derive_cov(ff), 2, 3) ; 

  ricci_star.inc_dzpuis() ;   // dzpuis : 3 --> 4

  // des_profile (ricci_star(1,1), 1, 8, 1,1, "riccistar");  // A enlever

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += dhh(i,k,l) * dhh(j,l,k) ; 
	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
  ricci_star -= sym_tmp ;
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  for (int m=1; m<=3; m++) {
	    for (int n=1; n<=3; n++) {
	      
	      tmp += 0.5 * tgam_uu(i,k)* tgam_uu(j,l) 
		* dhh(m,n,k) * dtgam(m,n,l)
		+ tgam_dd(n,l) * dhh(m,n,k) 
		* (tgam_uu(i,k) * dhh(j,l,m) + tgam_uu(j,k) *  dhh(i,l,m) )
		- tgam_dd(k,l) *tgam_uu(m,n) * dhh(i,k,m) * dhh(j,l,n) ;
	    }
	  } 
	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
  ricci_star += sym_tmp ;

  ricci_star = 0.5 * ricci_star ; 
        
  // Curvature scalar of conformal metric :
  // -------------------------------------
        
  Scalar tricci_scal = 
    0.25 * contract(tgam_uu, 0, 1,
		    contract(dhh, 0, 1, dtgam, 0, 1), 0, 1 ) 
    - 0.5  * contract(tgam_uu, 0, 1,
		      contract(dhh, 0, 1, dtgam, 0, 2), 0, 1 ) ;  
                                                       
  // Full quadratic part of source for h : S^{ij}
  // --------------------------------------------
        
  Sym_tensor ss(hij.get_mp(), CON, otriad) ;  // Source secondaire 
        
  sym_tmp = lapse * (ricci_star + 8.* tdln_psi_u * tdln_psi_u)
    + 4.*( tdln_psi_u * tdnn_u + tdnn_u * tdln_psi_u ) 
    - 0.3333333333333333 * 
    ( lapse * (tricci_scal  + 8.* contract(dln_psi, 0, tdln_psi_u, 0) )
      + 8.* contract(dln_psi, 0, tdnn_u, 0) ) *tgam_uu ;

  ss = sym_tmp / psi4  ;
        
  sym_tmp = contract(tgam_uu, 1, 
		     contract(tgam_uu, 1, dqq.derive_cov(ff), 0), 1) ;
                            
  sym_tmp.inc_dzpuis() ; // dzpuis : 3 --> 4
  //  des_profile (sym_tmp(1,1), 1, 8, 1,1, "sym_tmp");  // A enlever
  

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += ( hh(i,k)*dhh(l,j,k) + hh(k,j)*dhh(i,l,k)
		   - hh(k,l)*dhh(i,j,k) ) * dqq(l) ; 
	}
      }
      sym_tmp.set(i,j) += 0.5 * tmp ; 
    }
  }
        
  tmp = qq.derive_con(tgam).divergence(tgam) ; 
  tmp.inc_dzpuis() ; // dzpuis : 3 --> 4  // reverifier pourquoi
        
  //  des_profile (tmp, 1, 8, 1,1, "tmp");  // A enlever

  sym_tmp -= 0.3333333333333333 * tmp *tgam_uu ; 
                    
  ss -= sym_tmp / (psi4*conf_fact*conf_fact) ;  // Voir dans quel sens sont construits psi et psi4 (eviter les multiplications d'erreurs)

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
        
  tmp = psi4 * strain_tens.trace(tgam) ; // S = S_i^i  

  ss += (2.*lapse) * ( sym_tmp);// - qpig*( psi4* strain_tens 
  //    - 0.3333333333333333 * tmp * tgam_uu );
  Sym_tensor ss2 =2.*lapse*( qpig*(psi4*strain_tens - 0.33333333333333 * tmp * tgam_uu));
  ss2.inc_dzpuis(4); // A retirer!
 
  //  des_profile (ss2(1,1), 1, 8, 1,1, "ss2");  // A enlever

  // cout << zone << endl; 

  //  ss2.annule_domain(nz-1); 

  ss += -ss2; // ATTENTION!!!! A RETABLIR!!!!

  // maxabs(ss, "ss tot") ; 
  
  // Source for h^{ij} 
  // -----------------
                 
  Sym_tensor lbh = hh.derive_lie(beta) ; 

  source_hh =// (lapse*lapse/psi4 - 1.) 
    // * hh.derive_con(ff).divergence(ff) 
    + 2.* hh_point.derive_lie(beta); // - lbh.derive_lie(beta) ; // La double derivée de
  // Lie en Beta est retirée (prise en charge dans tensorelliptic.C)

  source_hh.inc_dzpuis() ; 
        
  //  des_profile (source_hh(1,1), 1, 8, 1,1, "sourcehh");  // A enlever

  source_hh += 2.* lapse * ss ;
              
  //## Provisory: waiting for the Lie derivative to allow
  //  derivation with respect to a vector with dzpuis != 0
  Vector vtmp = beta_point ; 
  // vtmp.dec_dzpuis(2) ;  // A remettre si jamais beta point est non nul;

  sym_tmp = hh.derive_lie(vtmp) ; 
  sym_tmp.inc_dzpuis(2) ;             

  //  des_profile (sym_tmp(1,1), 1, 8, 1,1, "sym_tmp");  // A enlever

  source_hh += sym_tmp 
    + 1.3333333333333333 * div_beta* (hh_point - lbh)
    + 2. * (lapse_point - lapse.derive_lie(beta)) * aa  ;
              

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0. ; 
      for (int k=1; k<=3; k++) {
	tmp += ( hh.derive_con(ff)(k,j,i) 
		 + hh.derive_con(ff)(i,k,j) 
		 - hh.derive_con(ff)(i,j,k) ) * dqq(k) ;
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
            
  source_hh -= lapse / (psi4*conf_fact*conf_fact) * sym_tmp ; 
         
  tmp =  beta_point.divergence(ff) - div_beta.derive_lie(beta) ; 
  tmp.inc_dzpuis() ; 

  //  des_profile (tmp, 1, 8, 1,1, "tmp");  // A enlever

  source_hh += 0.6666666666666666* 
    ( tmp - 0.6666666666666666* div_beta * div_beta ) * hh ; 
               
        
  // Term (d/dt - Lie_beta) (L beta)^{ij}--> sym_tmp        
  // ---------------------------------------
  Sym_tensor l_beta = beta.ope_killing_conf(ff) ; // Attention aux headers a inclure

  sym_tmp = beta_point.ope_killing_conf(ff) - l_beta.derive_lie(beta) ;
  
  sym_tmp.inc_dzpuis() ; 
  
  // Final source:
  // ------------
  source_hh += 0.6666666666666666* div_beta * l_beta - sym_tmp ; 

  // Invert it (because the right source is for a Laplace operator, not a Dalembertain operator as it is initially:

  source_hh = -source_hh;
 
   
  // Annulation de la source
//   for (int ii=1; ii<=3; ii++)
//     for (int jj=1; jj<=3; jj++){
//       source_hh.set(ii,jj).annule_hard();
//     }
//   source_hh.std_spectral_base();
 
  return; 
 	                       
}	       



}
