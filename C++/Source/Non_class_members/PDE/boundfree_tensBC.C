
// Header Lorene:
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "math.h"
#include "metric.h"
#include "param.h"
#include "param_elliptic.h"
#include "tensor.h"
#include "diff.h"
#include "proto.h"

namespace Lorene {
// Resolution of tensorial equation (N^2/Psi^4)Delta(hij) - LbLbhij = Sij, using degenerate elliptic solver.
// Here assumption is made that no boundary condition has to be enforced, mainly Beta^i*s_i = N/psi^2

Sym_tensor boundfree_tensBC ( Sym_tensor source, Vector Beta, Scalar Psi, Scalar Nn, Sym_tensor hij_guess, double precision, int loopmax) {
  cout << "================================================================" << endl;
    cout << "STARTING THE SUBITERATION FOR THE CONFORMAL METRIC" << endl;
    cout << "        iteration parameters are the following:        " << endl;
    cout << "        convergence precision required:" << precision << endl;
    cout << "        max number of global steps    :" << loopmax << endl;
    cout << "================================================================" << endl;

  // Construction of a multi-grid (Mg3d)
  // -----------------------------------
  const int nz = (*source.get_mp().get_mg()).get_nzone(); 	// Number of domains
  int nr = (*source.get_mp().get_mg()).get_nr(1); 	// Number of collocation points in r in each domain
  const Coord& rr = source.get_mp().r;
   Scalar rrr (source.get_mp()) ; 
  rrr = rr ; 
  rrr.std_spectral_base();  
 
  const Metric_flat& mets = (source.get_mp()).flat_met_spher() ;
 
  //// Initialisation of iteration variables.

  Sym_tensor hij  = hij_guess;    
  Sym_tensor hij_new = hij;

  Scalar n2sp4 = (Nn/(Psi*Psi))*(Nn/(Psi*Psi)) ; // Scale factor in front of Poisson Equation.
  n2sp4.std_spectral_base();

  // Resolution variables
  Scalar khi = hij(1,1);
   if (khi.get_etat() == ETATZERO) {
	khi.annule_hard() ;
	khi.set_dzpuis(4) ;
	khi.std_spectral_base() ;
      }
      khi.set_spectral_va().ylm();
  khi.mult_r();
  khi.mult_r_dzpuis(1);
  Scalar mmu = hij.mu();
   if (mmu.get_etat() == ETATZERO) {
	mmu.annule_hard() ;
	mmu.std_spectral_base() ;
      }
      mmu.set_spectral_va().ylm();
  Scalar etta = hij.eta();
   if (etta.get_etat() == ETATZERO) {
	etta.annule_hard() ;
	etta.std_spectral_base() ;
      }
      etta.set_spectral_va().ylm();
  Scalar Aa = hij.compute_A();
   if (Aa.get_etat() == ETATZERO) {
	Aa.annule_hard() ;
	Aa.std_spectral_base() ;
      }
      Aa.set_spectral_va().ylm();
  Scalar Bt = hij.compute_tilde_B();
   if (Bt.get_etat() == ETATZERO) {
	Bt.annule_hard() ;
	Bt.std_spectral_base() ;
      }
      Bt.set_spectral_va().ylm();
  Scalar hh = hij.trace(mets);
   if (hh.get_etat() == ETATZERO) {
	hh.annule_hard() ;
	hh.std_spectral_base() ;
      }

  //Fitting scalar
  Scalar fit1(source.get_mp()); fit1.set_etat_qcq(); 

  // For storing the result of inversion
  Scalar Aanew (source.get_mp()); Aanew.annule_hard(); Aanew.std_spectral_base();
  Scalar Btnew (source.get_mp()); Btnew.annule_hard(); Btnew.std_spectral_base();


  // Construction of sources for the next iteration

  Sym_tensor  LbLbhij = (hij.derive_lie(Beta)).derive_lie(Beta);
  hij.annule_domain(0);
  LbLbhij.annule_domain(0);
  LbLbhij.inc_dzpuis(1);

  Sym_tensor source_hij = source/n2sp4;
  Scalar sourcetrace = source_hij.trace(mets);
  Scalar Bttrace = source_hij.compute_tilde_B(); 
        
  Sym_tensor source_hij2 = LbLbhij/n2sp4;

  Scalar Btsource2 = source_hij2.compute_tilde_B();

  Scalar source_Bt2 = Bttrace + Btsource2;

  source_hij = source_hij + source_hij2;

  Scalar r2LbLbrr = LbLbhij(1,1);
  r2LbLbrr.mult_r();
  r2LbLbrr.mult_r();
  Scalar LbLbmu = LbLbhij.mu();

  Scalar source_khi2 = source_hij(1,1);
  source_khi2.mult_r();
  source_khi2.mult_r();
 
  Scalar source_mu2 = source_hij.mu();
  Scalar source_eta2 = source_hij.eta();

  Scalar source_A2 = source_hij.compute_A();
 

  source_khi2.annule_domain(0);
  source_mu2.annule_domain(0);
  source_A2.annule_domain(0);
  source_Bt2.annule_domain(0);
  source_eta2.annule_domain(0);
  
  //  source_A2.set_spectral_va().set_base( Aa.set_spectral_va().get_base());
  //  source_Bt2.set_spectral_va().set_base( Bt.set_spectral_va().get_base());

 

 //////////////// // Approximation of (Beta/(N^²/psi^4))^2, for input in degenerate operator parameters.
  Scalar Betacarre = (Beta(1)*Beta(1))/n2sp4 ;

  double  fitd1 = (Betacarre.val_grid_point(1,0,0,nr-1) - Betacarre.val_grid_point(1,0,0,0))/(rrr.val_grid_point(1,0,0,nr-1) - rrr.val_grid_point(1,0,0,0)) ;

//   double error = 0.; // Voluntary error on fit.  
//   fitd1 += error;
  

   int nrint = (nr-1)/2 ;
  double ampl_r = (rrr.val_grid_point(1,0,0, nr -1) - rrr.val_grid_point(1,0,0 ,0))/2.;
  
  Scalar approx(source.get_mp());
  approx.annule_hard();
  approx.std_spectral_base();
  


  fit1 = fitd1*(rrr-1.) +1.;
 // First order approximation
  approx.set_domain(1)= fit1.set_domain(1); // MAKE PARTICULAR FOR ETATUN; DECLARE FIT1?2?3 FIRST TO BE CLEAN.

 
 
  //Second order approximation
  Scalar firststep = Betacarre - approx;
  
  double ampli = firststep.val_grid_point(1,0,0,nrint);
  double  fit2d1 = - ampli/(ampl_r* ampl_r);
 

 
  approx.set_domain(1) += (fit2d1*(rrr - 1.)*(rrr - rrr.val_grid_point(1,0,0,nr-1))).set_domain(1);


  double fit0d2 = approx.val_grid_point(1,0,0,nr -1); 
  double fit1d2 = (Betacarre.val_grid_point(2,0,0,nr-1) - fit0d2)/(rrr.val_grid_point(2,0,0,nr-1)- rrr.val_grid_point(2,0,0,0));
  double fit0d3 = Betacarre.val_grid_point(3,0,0,0);
  double fit1d3 = ( - fit0d3)/(rrr.val_grid_point(3,0,0,nr-1)- rrr.val_grid_point(3,0,0,0));
  
  approx.set_domain(2) = (fit0d2 + fit1d2*(rrr - rrr.val_grid_point(2,0,0,0))).set_domain(2);
  approx.set_domain(3) = (fit0d3 + fit1d3*(rrr - rrr.val_grid_point(3,0,0,0))).set_domain(3);
  


  for(int ii=1; ii<=3; ii++){

   source_khi2.set_domain(ii) += (-approx*(khi.dsdr().dsdr())).set_domain(ii);

   source_mu2.set_domain(ii) += (-approx*(mmu.dsdr().dsdr())).set_domain(ii);

   source_eta2.set_domain(ii) += (-approx*(etta.dsdr().dsdr())).set_domain(ii);

   source_A2.set_domain(ii) += (-approx*(Aa.dsdr().dsdr())).set_domain(ii);

   source_Bt2.set_domain(ii) += (-approx*(Bt.dsdr().dsdr())).set_domain(ii);

  } 
  


  // Convergence markers
    Scalar Aa_old(source.get_mp());
    Scalar Bt_old(source.get_mp());
 
       
  // Parameters for the iteration
    cout <<"==================================================================================" << endl;
    cout << "amplitude for the tensor equation source (used as scaling for convergence marker)" << endl;
   cout <<"==================================================================================" << endl;
  double scale = max(maxabs(source));
  double diff_ent = 0.15 ; // Initialisation of the difference marker between two iterations on some value

/////////////////////////////////////////////////////////////////////////////////////
/////////////////   ITERATION ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
  
  for (int mer = 0; (diff_ent>precision) && (mer< loopmax) ; mer++) {
    double relax = 0.15; // Global relaxation parameter

// Résolution for variables A and tilde(B)
    tensorelliptic (source_A2, Aanew , fitd1, fit2d1, fit0d2, fit1d2, fit0d3, fit1d3);
    tensorellipticBt (source_Bt2, Btnew, fitd1, fit2d1, fit0d2, fit1d2, fit0d3, fit1d3);
 

    Aa = Aanew ;
    Bt = Btnew ; // No relaxation on these variables; it will be done globally on the tensor.
 
    /////////////////////////////////////////////////////////////////////////////////    
    //// Resolution tests

      Scalar essaiA = Aanew.laplacian() - approx*(Aanew.dsdr().dsdr());
      essaiA.annule_domain(0);
      Scalar bobofA = essaiA -source_A2 ;
      bobofA.set_spectral_va().ylm() ;
      //    maxabs (bobofA);
      //      bobofA.spectral_display();
      bobofA.set_spectral_va().ylm_i();
	 
      ////////////////////////////////////////////////////////////////////////////////////

      //-------------------------------------------
     // Retrieving the tensor hij by Dirac inversion
      //------------------------------------------

      /////////////
    // Magnetic part

     Scalar mmuA = mmu;
     Scalar mmuAsr = mmuA; mmuAsr.div_r();
     Scalar xxA(source.get_mp()); xxA.annule_hard(); xxA.std_spectral_base(); 
  
     const Scalar AA = Aa;
     Param* par_bc = 0x0;
     Sym_tensor_trans hijt(source.get_mp(), source.get_mp().get_bvect_spher(), mets);
     hijt = hij;
     hijt.sol_Dirac_A2 ( AA , mmuAsr, xxA, source_mu2, par_bc);


     // Monitoring  A and Hmu;
      Scalar musrsr = mmuAsr; musrsr.div_r_dzpuis(2);
      Aanew = xxA.dsdr() - musrsr;
      Scalar xxsr = xxA; xxsr.div_r_dzpuis(2);     
      //  Scalar Hmut = 3.*musrsr + mmuAsr.dsdr() + 2.*xxsr + xxsr.lapang(); 

      ////////////
    //Electric part

     Scalar hrrBC (source.get_mp()); hrrBC.annule_hard();
     Scalar wwBC(source.get_mp()); wwBC.annule_hard();
     Scalar tilde_etaBC(source.get_mp()); tilde_etaBC.annule_hard();

     Scalar source_khi3 = source_khi2;
     source_khi3.annule_domain(nz-1);
     source_khi3 += 2*hh;
     source_khi3.set_spectral_va().ylm();

      hijt.sol_Dirac_BC3( Bt, hh, hrrBC , tilde_etaBC , wwBC , source_khi3, source_eta2, par_bc, par_bc)  ;
  
      // Tensor reconstruction
      hij_new.set_auxiliary(hrrBC, tilde_etaBC, mmuAsr, wwBC, xxA, hh -hrrBC);
      
      hij= relax*hij_new + (1 - relax)*hij ; // Global relaxation (opposite to relaxation on resolution variables).
      
      // Calculation of updated trace.
    hh = (1 + hij(1,1))*( hij(2,3)*hij(2,3) - hij(2,2)*hij(3,3) )
        + hij(1,2)*hij(1,2)*(1 + hij(3,3)) 
        + hij(1,3)*hij(1,3)*(1 + hij(2,2)) 
        - hij(1,1)*(hij(2,2) + hij(3,3)) - 2*hij(1,2)*hij(1,3)*hij(2,3) ;

      khi = hij(1,1);
      khi.mult_r();
      khi.mult_r_dzpuis(0);    
      mmu = hij.mu();
      etta = hij.eta();
         
      Aa = hij.compute_A();  
      Bt = hij.compute_tilde_B();
 
      if (mer >=1){
	Aa.set_spectral_va().ylm();
	Bt.set_spectral_va().ylm();
	Aa_old.set_spectral_va().ylm();
	Bt_old.set_spectral_va().ylm();
	
	diff_ent = max(maxabs (Bt - Bt_old))/scale; // Convergence marker
      }
      Aa_old = Aa;
      Bt_old = Bt;
      

    // Update of sources for the next loop
     LbLbhij = (hij.derive_lie(Beta)).derive_lie(Beta);
     LbLbhij.inc_dzpuis(1);
     r2LbLbrr = LbLbhij(1,1);
     r2LbLbrr.mult_r();
     r2LbLbrr.mult_r();
     LbLbmu = LbLbhij.mu();
   
         
  source_hij = source/n2sp4;

  Bttrace = source_hij.compute_tilde_B();
  
  source_hij2 = LbLbhij/n2sp4;

  Btsource2 = source_hij2.compute_tilde_B();

  source_Bt2 = Bttrace + Btsource2;
  
  source_hij = source_hij + source_hij2;

  source_khi2 = source_hij(1,1);
  source_khi2.mult_r();
  source_khi2.mult_r();
  source_mu2 = source_hij.mu();
  source_eta2 = source_hij.eta();
  source_A2= source_hij.compute_A();
 
  
  source_A2.set_spectral_va().set_base( Aa.set_spectral_va().get_base());
  source_Bt2.set_spectral_va().set_base( Bt.set_spectral_va().get_base());
  
  for(int ii=1; ii<=3; ii++){
    source_khi2.set_domain(ii) += (-approx*(khi.dsdr().dsdr())).set_domain(ii);
    source_mu2.set_domain(ii) += (-approx*(mmu.dsdr().dsdr())).set_domain(ii);
    source_eta2.set_domain(ii) += (-approx*(etta.dsdr().dsdr())).set_domain(ii);
    source_A2.set_domain(ii) += (-approx*(Aa.dsdr().dsdr())).set_domain(ii);
    source_Bt2.set_domain(ii) += (-approx*(Bt.dsdr().dsdr())).set_domain(ii);
  } 
  source_khi2.annule_domain(0);  
  source_mu2.annule_domain(0);
  source_eta2.annule_domain(0);
  source_A2.annule_domain(0);
  source_Bt2.annule_domain(0);
  
  
  
//   cout << "real A resolution?" << endl;
//   Sym_tensor mucorrect  = LbLbhij + source;
//   mucorrect = mucorrect/n2sp4;
//   Scalar AAs = mucorrect.compute_A();
//   AAs.annule_domain(nz-1);
//   Aa.annule_domain(nz-1);
//   Scalar Aa2 = Aa.laplacian();
//   Scalar voirA = AAs - Aa2;
  
//   voirA.set_spectral_va().ylm_i();
//   //  des_meridian(voirA, 1.00001*Rhor, Rout,"diffA", 1); 
//   voirA.set_spectral_va().ylm();
  
  
//   maxabs (voirA);
//   voirA.spectral_display("voirA", 1.e-9);
  
  
  cout << "diff_ent" << diff_ent << endl;
  cout << mer << endl;
  
  
  }
  
  return hij;
 }





    
}
