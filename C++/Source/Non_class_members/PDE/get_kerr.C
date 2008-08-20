
// Header Lorene:
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "math.h"
#include "metric.h"
#include "param.h"
#include "param_elliptic.h"
#include "vector.h"
#include "scalar.h"
#include "diff.h"
#include "proto.h"

void get_kerr(double rot_par, Sym_tensor hij, Scalar nn, Scalar ppsi, Vector bb, double bound_n, bool non_conf_flat, double precis) {
  // Implementation of FCF system with the solver sol_elliptic_boundary, and axisymmetric isolated horizon boundary conditions. Data are meant to be stationary, in maximal slicing gauge.
  // This version uses variables Psi, N*Psi, and Beta, in non conformally flat case. Actualizations of variables are done during the iteration. Here, extrinsic curvature is scaled in psi10.
  // Here is also implemented dealing of trace of hij, regarding det=1 constraint.
  // Construction of a multi-grid (Mg3d)
  // -----------------------------------
  const Map_af* map = dynamic_cast<const Map_af*>(&hij.get_mp()) ;

  const Mg3d* mgrid = (*map).get_mg();
	
  // Construct angular grid for h(theta,phi) 
  const Mg3d& g_angu = *(*mgrid).get_angu_1dom() ;
  

  const int nz = (*hij.get_mp().get_mg()).get_nzone(); 	// Number of domains
  int nt = (*hij.get_mp().get_mg()).get_nt(1); 	// Number of collocation points in theta in each domain
  const int np = 1; // Axisymmetry, then only one point in phi.
  const Coord& rr = hij.get_mp().r;
   Scalar rrr (hij.get_mp()) ; 
  rrr = rr ; 
  rrr.std_spectral_base();  
  assert((rrr.val_grid_point(1,0,0,0) - 1.) <= 1.e-9); // For now the code handles only horizons at r=1
  double r_limits2[] = {rrr.val_grid_point(1,0,0,0), rrr.val_grid_point(2,0,0,0)} ; 
  const Map_af map_2(g_angu, r_limits2);
 
  const Metric_flat& mets = (hij.get_mp()).flat_met_spher() ;
 

  Scalar logn (*map) ; logn = log(nn) ; 
  logn.std_spectral_base();


  Scalar logpsi(*map) ; logpsi = log(ppsi) ;
  logpsi.std_spectral_base();
  Scalar psi4 (*map) ; psi4 = ppsi*ppsi*ppsi*ppsi ;
  Scalar npsi (*map) ;  npsi =ppsi*nn ;
  

  Scalar ppsi_new(*map) ; ppsi_new.annule_hard(); ppsi_new.std_spectral_base(); 
  Scalar npsi_new(*map); npsi_new.annule_hard(); npsi_new.std_spectral_base();

  Vector bb_new (*map, CON, (*map).get_bvect_spher()); 
  for(int i=1; i<=3; i++){
    bb_new.set(i)=0;
  }
  bb_new.std_spectral_base();
 
  // Non conformally flat variables

 
  Sym_tensor gamtuu = mets.con() + hij; 
  Metric gamt(gamtuu);
  Metric gam(gamt.cov()*psi4) ;
  Sym_tensor gamma = gam.cov();
  

  // Extrinsic curvature variables

  Sym_tensor aa(*map, CON, (*map).get_bvect_spher());
  for (int iii= 1; iii<=3; iii++){ 
    for(int j=1; j<=3; j++){
      aa.set(iii,j)= 0;
    }
  }
  aa.std_spectral_base(); 
  Scalar aa_quad_scal(*map) ; aa_quad_scal = 0. ;

  Sym_tensor aa_hat(*map, CON, (*map).get_bvect_spher());
  for (int iii= 1; iii<=3; iii++){ 
    for(int j=1; j<=3; j++){
      aa_hat.set(iii,j)= 0;
    }
  }
 
  Sym_tensor kuu = aa/psi4 ;
  Sym_tensor kuu2 = aa_hat/(psi4*psi4*sqrt(psi4));
  Sym_tensor kdd = contract (gamma, 0, contract(gamma, 1, kuu, 0),1);
 

  // (2,1)-rank delta tensor: difference between ricci rotation coefficients. 

   Tensor delta = -0.5*contract( hij, 1, gamt.cov().derive_cov(mets), 2);
   Scalar tmp(*map);

   for (int i=1; i<=3; i++) {
     for (int j=1; j<=3; j++) {
       for (int k=1; k<=3; k++) {
	 tmp = 0.;
	 tmp =  -0.5 *(gamt.cov().derive_con(mets))(i,j,k);
	 for (int l=1; l<=3; l++) {
	   tmp += -0.5*( gamt.cov()(i,l)*(hij.derive_cov(mets))(k,l,j) + gamt.cov()(l,j)*(hij.derive_cov(mets))(k,l,i));
	 }
       	 delta.set(k,i,j) += tmp ; 
       }
     }
   }
   

   // Conformal Rstar scalar(eq 61, Bonazzola et al. 2003) 
   Scalar Rstar = 
    0.25 * contract(gamt.con(), 0, 1,
		    contract(hij.derive_cov(mets), 0, 1, gamt.cov().derive_cov(mets), 0, 1), 0, 1 ) 
    - 0.5  * contract(gamt.con(), 0, 1,
		      contract(hij.derive_cov(mets), 0, 1, gamt.cov().derive_cov(mets), 0, 2), 0, 1 ) ; 
  
   Scalar norm(*map);
   Scalar norm3(*map);

   
  // Parameters for the iteration 

  int mer_max = 5000 ; // 5000 iterations maximum
  double diff_ent = 1 ; // initialization of difference marker between two iterations; 

  int util = 0; // Tool used to stop tensorial iteration at any wished step "util"
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////// ITERATION   ////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  for(int mer=0 ;(diff_ent > precis) && (mer<mer_max) ; mer++) {    

   //Global relaxation coefficient

    double relax = 0.1;

    // Scalar variables linked to the norm of normal vector to horizon.
    norm = sqrt(1. + hij(1,1)); norm.std_spectral_base();
    norm3 = sqrt(1. + hij(3,3)); norm3.std_spectral_base();
    


     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



      ///////////////////////////
      // Solving for (Psi-1) //
      ///////////////////////////

  
      // Setting of the boundary

      double   diric= 0.; 
      double   neum = 1.; 

      Vector ssalt = rrr.derive_cov(gam);
      Vector ssaltcon = ssalt.up_down(gam);
      Scalar ssnormalt = sqrt(contract (ssalt,0, ssaltcon, 0));
      ssnormalt.std_spectral_base();
      
      ssalt.annule_domain(nz-1);
      ssalt.annule_domain(0);
      ssaltcon.annule_domain(nz-1);
      ssaltcon.annule_domain(0);
      
      ssalt = ssalt/ssnormalt;
      ssaltcon = ssaltcon/ssnormalt;
      Vector ssconalt = ssaltcon*ppsi*ppsi;    // \tilde{s} in the notations of Gourgoulhon and Jaramillo, 2006
      ssconalt.std_spectral_base();
      ssconalt.annule_domain(nz-1);
      Scalar bound3bis =   -((1./ppsi)*contract((contract(kdd,1,ssconalt,0)),0, ssconalt,0));
  
      bound3bis.annule_domain(nz-1);
      bound3bis += -ppsi*ssconalt.divergence(gamt);
      bound3bis.annule_domain(nz-1);     
      bound3bis = 0.25*bound3bis;
      bound3bis += -contract(ppsi.derive_cov(gamt),0,ssconalt,0) + ppsi.dsdr();
      bound3bis.annule_domain(nz-1);
      bound3bis.std_spectral_base();
      bound3bis.set_spectral_va().ylm();    

      Mtbl_cf *boundd3bis = bound3bis.set_spectral_va().c_cf;
            

     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



      // Computing the source
   
      Scalar source_ppsi(*map) ; source_ppsi=3. ; // Pour le fun... 
      source_ppsi.std_spectral_base();                            

      Scalar d2logpsi = contract(ppsi.derive_cov(mets).derive_cov(mets), 0, 1, hij, 0,1);  
      d2logpsi.inc_dzpuis(1);
 
      source_ppsi = -(0.125* aa_quad_scal )/(psi4*ppsi*ppsi*ppsi) +  ppsi* 0.125* Rstar - d2logpsi; 
      
      source_ppsi.std_spectral_base(); 
      if (source_ppsi.get_etat() == ETATZERO) {
	source_ppsi.annule_hard() ;
	source_ppsi.set_dzpuis(4) ;
	source_ppsi.std_spectral_base() ;
      }
      source_ppsi.set_spectral_va().ylm();
      
      
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // System inversion   

      Param_elliptic source11(source_ppsi);
      ppsi_new = source_ppsi.sol_elliptic_boundary(source11, *boundd3bis, diric , neum) + 1 ; // Resolution has been done for quantity Q-1, because our solver gives a vanishing solution at infinity! 
      

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
    // tests for resolution

    Scalar baba2 = (ppsi_new-1).laplacian();
//     cout << "psi+1-résolution" << endl;
//     maxabs (baba2 - source_ppsi);
   
    Scalar psinewbis = ppsi_new -1. ; psinewbis.annule_domain(nz -1);
    psinewbis.std_spectral_base();
    psinewbis = psinewbis.dsdr();
    Scalar psinewfin2 (map_2) ;
    psinewfin2.allocate_all(); 
    psinewfin2.set_etat_qcq();  
    psinewfin2.std_spectral_base();
    
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
	
	psinewfin2.set_grid_point(0, k, j, 0) = psinewbis.val_grid_point(1, k,j,0) - bound3bis.val_grid_point(1, k, j, 0);
	
      }
    
    //  maxabs (psinewfin2);
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Actualization during the loop

    ppsi = ppsi_new* (1-relax) + ppsi* relax ;
    psi4 = ppsi*ppsi*ppsi*ppsi;
    logpsi = log(ppsi) ; 
    logpsi.std_spectral_base();	   




     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    //////////////////////////    
    // Solving for (N*Psi -1)/ 
    //////////////////////////


    // Setting of the boundary  
    Scalar bound(*map); 
    bound = (bound_n)*ppsi -1;
    bound.annule_domain(nz -1);
    bound.std_spectral_base();
    bound.set_spectral_va().ylm();
    Mtbl_cf *boundd = bound.get_spectral_va().c_cf;

     diric =1; 
     neum = 0 ; 


    
     /////////////////////////////////////////////////////////////////////////////////////////////////:

    // Computing the source ...      
         Scalar d2lognpsi = contract(npsi.derive_cov(mets).derive_cov(mets), 0, 1, hij, 0,1);
       d2lognpsi.inc_dzpuis(1); //  dzpuis correction.
  
       Scalar source_npsi = npsi*(aa_quad_scal*(7./8.)/(psi4*psi4) + Rstar/8.) - d2lognpsi; // ?
    source_npsi.std_spectral_base();
    if (source_npsi.get_etat() == ETATZERO) {
      source_npsi.annule_hard() ;
      source_npsi.set_dzpuis(4) ;
      source_npsi.std_spectral_base() ;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////::

    // Inversion of the operator
    Param_elliptic source1 (source_npsi); 
    npsi_new = source_npsi.sol_elliptic_boundary(source1, *boundd, diric, neum) ;

    npsi_new = npsi_new +1;

  
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:: 

    // Resolution tests in npsi
    Scalar baba = npsi_new.laplacian();
//       cout << "resolution_npsi" << endl;
//      maxabs (baba - source_npsi);


    //  cout << "bound_npsi" << endl;
        Scalar npsibound2 (map_2) ;
    npsibound2.allocate_all(); 
    npsibound2.set_etat_qcq();  
    npsibound2.std_spectral_base();		 
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
		     
	npsibound2.set_grid_point(0, k, j, 0) = npsi_new.val_grid_point(1, k,j,0) - bound.val_grid_point(1, k, j, 0) -1.;
		     
      }

    //   maxabs (npsibound2);
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:: 

    // Actualisation during the loop

      
      npsi = npsi_new*(1-relax) + npsi* relax; 
       nn = npsi/ppsi; 
       logn = log(nn);
       logn.std_spectral_base(); 

 

      
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:: 
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:: 
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:: 

      ///////////////////////
      //Resolution in Beta //
     ////////////////////////


      // Setting of the boundary  
           
      bound = (bound_n)/(ppsi*ppsi) ;
      bound.annule_domain(nz -1);
 

      Scalar Omega(*map); Omega.annule_hard(); Omega = rot_par; // Rotation parameter for spacetime
      Omega.std_spectral_base() ; Omega.mult_rsint();
      Omega.annule_domain(nz -1);

      Vector limit = bb_new;
      Vector ephi(*map, CON, (*map).get_bvect_spher());
      ephi.set(1).annule_hard();
      ephi.set(2).annule_hard();
      ephi.set(3) = 1;
      ephi.std_spectral_base();
      ephi.annule_domain(nz -1);

      limit = bound*ssconalt + Omega*ephi;
      limit.std_spectral_base(); // Boundary is fixed by value of 3 components of a vector (rather than value of potentials)   

      Scalar Vrb = limit(1);  Vrb.set_spectral_va().ylm();
      Scalar mmub = limit.mu(); mmub.set_spectral_va().ylm(); 
      Scalar etab = limit.eta(); etab.set_spectral_va().ylm();



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////::  
 
     // Computing the source
       Vector deltaA =  - 2*nn*contract(delta, 1,2, aa, 0,1);
       Vector hijddb =  - contract (bb.derive_cov(mets).derive_cov(mets), 1,2, hij, 0,1) ;
       Vector hijddivb =  - 0.3333333333333* contract (bb.divergence(mets).derive_cov(mets),0, hij,1);
       hijddb.inc_dzpuis(); // dzpuis fixing patch... 
       hijddivb.inc_dzpuis(); 
       
       Vector sourcevect2(*map,CON, (*map).get_bvect_spher()); 
      sourcevect2= (2.* contract(aa, 1, nn.derive_cov(mets),0)-12*nn*contract(aa, 1, logpsi.derive_cov(mets), 0))  + deltaA + hijddb + hijddivb ; 
       
      sourcevect2.set(1).set_dzpuis(4);
      sourcevect2.set(2).set_dzpuis(4);
      sourcevect2.set(3).set_dzpuis(4);
      sourcevect2.std_spectral_base(); 
      if(sourcevect2.eta().get_etat() == ETATZERO)
	{ sourcevect2.set(2).annule_hard();}
 
      
      double lam = (1./3.);    
    
   

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // System inversion
      sourcevect2.poisson_boundary2(lam, bb_new, Vrb, etab, mmub, 1., 0., 1. ,0. ,1. ,0.) ;   



  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////::       
    
    // resolution tests
    Vector source2 = contract(bb_new.derive_con(mets).derive_cov(mets), 1,2) + lam* contract(bb_new.derive_cov(mets), 0,1).derive_con(mets);
    source2.inc_dzpuis(1);
    //  maxabs (source2 - sourcevect2);   

    Scalar mufin = bb_new.mu();
    mufin.set_spectral_va().coef();
		 
    Scalar mufin2 (map_2) ;
    mufin2.allocate_all(); 
    mufin2.set_etat_qcq();  
    mufin2.std_spectral_base();
		 
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
		     
	mufin2.set_grid_point(0, k, j, 0) = mufin.val_grid_point(1, k,j,0) - mmub.val_grid_point(1, k, j, 0);
		     
      }

    //  maxabs (mufin2);


    Scalar brfin = bb_new(1);
    brfin.set_spectral_va().coef();
		 
    Scalar brfin2 (map_2) ;
    brfin2.allocate_all(); 
    brfin2.set_etat_qcq();  
    brfin2.std_spectral_base();
		 
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++) {
		     
	brfin2.set_grid_point(0, k, j, 0) = brfin.val_grid_point(1, k,j,0) - Vrb.val_grid_point(1, k, j, 0);
		     
      }
    //  maxabs (brfin2);



  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////::       

    // Actualisation during the loop
    for (int ii=1; ii <=3; ii++){
	 bb.set(ii) = bb_new(ii)*(1-relax) + bb(ii)* relax;
    }
       

    diff_ent = max(maxabs(npsi_new - npsi )); // Convergence parameter (discutable relevance...) 



  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////::       
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////::       
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////::       


  
    ////////////////////////////////////
    // Tensor hij resolution        ///
    ///////////////////////////////////

    if (non_conf_flat){
   
    if (diff_ent <=5.e-7) { // No resolution until we are close to the result.

      util = util+1; // Loop marker for NCF equation.
     

  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////// ITERATION   ////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

	 // Local convergence can be asked for hij equation. Allows to satisfy integrability conditions.
  int mer_max2 = 1 ; // maximum iterations allowed
  double precis2 =  1.e5*precis ; // Here we ask for a local convergence; this can be improved later.
  double diff_ent2 = 1 ; // Local convergence marker
  double relax2 = 1.; // Local relaxation parameter. If not 1, the determinant condition won't be satisfied on a particular iteration.




  Sym_tensor sourcehij = hij; // Random initialization...

  for(int mer2=0 ;(diff_ent2 > precis2) && (mer2<mer_max2) ; mer2++) {    

   /////////////////////////////////////////////////////////////////////////////////////////////////

    // Calculation of the source

    // ATTENTION!!!!!
    // Sign in front of the source comes from the fact that it is calculated for a poisson-like operator, and not dalembertian.  
    // The double Lie derivative term is taken care of in the subroutine. 


 
        sourcehij =  -secmembre_kerr (hij, aa, nn, ppsi, bb);   
   
   /////////////////////////////////////////////////////////////////////////////////////////////////
	//System inversion (note that no boundary condition is imposed)


     Sym_tensor hij_new = hij;
 
     hij_new = boundfree_tensBC (sourcehij, bb , ppsi, nn, hij);
     diff_ent2 = max(maxabs(hij - hij_new));
   
      hij = relax2*hij_new + (1 - relax2)*hij;


      cout << "mer2, diffent2" << endl;
 
      cout << mer2 << endl;
      cout << diff_ent2 << endl;


  }

   /////////////////////////////////////////////////////////////////////////////////////////////////
   // Resolution tests


       Sym_tensor gammatilde = mets.con() + hij;   
       Metric gammatilde2(gammatilde); Scalar detgam = gammatilde2.determinant();
 //       cout << "determinant of result" << endl;
//       maxabs (detgam-1.);
   
//      cout << "comment l'equation en hij est elle vérifiée?" << endl;
        
      Sym_tensor test =contract (hij.derive_cov(mets).derive_con(mets), 2,3);
      test.annule(nz-1, nz-1);
      test = test - hij.derive_lie(bb).derive_lie(bb)/ ((nn/(ppsi*ppsi))*(nn/(ppsi*ppsi)));        
      test.annule(nz-1, nz-1);
      Sym_tensor youps = test - sourcehij/((nn/(ppsi*ppsi))*(nn/(ppsi*ppsi)));
//       maxabs (youps); 
//       maxabs((youps).trace(mets));
//       cout << " AAABBB" << endl;
//       maxabs((youps).compute_A());
//       maxabs((youps).compute_tilde_B());

	
    }
    }
 



   /////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////  

    

     ////////////////////////////////////////////////////
     // Global variable update after an entire loop ////
     ///////////////////////////////////////////////////

   
       gamtuu = mets.con() + hij;
       gamt = gamtuu; // Métrique.
       gam = gamt.cov()*psi4;
       gamma = gam.cov();
       

       for (int i=1; i<=3; i++) {
	 for (int j=1; j<=3; j++) {
	   
	   tmp = 0;       
	   tmp = ((bb.derive_con(mets))(i,j) + (bb.derive_con(mets))(j,i) - (2./3.)*bb.divergence(mets) * mets.con()(i,j))*(1./(2.*nn));	   
	   aa.set(i,j) = tmp ; 
	   
	 }
       }

       aa = aa -  (hij.derive_lie(bb) + (2./3.)*bb.divergence(mets)*hij)*(1./(2.*nn)); //Non conformally flat correction; we suppose here dhij/dt = 0.
  
       aa_hat = aa*psi4*sqrt(psi4); // Rescaling of traceless exrinsic curvature.
       aa_hat.std_spectral_base();

       Sym_tensor aaud = aa.up_down(gamt);
       Sym_tensor aaud_hat = aa_hat.up_down(gamt);
       aa_quad_scal =  contract(contract (aa_hat, 0, aaud_hat, 0), 0,1);

       delta = -0.5*contract( hij, 1, gamt.cov().derive_cov(mets), 2);
 

   for (int i=1; i<=3; i++) {
     for (int j=1; j<=3; j++) {
       for (int k=1; k<=3; k++) {
	 tmp = 0.;
	 tmp =  -0.5 *(gamt.cov().derive_con(mets))(i,j,k);
	 for (int l=1; l<=3; l++) {
	   tmp += -0.5*( gamt.cov()(i,l)*(hij.derive_cov(mets))(k,l,j) + gamt.cov()(l,j)*(hij.derive_cov(mets))(k,l,i));
	 }
       	 delta.set(k,i,j) += tmp ; 
       }
     }
   } 
   
  Rstar = 
    0.25 * contract(gamt.con(), 0, 1,
		    contract(gamt.con().derive_cov(mets), 0, 1, gamt.cov().derive_cov(mets), 0, 1), 0, 1 ) 
    - 0.5  * contract(gamt.con(), 0, 1,
		      contract(gamt.con().derive_cov(mets), 0, 1, gamt.cov().derive_cov(mets), 0, 2), 0, 1 ) ;  

    kuu = aa/(psi4); 
    kdd =  contract (gamma, 0, contract(gamma, 1, kuu, 0),1);  
 


  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Convergence markers
    cout << "diffent" << endl; 
    cout<< diff_ent << endl;   
    cout <<"mer" << mer << endl;   
   



   /////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////  
     
  //------------------------------------------------
  //     Check of Einstein equations (3+1 form)
  //------------------------------------------------


  // Lapse 
  //------
  Scalar lapse(*map) ;
  lapse = nn ;
  lapse.std_spectral_base();

  // 3-metric
  //---------
  //   const Metric gam(mets.cov()*psi4) ;
   
  Sym_tensor gamt2(*map, COV, (*map).get_bvect_spher()); 
  for (int i=1; i<=3; i++)
    for (int j=1; j<=3; j++) 
      { gamt2.set(i,j)=gam.cov()(i,j); 
      }
  //Shift
  //-----
  Vector beta = bb ;
  
  // Extrinsic curvature
  //--------------------

  Scalar TrK3(*map);
  Sym_tensor k_uu = aa/(psi4) ;
  k_uu.dec_dzpuis(k_uu(1,1).get_dzpuis()); 
  Sym_tensor k_dd = k_uu.up_down(gam); // Another way of computing the same thing, just to be sure... 

  TrK3 = k_uu.trace(gam);
  //  TrK3.spectral_display("TraceKvraie", 1.e-10);
 
  // Hamiltonian constraint
  //-----------------------
  Scalar ham_constr = gam.ricci_scal() ;
  ham_constr.dec_dzpuis(3) ;
  ham_constr +=  TrK3*TrK3 - contract(k_uu, 0, 1, k_dd, 0, 1) ;
  // maxabs(ham_constr, "Hamiltonian constraint: ") ;

  ham_constr.set_spectral_va().ylm();
  //  ham_constr.spectral_display("ham_constr", 1.e-9);
  

 
  // Momentum constraint
  //-------------------
  Vector mom_constr = k_uu.divergence(gam)  - TrK3.derive_con(gam) ;
  mom_constr.dec_dzpuis(2) ;
//   maxabs(mom_constr, "Momentum constraint: ") ;
//   mom_constr(1).spectral_display("mom1", 1.e-9) ;
//  mom_constr(2).spectral_display("mom2", 1.e-9) ;
//  mom_constr(3).spectral_display("mom3", 1.e-9) ;
  
 
  // Evolution equations
  //--------------------
  Sym_tensor evol_eq = lapse*gam.ricci() 
    - lapse.derive_cov(gam).derive_cov(gam);
  evol_eq.dec_dzpuis() ;
  evol_eq += k_dd.derive_lie(beta) ;
  evol_eq.dec_dzpuis(2) ;
  evol_eq += lapse*(TrK3*k_dd - 2*contract(k_dd, 1, k_dd.up(0, gam), 0) ) ;
//   maxabs(evol_eq, "Evolution equations: ") ;
  
//   evol_eq.trace(gam).spectral_display("evoltrace", 1.e-10);
//   maxabs (evol_eq.trace(gam));
 
 
  //  evol_eq.spectral_display("evol", 1.e-10);
  


  

  }

 
  return; 
}
 

 

  
