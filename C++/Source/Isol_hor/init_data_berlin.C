/*
 *  Method of class Hor_isol to compute initial data following the Berlin
 *  prescription for the shift 
 *   
 *
 *    (see file isol_hor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
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

char init_data_berlin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/03/04 17:06:26  jl_jaramillo
 * Addition of the boost to the shift after solving the shift eqaution
 *
 * Revision 1.2  2005/03/03 10:04:28  f_limousin
 * The boundary conditions for the lapse, psi and shift are now
 * parameters (in file par_hor.d).
 *
 * Revision 1.1  2005/02/08 11:00:42  jl_jaramillo
 * Function to compute a single black hole with berlin boundary condition.
 *
 * Revision 1.4  2004/11/05 10:56:33  f_limousin
 * Delete arguments ener_dens, mom_dens and trace_stress in the function
 * init_data.
 *
 * Revision 1.1  2004/10/29 12:54:53  jl_jaramillo
 * First version
 *
 * Revision 1.4  2004/10/01 16:47:51  f_limousin
 * Case \alpha=0 included
 *
 * Revision 1.3  2004/09/28 16:10:05  f_limousin
 * Many improvements. Now the resolution for the shift is working !
 *
 * Revision 1.1  2004/09/09 16:41:50  f_limousin
 * First version
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
#include "cmp.h"
#include "tenseur.h"
#include "utilitaires.h"

void Isol_hor::init_data_berlin(int bound_psi, 
				int bound_beta, double precis, double relax, 
				int niter) {

    using namespace Unites ;
   
    // Initialisations
    // ---------------
    double ttime = the_time[jtime] ;    

    const Map& map = nn().get_mp() ; 
    const Base_vect& triad = *(nn().get_triad()) ;
    
    Scalar tmp(map) ;
    Scalar tmp_scal(map) ; 
    Vector tmp_vect(map, CON, triad) ;
    Sym_tensor tmp_sym(map, CON, triad) ;

    // Reset of quantities depending on K:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 

    set_aa( ( beta().ope_killing_conf(tgam()) +  gamt_point) / (2.* nn()) ) ; 

    
    ofstream conv("resconv.d") ; 
    conv << " # diff_nn   diff_psi   diff_beta " << endl ;
    
    // Iteration
    // ---------
    for (int mer=0; mer<niter; mer++) {
	
	//============================================
	//          Boundary conditions
	//============================================
      
	// Boundary for N	
	//---------------
	
	Valeur nn_bound (boundary_nn_Dir(0.2)) ;
	
	// Boundary for Psi (Neumann)
	//---------------------------
	
	Valeur psi_bound (boundary_psi_Neu_spat()) ;

      //=============================================
      // Resolution of elliptic equations
      //=============================================

	/*
	   
      // Resolution of the Poisson equation for the lapse
      // ------------------------------------------------
     
	Scalar nn_jp1 = source_nn().poisson_dirichlet(nn_bound, 0) + 1. ; 

	// Test:
	maxabs(nn_jp1.laplacian() - source_nn(),
	       "Absolute error in the resolution of the equation for N") ;  
	
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	//-----------
	if (mer>0)
	    nn_jp1 = relax * nn_jp1 + (1 - relax) * nn() ;
	double diff_nn = max( diffrel(nn(), nn_jp1) ) ;   
	n_evol.update(nn_jp1, jtime, ttime) ; 
	
	// Resolution of the Poisson equation for Psi
	// ------------------------------------------
	
	Scalar psi_jp1 = source_psi().poisson_neumann(psi_bound, 0) + 1. ;
	
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	//-----------
	
	// Test:
	maxabs(psi_jp1.laplacian() - source_psi(),
	       "Absolute error in the resolution of the equation for Psi") ;  
	

	psi_jp1 = relax * psi_jp1 + (1 - relax) * psi() ;
           
		
	*/


	// Resolution of the vector Poisson equation for the shift
	//---------------------------------------------------------



	Vector source_beta_old (source_beta() ) ;
	cout<<"source_beta(1) ="<<source_beta()(1).get_dzpuis()<<endl ;
	
	Vector beta_old(beta()) ; //in order to calculate the convergence with
	                          //the previous step


	// Resolution of the shift at THIS  time step
	//-------------------------------------------

	for (int int_b=0; int_b<niter; int_b++) {

	  Vector beta_int(beta()) ;

	  Scalar b_tilde_int =  b_tilde() ;

	  des_profile(b_tilde_int, 1.00001, 10, M_PI/2., 0., "b_tilde") ;
    

	  // Resolution for b_tilde
	  //-----------------------

	  //Corrected source of beta
	  //------------------------
	  double lambda = 0. ;
	  Vector source_reg = - (1./3. - lambda) * beta_int.divergence(ff).derive_con(ff) ;
	  source_reg.inc_dzpuis() ;
	  source_beta_old = source_beta_old + source_reg ;
	  



	  // Source for b_tilde
	  //-------------------
	  Scalar source_scal (map );


	  source_scal = contract( source_beta_old, 0,  met_gamt.radial_vect().down(0, met_gamt), 0) ;

	  tmp = - 2*contract (met_gamt.radial_vect().down(0, met_gamt), 0, 
			    contract( b_tilde_int.derive_cov(ff), 0, 
				      met_gamt.radial_vect().derive_con(ff), 1), 0) ;

	  source_scal += tmp ;

	  tmp = - b_tilde_int * contract(met_gamt.radial_vect().down(0, met_gamt), 0,
				     contract(met_gamt.radial_vect().derive_con(ff).
					      derive_cov(ff), 1, 2), 0) ;
	  tmp.inc_dzpuis() ;
	  source_scal += tmp ;
	  

	  // Boundary for b_tilde
	  //---------------------
	  Valeur b_tilde_bound =  boundary_b_tilde_Neu() ;

	  Scalar b_tilde_jp1 = source_scal.poisson_neumann(b_tilde_bound, 0)   ;

	  Scalar b_tilde_test = b_tilde_jp1 ;

	  // Test:
	  maxabs(b_tilde_test.laplacian() - source_scal,
		 "Absolute error in the resolution of the equation for b_tilde") ;  
	
	  
       

	   // Resolution for V^i
	  //-------------------
	  Vector vv_jp1 (beta()) ;

	  Vector vv_int = - b_tilde() * met_gamt.radial_vect() + beta() ;
	  

	  // Source for V^i
	  //---------------
	  
	  
	  //Effective source of V^i 
	  //-----------------------
	  Vector source_vv_eff = source_beta_old ;

	  
	  source_vv_eff -=  2*contract( b_tilde_jp1.derive_cov(ff), 0, 
				  met_gamt.radial_vect().derive_con(ff), 1) ;

	  tmp_vect = b_tilde_jp1 * contract(met_gamt.radial_vect().derive_con(ff).
					    derive_cov(ff), 1, 2) ;
	  tmp_vect.inc_dzpuis() ;

	  source_vv_eff -= tmp_vect ;


	  Scalar source_vv_eff_rad = contract(source_vv_eff, 0, 
					      met_gamt.radial_vect().down(0, met_gamt), 0) ;

	  source_vv_eff -= source_vv_eff_rad * met_gamt.radial_vect() ;

	  Vector source_vector (source_vv_eff) ;


	  Scalar sv_rad = contract( source_vector, 0,  met_gamt.radial_vect().down(0, met_gamt), 0) ;
	  sv_rad.dec_dzpuis(4) ;
	  cout<<"radial component of source of V^i"<<max(sv_rad)<<endl ;
	  //des_profile(sv_rad, 1.00001, 10, M_PI/2., 0., "sv_rad") ;
	  

	 
	  // Boundary values for V^i
	  //------------------------

	Valeur boundary_x (mp.get_mg()-> get_angu()) ;
	Valeur boundary_y (mp.get_mg()-> get_angu()) ;
	Valeur boundary_z (mp.get_mg()-> get_angu()) ;

	switch (bound_beta) {
	    
	    case 2 : {
		boundary_x = boundary_vv_x(omega) ;
		boundary_y = boundary_vv_y(omega) ;
		boundary_z = boundary_vv_z(omega) ;
		break ; 
	    }
	    default : {
		cout << "Unexpected type of boundary conditions for psi!" 
		     << endl 
		     << "  bound_psi = " << bound_psi << endl ; 
		abort() ;
		break ; 
	    }
		
	} // End of switch  

	if (boost_x != 0.) 
	    boundary_x -= beta_boost_x() ;
	if (boost_z != 0.) 
	    boundary_z -= beta_boost_z() ;
	

	  // Resolution
	  //-----------
	
	  double precision = 1e-8 ;
	  poisson_vect_boundary(lambda, source_vector, vv_jp1, boundary_x, 
			      boundary_y, boundary_z, 0, precision, 20) ;


	Vector boost_vect(mp, CON, mp.get_bvect_cart()) ;
	if (boost_x != 0.) {
	  boost_vect.set(1) = boost_x ;
	  boost_vect.set(2) = 0. ;
	  boost_vect.set(3) = 0. ;
	  boost_vect.std_spectral_base() ;
	  boost_vect.change_triad(mp.get_bvect_spher()) ;
	  vv_jp1 = vv_jp1 + boost_vect ;
	}
	  
	if (boost_z != 0.) {
	  boost_vect.set(1) = boost_z ;
	  boost_vect.set(2) = 0. ;
	  boost_vect.set(3) = 0. ;
	  boost_vect.std_spectral_base() ;
	  boost_vect.change_triad(mp.get_bvect_spher()) ;
	  vv_jp1 = vv_jp1 + boost_vect ;
	}


	  Vector vv_test (vv_jp1) ;
	
	  Scalar vv_rad = contract( vv_jp1, 0, met_gamt.radial_vect().down(0, met_gamt), 0) ;
	  cout<<"radial component of V^i"<<max(sv_rad)<<endl ;
	  des_profile(vv_rad, 1.00001, 10, M_PI/2., 0., "vv_rad") ;
	  


	  // Check of the resolution of V^i
	  // ------------------------------
	
	  source_vector.dec_dzpuis() ;
	  maxabs(vv_jp1.derive_con(ff).divergence(ff) 
		 + lambda * vv_jp1.divergence(ff)
		 .derive_con(ff) - source_vector,
		 "Absolute error in the resolution of the equation for V^i") ;  
	  cout << endl ;
	  
	 	 
	  // Relaxation (relax=1 -> new ; relax=0 -> old )  
	  //-----------
	  if (mer>0)
	    b_tilde_jp1 = relax * b_tilde_jp1 + (1 - relax) * b_tilde_int ;
	  
	  vv_jp1 = relax * vv_jp1 + (1 - relax) * vv_int ;
	  

	  // Intermediate update of beta 
	  //----------------------------

	  Vector beta_jp1 =  b_tilde_jp1 * met_gamt.radial_vect() + vv_jp1 ;
	  
	  Vector beta_test =  b_tilde_test * met_gamt.radial_vect() + vv_test ;

	 

	  //	  des_profile(beta_test(1), 1.00001, 10, M_PI/2., 0., "beta(1)") ;
	  //	  des_profile(beta_test(2), 1.00001, 10, M_PI/2., 0., "beta(2)") ;
	  //	  des_profile(beta_test(3), 1.00001, 10, M_PI/2., 0., "beta(3)") ;

	  // Check I of the resolution of beta
	  // -------------------------------
	
	  cout<<"------------------------------------------------"<<endl ;
	  source_beta_old.dec_dzpuis() ;
	  maxabs(contract(beta_test.derive_con(ff).derive_cov(ff), 1, 2)
		 //+ lambda * beta_test.divergence(ff).derive_con(ff) 
		 - source_beta_old,
		 "Absolute error in the resolution of the equation for beta") ;  
	  cout <<"------------------------------------------------"<< endl ;



	  // Check II of the resolution of beta
	  // ----------------------------------
	  Vector test_source ( source_vector) ;
	  
	  
	  test_source.inc_dzpuis() ;
	  test_source += 2*contract( b_tilde_jp1.derive_cov(ff), 0,  
				     met_gamt.radial_vect().derive_con(ff), 1) ;
	  test_source += source_scal * met_gamt.radial_vect() ;
	  tmp_vect = b_tilde_jp1 *
	    contract(met_gamt.radial_vect().derive_con(ff).derive_cov(ff), 1, 2) ;
	  tmp_vect.inc_dzpuis() ;

	  test_source +=tmp_vect ;
	  source_beta_old.inc_dzpuis() ;
	  cout<<"------------------------------------------------"<<endl ;
	  maxabs(test_source - source_beta_old ,
		 "Absolute error in the source") ;  
	  cout <<"------------------------------------------------"<< endl ;

	  arrete() ;

 





	  
	  //=============================================
	  //      Convergence control at THIS time step
	  //=============================================
	  cout<<""<<endl;
	  cout<<""<<endl;
	  double diff_b_tilde = max( diffrel(b_tilde_int, b_tilde_jp1) ) ;   
	  double diff_vv = max( maxabs(vv_jp1 - vv_int) ) ; 
	  double diff_beta = max( maxabs(beta_jp1 - beta_int) ) ; 
	  
	  cout << "step int= " << int_b << " :  diff_b_tilde = " << diff_b_tilde 
	       << "  diff_vv_int = " << diff_vv << endl  << "  diff_beta = " 
	       << diff_beta << endl  ;
	  cout << "----------------------------------------------" << endl ;
	  if ( (diff_b_tilde < precis) && (diff_vv < precis) )
	    break ; 
	  
	  //=============================================
	  //      Updates for next step at THIS time step
	  //=============================================
	
	  beta_evol.update(beta_jp1, jtime, ttime) ; 
	
	  // New value of A^{ij}:
	  
	  Sym_tensor aa_jp1 (map, CON, map.get_bvect_spher()) ;
	  int nnp = map.get_mg()->get_np(1) ;
	  int nnt = map.get_mg()->get_nt(1) ;
	  
	  int check ;
	  check = 0 ;
	  for (int k=0; k<nnp; k++)
	    for (int j=0; j<nnt; j++){
	      if (nn().val_grid_point(1, k, j , 0) < 1e-12){
		check = 1 ;
		break ;
	      }
	    }
	  
	  if (check == 0)
	    aa_jp1 = ( beta().ope_killing_conf(tgam()) + gamt_point ) 
	      / (2.* nn()) ;            
	  else {
	    Scalar nn_sxpun (division_xpun (Cmp(nn()), 0)) ;
	    
	    Sym_tensor aa_sxpun = beta().ope_killing_conf(tgam())+gamt_point ;
	    
	    for(int i=1; i<=3; i++)
	      for(int j=1; j<=i; j++){
		aa_sxpun.set(i,j).set_inner_boundary(1, 0) ;
		aa_sxpun.set(i,j) = Scalar (division_xpun 
					    (Cmp(aa_sxpun(i,j)), 0)) ;
	      }
	    aa_jp1 = aa_sxpun / (2*nn_sxpun) ;
	  }
	  
	  cout << "check = " << check << endl ;
	  
	  set_aa(aa_jp1) ; 
	  
	  




	 

	}
	/*
	// Source
	//-------
	
	Vector source_vector ( source_beta() ) ;
	double lambda = 0. ;
	Vector source_reg = - (1./3. - lambda) * beta().divergence(ff)
	    .derive_con(ff) ;
	source_reg.inc_dzpuis() ;
	source_vector = source_vector + source_reg ;
	
	// Boundary values
	//----------------
	
	Valeur boundary_x ( boundary_beta_x(ang_vel) ) ;
	Valeur boundary_y ( boundary_beta_y(ang_vel) ) ;
	Valeur boundary_z ( boundary_beta_z(ang_vel) ) ;
	
	// Resolution
	//-----------
	
	double precision = 1e-8 ;
	poisson_vect_boundary(lambda, source_vector, beta_jp1, boundary_x, 
			      boundary_y, boundary_z, 0, precision, 20) ;
	
	// Check of the resolution
	// ------------------------
	
	source_vector.dec_dzpuis() ;
	maxabs(beta_jp1.derive_con(ff).divergence(ff) 
	       + lambda * beta_jp1.divergence(ff)
	       .derive_con(ff) - source_vector,
	       "Absolute error in the resolution of the equation for beta") ;  
	cout << endl ;
	
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	//-----------
	
	beta_jp1 = relax * beta_jp1 + (1 - relax) * beta() ;
	
	//===========================================
	//      Convergence control
	//===========================================
	
//      double diff_nn = max( diffrel(nn(), nn_jp1) ) ;   
	double diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
	double diff_beta = max( maxabs(beta_jp1 - beta()) ) ; 
	
	cout << "step = " << mer << " :  diff_psi = " << diff_psi 
	     << "  diff_nn = " << diff_nn 
	     << "  diff_beta = " << diff_beta << endl ;
	cout << "----------------------------------------------" << endl ;
	if ( (diff_psi < precis) && (diff_nn < precis) && (diff_beta < precis) )
	    break ; 
	
	conv << mer << "  " << log(diff_nn) << " " << log(diff_psi) 
	     << " " << log(diff_beta) << endl ;
	
	//=============================================
	//      Updates for next step 
	//=============================================
	
	set_psi_del_q(psi_jp1) ; 
	n_evol.update(nn_jp1, jtime, ttime) ; 
	beta_evol.update(beta_jp1, jtime, ttime) ; 
	
	// New value of A^{ij}:
	
	Sym_tensor aa_jp1 (map, CON, map.get_bvect_spher()) ;
	int nnp = map.get_mg()->get_np(1) ;
	int nnt = map.get_mg()->get_nt(1) ;
	
	int check ;
	check = 0 ;
	for (int k=0; k<nnp; k++)
	    for (int j=0; j<nnt; j++){
		if (nn().val_grid_point(1, k, j , 0) < 1e-12){
		    check = 1 ;
		    break ;
		}
	    }
	
	if (check == 0)
	    aa_jp1 = ( beta().ope_killing_conf(tgam()) + gamt_point ) 
		/ (2.* nn()) ;            
	else {
	    Scalar nn_sxpun (division_xpun (Cmp(nn()), 0)) ;
	    
	    Sym_tensor aa_sxpun = beta().ope_killing_conf(tgam())+gamt_point ;
	    
	    for(int i=1; i<=3; i++)
		for(int j=1; j<=i; j++){
		    aa_sxpun.set(i,j).set_inner_boundary(1, 0) ;
		    aa_sxpun.set(i,j) = Scalar (division_xpun 
						(Cmp(aa_sxpun(i,j)), 0)) ;
		}
	    aa_jp1 = aa_sxpun / (2*nn_sxpun) ;
	}

	cout << "check = " << check << endl ;
	
	set_aa(aa_jp1) ; 
	
    }
    
	

    conv.close() ;

*/
    
    }

}



