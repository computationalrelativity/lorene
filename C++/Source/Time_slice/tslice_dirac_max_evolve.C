/*
 *  Method of class Tslice_dirac_max for time evolution 
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon & Jerome Novak
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

char tslice_dirac_max_evolve_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.19  2012/02/06 12:59:07  j_novak
 * Correction of some errors.
 *
 * Revision 1.18  2011/07/22 13:21:02  j_novak
 * Corrected an error on BC treatment.
 *
 * Revision 1.17  2010/10/20 07:58:10  j_novak
 * Better implementation of the explicit time-integration. Not fully-tested yet.
 *
 * Revision 1.16  2008/12/04 18:22:49  j_novak
 * Enhancement of the dzpuis treatment + various bug fixes.
 *
 * Revision 1.15  2008/12/02 15:02:22  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.13  2004/06/14 20:51:37  e_gourgoulhon
 * Method solve_hij has now argument method_poisson.
 * Its value is set **provisory** to 1 (instead of method_poisson_vect !).
 *
 * Revision 1.12  2004/05/31 09:09:59  e_gourgoulhon
 * Added monitoring of khi and mu.
 * Added writing of whole configuration in file (via Time_slice::save).
 *
 * Revision 1.11  2004/05/24 20:58:05  e_gourgoulhon
 * Added graphical output of khi, mu and trh.
 *
 * Revision 1.10  2004/05/20 20:32:01  e_gourgoulhon
 * Added arguments check_mod and save_mod.
 * Argument graph_device passed to des_evol.
 *
 * Revision 1.9  2004/05/17 19:55:10  e_gourgoulhon
 * Added arguments method_poisson_vect, nopause and graph_device
 *
 * Revision 1.8  2004/05/13 21:35:30  e_gourgoulhon
 * Added monitoring of various quantities (as Evolution_full<Tbl>).
 * Added function monitor_scalar.
 *
 * Revision 1.7  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.6  2004/05/11 20:15:10  e_gourgoulhon
 * Added Evolution_full's for ADM mass and checks of the constraint,
 * as well as the corresponding plots and write to files.
 *
 * Revision 1.5  2004/05/10 09:19:27  e_gourgoulhon
 * Added a call to del_deriv() after set_khi_mu.
 *
 * Revision 1.4  2004/05/09 20:59:06  e_gourgoulhon
 * Change of the time scheme: first solve d'Alembert equations,
 * then psuh forward in time and solve the elliptic equation
 * on the new slice.
 *
 * Revision 1.3  2004/05/06 15:26:29  e_gourgoulhon
 * No longer necessary to initialize khi and mu.
 *
 * Revision 1.2  2004/05/05 14:39:32  e_gourgoulhon
 * Added graphical outputs.
 *
 * Revision 1.1  2004/05/03 14:49:10  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "time_slice.h"
#include "param.h"
#include "graphique.h"
#include "utilitaires.h"
#include "proto.h"

const Tbl& monitor_scalar(const Scalar& uu, Tbl& resu) ;

void Tslice_dirac_max::evolve(double pdt, int nb_time_steps,
                              int niter_elliptic, double relax, 
                              int check_mod, int save_mod,
                              int method_poisson_vect, int nopause,  
                              const char* graph_device, bool verbose,
			      const Scalar* ener_euler,
			      const Vector* mom_euler, const Scalar* s_euler,
			      const Sym_tensor* strain_euler) {

    // Intermediate quantities
    // -----------------------
    const Map& map = nn().get_mp() ; 
    const Base_vect& triad = *(beta().get_triad()) ;
    assert( triad == map.get_bvect_spher() ) ;

    // For graphical outputs:
    int ngraph0 = 20 ;  // index of the first graphic device to be used
    int ngraph0_mon = 70 ;  // for monitoring global quantities
    int nz = map.get_mg()->get_nzone() ; 
    int nz_bound = nz - 2 ;
    int nt = map.get_mg()->get_nt(0) ;
    int np = map.get_mg()->get_np(0) ;
    int np2 = np+2 ;
    Scalar tmp(map) ; tmp.std_spectral_base() ;
    Base_val base_ref = tmp.get_spectral_base() ;
    Base_val base_pseudo = base_ref ;
    base_pseudo.mult_cost() ;
    base_pseudo.mult_x() ;
#ifndef NDEBUG
    for (int lz=1; lz<nz; lz++) {
	assert( map.get_mg()->get_np(lz) == np ) ;
	assert( map.get_mg()->get_nt(lz) == nt ) ;
    }
    assert (depth > 2) ; //## for the moment, a more flexible test should be put
    for (int it=0; it<depth; it++) {
	assert ( hh_evol.is_known( jtime - it ) ) ;
	assert ( hata_evol.is_known( jtime - it ) ) ;
	assert ( A_hata_evol.is_known( jtime - it ) ) ;
	assert ( B_hata_evol.is_known( jtime - it ) ) ;
	assert ( A_hh_evol.is_known( jtime - it ) ) ;
	assert ( B_hh_evol.is_known( jtime - it ) ) ;
    }	
#endif
    
    // Initialization of the TT fields
    //--------------------------------
    Sym_tensor_tt hij_tt( map, triad, ff ) ;
    hij_tt.set_A_tildeB( A_hh_evol[jtime], B_hh_evol[jtime] ) ;
    Sym_tensor_tt hijtt_old = hij_tt ;
    for (int i=1; i<=3; i++)
	for (int j=i; j<=3; j++)
	    if ( hijtt_old(i,j).get_etat() == ETATZERO )
		hijtt_old.set( i, j ).annule_hard() ;
    hijtt_old.annule(0, nz_bound) ;

    Sym_tensor_tt hata_tt( map, triad, ff ) ;
    hata_tt.set_A_tildeB( A_hata_evol[jtime], B_hata_evol[jtime] ) ;
    hata_tt.inc_dzpuis(2) ;
    Sym_tensor_tt hatatt_old = hata_tt ;
    for (int i=1; i<=3; i++)
	for (int j=i; j<=3; j++)
	    if ( hatatt_old(i,j).get_etat() == ETATZERO )
		hatatt_old.set( i, j ).annule_hard() ;
    hatatt_old.annule(0, nz_bound) ;

    // Declaration / initialization of mu and khi for hh and hata
    //-----------------------------------------------------------
    Evolution_std<Scalar> khi_hh_evol(depth) ;
    Evolution_std<Scalar> mu_hh_evol(depth) ;
    Evolution_std<Scalar> khi_a_evol(depth) ;
    Evolution_std<Scalar> mu_a_evol(depth) ;
    Sym_tensor_trans Tij(map, map.get_bvect_spher(), ff) ;
    for (int j=jtime-depth+1; j<=jtime; j++) {
	tmp = hij_tt(1,1) ;
	tmp.mult_r() ; tmp.mult_r() ;
	khi_hh_evol.update( tmp, j, the_time[j] ) ;
	mu_hh_evol.update( hij_tt.mu(), j, the_time[j] ) ;
	tmp = hata_tt(1,1) ;
	tmp.mult_r() ; tmp.mult_r() ;
	khi_a_evol.update( tmp, j, the_time[j] ) ;
	mu_a_evol.update( hata_tt.mu(), j, the_time[j] ) ;
    }

    double Rmax = map.val_r(nz-2, 1., 0., 0.) ; // outermost radius
    double ray_des = 1.25 * Rmax ; // for plots

    // Parameters for the evolution equations
    //---------------------------------------
    double an = 23./12. ; 
    double anm1 = -4./3. ; 
    double anm2 = 5./12. ; 

    Param par_A ;
    double l_val_A = 1./Rmax ;
    double l_der_A = 1. ;
    par_A.add_int(nz_bound, 0) ;
    par_A.add_int(2, 1) ; //matching of function and derivative
    par_A.add_int(0, 2) ;// no shift in l 
    par_A.add_int(2, 3) ; // only for l>=2
    par_A.add_double_mod(l_val_A, 0) ;
    par_A.add_double_mod(l_der_A, 1) ;
    Tbl* tmp_Acc = new Tbl(np2, nt) ;
    Tbl& Acc = *tmp_Acc ;
    Acc.annule_hard() ;
    par_A.add_tbl_mod(Acc) ;
    Param par_mat_A_hh ;

    Param par_B ;
    double l_val_B = 1./Rmax ;
    double l_der_B = 1. ; 
    par_B.add_int(nz_bound, 0) ;
    par_B.add_int(2, 1) ; //matching of function and derivative
    par_B.add_int(-1, 2) ;// shift in l for tilde{B}
    par_B.add_int(2, 3) ; // only for l>=2
    par_B.add_double_mod(l_val_B, 0) ;
    par_B.add_double_mod(l_der_B, 1) ;
    Tbl* tmp_Bcc = new Tbl(np2, nt) ;
    Tbl& Bcc = *tmp_Bcc ;
    Bcc.annule_hard() ;
    par_B.add_tbl_mod(Bcc) ;
    Param par_mat_B_hh ;

    Tbl xij_b(np2, nt) ;
    xij_b.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, B_hh_evol[jtime] , B_hata_evol[jtime], xij_b) ;
    Tbl xijm1_b(np2, nt) ;
    xijm1_b.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, B_hh_evol[jtime-1] , 
			   B_hata_evol[jtime-1], xijm1_b) ;
    Tbl xij_a(np2, nt) ;
    xij_a.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, A_hh_evol[jtime] , A_hata_evol[jtime], xij_a) ;
    Tbl xijm1_a(np2, nt) ;
    xijm1_a.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, A_hh_evol[jtime-1] , 
			   A_hata_evol[jtime-1], xijm1_a) ;

    // Parameters for the Dirac systems
    //---------------------------------

    Param par_bc_hh ;
    par_bc_hh.add_int(nz_bound) ;
    Tbl* cf_b_hh = new Tbl(10) ;
    cf_b_hh->annule_hard() ;
    cf_b_hh->set(0) = 11*Rmax + 12*pdt ; // mu
    cf_b_hh->set(1) = 6*Rmax*pdt ; // d mu / dr
    cf_b_hh->set(2) = 0 ; // X
    cf_b_hh->set(3) = 0 ; // d X / dr
    cf_b_hh->set(4) = 11*Rmax*Rmax + 18*Rmax*pdt ; // h^rr
    cf_b_hh->set(5) = 6*Rmax*Rmax*pdt ;  // d h^rr / dr
    cf_b_hh->set(6) = 0 ; //eta
    cf_b_hh->set(7) = 0 ; //d eta / dr
    cf_b_hh->set(8) = 0 ; //W
    cf_b_hh->set(9) = 0 ; //d W / dr
    par_bc_hh.add_tbl_mod(*cf_b_hh, 0) ;
    Tbl* kib_hh = new Tbl(np2, nt) ;
    Tbl& khib_hh = *kib_hh ;
    khib_hh.annule_hard() ;
    par_bc_hh.add_tbl_mod(khib_hh,1) ;
    Tbl* mb_hh = new Tbl(np2, nt) ;
    Tbl& mub_hh = *mb_hh ;
    mub_hh.annule_hard() ;
    par_bc_hh.add_tbl_mod(mub_hh, 2) ;

    Param par_mat_hh ;

    Tbl xij_mu_hh(np2, nt) ;
    xij_mu_hh.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, mu_hh_evol[jtime] , mu_a_evol[jtime], xij_mu_hh) ;
    Tbl xijm1_mu_hh(np2, nt) ;
    xijm1_mu_hh.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, mu_hh_evol[jtime-1] , mu_a_evol[jtime-1], 
			   xijm1_mu_hh) ;

    Tbl xij_ki_hh(np2, nt) ;
    xij_ki_hh.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, khi_hh_evol[jtime] , khi_a_evol[jtime], xij_ki_hh) ;
    Tbl xijm1_ki_hh(np2, nt) ;
    xijm1_ki_hh.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, khi_hh_evol[jtime-1] , khi_a_evol[jtime-1], 
			   xijm1_ki_hh) ;

    Param par_bc_hata ;
    par_bc_hata.add_int(nz_bound) ;
    Tbl* cf_b_hata = new Tbl(10) ;
    cf_b_hata->annule_hard() ;
    cf_b_hata->set(0) = 11*Rmax + 12*pdt ; // mu
    cf_b_hata->set(1) = 6*Rmax*pdt ; // d mu / dr
    cf_b_hata->set(2) = 0 ; // X
    cf_b_hata->set(3) = 0 ; // d X / dr
    cf_b_hata->set(4) = 11*Rmax*Rmax + 18*Rmax*pdt ; // h^rr
    cf_b_hata->set(5) =  6*Rmax*Rmax*pdt ;  // d h^rr / dr
    cf_b_hata->set(6) = 0 ; //eta
    cf_b_hata->set(7) = 0 ; //d eta / dr
    cf_b_hata->set(8) = 0 ; //W
    cf_b_hata->set(9) = 0 ; //d W / dr
    par_bc_hata.add_tbl_mod(*cf_b_hata, 0) ;
    Tbl* kib_hata = new Tbl(np2, nt) ;
    Tbl& khib_hata = *kib_hata ;
    khib_hata.annule_hard() ;
    par_bc_hata.add_tbl_mod(khib_hata,1) ;
    Tbl* mb_hata = new Tbl(np2, nt) ;
    Tbl& mub_hata = *mb_hata ;
    mub_hata.annule_hard() ;
    par_bc_hata.add_tbl_mod(mub_hata, 2) ;

    Param par_mat_hata ;

    Tbl xij_mu_a(np2, nt) ;
    xij_mu_a.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, mu_a_evol[jtime] , 
			   mu_a_evol.time_derive(jtime, 2), xij_mu_a) ;
    Tbl xijm1_mu_a(np2, nt) ;
    xijm1_mu_a.set_etat_qcq() ;
    tmp = ( mu_a_evol[jtime] - mu_a_evol[jtime-2] ) / (2.*pdt) ;
    initialize_outgoing_BC(nz_bound, mu_a_evol[jtime-1] , tmp, xijm1_mu_a) ;

    Tbl xij_ki_a(np2, nt) ;
    xij_ki_a.set_etat_qcq() ;
    initialize_outgoing_BC(nz_bound, khi_a_evol[jtime] , 
			   khi_a_evol.time_derive(jtime, 2), xij_ki_a) ;
    Tbl xijm1_ki_a(np2, nt) ;
    xijm1_ki_a.set_etat_qcq() ;
    tmp = ( khi_a_evol[jtime] - khi_a_evol[jtime-2] ) / (2.*pdt) ;
    initialize_outgoing_BC(nz_bound, khi_a_evol[jtime-1] , tmp, xijm1_ki_a) ;

    // Quantities at new time-step
    //----------------------------
    Scalar n_new(map) ; 
    Scalar psi_new(map) ; 
    Scalar npsi_new(map) ; 
    Vector beta_new(map, CON, triad) ; 
    Scalar A_hh_new(map) ; 
    Scalar B_hh_new(map) ; 
    Scalar A_hata_new(map) ; 
    Scalar B_hata_new(map) ; 
    
    // Successive values of various quantities:
    // ---------------------------------------
    Evolution_full<double> m_adm(adm_mass(), jtime, the_time[jtime]) ; 
    Evolution_full<double> test_ham_constr ; 
    Evolution_full<double> test_mom_constr_r ; 
    Evolution_full<double> test_mom_constr_t ; 
    Evolution_full<double> test_mom_constr_p ; 
    Evolution_full<Tbl> nn_monitor ;
    Evolution_full<Tbl> psi_monitor ;
    Evolution_full<Tbl> A_h_monitor ;
    Evolution_full<Tbl> B_h_monitor ;
    Evolution_full<Tbl> trh_monitor ;
    Evolution_full<Tbl> beta_monitor_maxabs ;
    Evolution_full<Tbl> hh_monitor_central ;
    Evolution_full<Tbl> hh_monitor_maxabs ;
    Evolution_full<Tbl> hata_monitor_central ;
    Evolution_full<Tbl> hata_monitor_maxabs ;
    Evolution_full<Tbl> check_evol ;
    Tbl select_scalar(6) ; 
    Tbl select_tens(6) ;

    Vector zero_vec( map, CON, map.get_bvect_spher() ) ;
    zero_vec.set_etat_zero() ;
    const Vector& hat_S = ( mom_euler == 0x0 ? zero_vec : *mom_euler ) ;
    Scalar lapB(map) ;
    Scalar lapBm1 = source_B_hata_evol[jtime-1] ;
    Scalar lapBm2 = source_B_hata_evol[jtime-2] ;

    // Evolution loop
    // --------------

    for (int jt = 0; jt < nb_time_steps; jt++) {
    
        double ttime = the_time[jtime] ; 
	k_dd() ;
            
        if (jt%check_mod == 0) {
	  cout << 
	    "==============================================================\n"
	       << "  step: " << jtime << "   time = " << the_time[jtime] << endl  
	       << " ADM mass : " << adm_mass() 
	       << ", Log of central lapse: " << log(nn().val_grid_point(0,0,0,0)) << endl 
	       << "==============================================================\n" ;
	  
	  // Monitoring
	  // ---------- 
	  m_adm.update(adm_mass(), jtime, the_time[jtime]) ;
	  if (jt > 0) des_evol(m_adm, "ADM mass", "Variation of ADM mass", 
			       ngraph0_mon, graph_device) ;          
	  
	  
	  nn_monitor.update(monitor_scalar(nn(), select_scalar), 
                            jtime, the_time[jtime]) ; 
	  
	  psi_monitor.update(monitor_scalar(psi(), select_scalar), 
			     jtime, the_time[jtime]) ; 
	  
	  A_h_monitor.update(monitor_scalar(A_hh(), select_scalar), 
			     jtime, the_time[jtime]) ; 
	  
	  B_h_monitor.update(monitor_scalar(B_hh(), select_scalar), 
			     jtime, the_time[jtime]) ; 
        
	  trh_monitor.update(monitor_scalar(trh(), select_scalar), 
			     jtime, the_time[jtime]) ; 
	  
	  beta_monitor_maxabs.update(maxabs_all_domains(beta(), -1, 0x0, cout, verbose), 
				     jtime, the_time[jtime]) ; 
        
	  hh_monitor_central.update(central_value(hh()), 
                                    jtime, the_time[jtime]) ; 
        
	  hh_monitor_maxabs.update(maxabs_all_domains(hh(), -1, 0x0, cout, verbose), 
				   jtime, the_time[jtime]) ; 
        
	  hata_monitor_central.update(central_value(hata()), 
				      jtime, the_time[jtime]) ; 
        
	  hata_monitor_maxabs.update(maxabs_all_domains(hata(), -1, 0x0, cout, verbose), 
				     jtime, the_time[jtime]) ; 
        

	  int jt_graph = jt / check_mod ; 
            
	  Tbl tham = check_hamiltonian_constraint(0x0, cout, verbose) ; 
	  double max_error = tham(0,0) ; 
	  for (int l=1; l<nz-1; l++) {    // all domains but the last one
	    double xx = fabs(tham(0,l)) ;  
	    if (xx > max_error) max_error = xx ; 
	  }
	  test_ham_constr.update(max_error, jt_graph, the_time[jtime]) ; 
	  if (jt > 0) des_evol(test_ham_constr, "Absolute error", 
			       "Check of Hamiltonian constraint", 
			       ngraph0_mon+1, graph_device) ; 
	  
	  Tbl tmom = check_momentum_constraint(0x0, cout, verbose) ; 
            max_error = tmom(0,0) ;
            for (int l=1; l<nz-1; l++) {    // all domains but the last one
                double xx = fabs(tmom(0,l)) ;  
                if (xx > max_error) max_error = xx ; 
            }
            test_mom_constr_r.update(max_error, jt_graph, the_time[jtime]) ; 
            if (jt > 0) des_evol(test_mom_constr_r, "Absolute error", 
                "Check of momentum constraint (r comp.)", ngraph0_mon+2, 
                graph_device) ; 

            max_error = tmom(1,0) ;
            for (int l=1; l<nz-1; l++) {    // all domains but the last one
                double xx = fabs(tmom(1,l)) ;  
                if (xx > max_error) max_error = xx ; 
            }
            test_mom_constr_t.update(max_error, jt_graph, the_time[jtime]) ; 
            if (jt > 0) des_evol(test_mom_constr_t, "Absolute error", 
                "Check of momentum constraint (\\gh comp.)", ngraph0_mon+3,
                 graph_device) ; 

            max_error = tmom(2,0) ;
            for (int l=1; l<nz-1; l++) {    // all domains but the last one
                double xx = fabs(tmom(2,l)) ;  
                if (xx > max_error) max_error = xx ; 
            }
            test_mom_constr_p.update(max_error, jt_graph, the_time[jtime]) ; 
            if (jt > 0) des_evol(test_mom_constr_p, "Absolute error", 
                "Check of momentum constraint (\\gf comp.)", ngraph0_mon+4, 
                graph_device) ; 
               
	    if (jt>2) {
	      Tbl tevol = check_dynamical_equations(0x0, 0x0, cout, verbose) ;
	    Tbl evol_check(6) ; evol_check.set_etat_qcq() ;
	    for (int i=1; i<=3; i++) 
		for(int j=1; j<=i; j++) {
		    max_error = tevol(i, j, 0) ;
		    for (int l=1; l<nz-1; l++) {
			double xx = fabs(tevol(i,j,l)) ;
			if (xx > max_error) max_error = xx ;
		    }
		    evol_check.set(i) = max_error ;
		}
	    check_evol.update(evol_check,  jtime, the_time[jtime]) ;
	    }
        }

        if (jt%save_mod == 0) { 
            m_adm.save("adm_mass.d") ; 
            nn_monitor.save("nn_monitor.d") ;
            psi_monitor.save("psi_monitor.d") ;
            A_h_monitor.save("potA_monitor.d") ;
            B_h_monitor.save("potB_monitor.d") ;
            trh_monitor.save("trh_monitor.d") ;
            beta_monitor_maxabs.save("beta_monitor_maxabs.d") ; 
            hh_monitor_central.save("hh_monitor_central.d") ; 
            hh_monitor_maxabs.save("hh_monitor_maxabs.d") ; 
            hata_monitor_central.save("hata_monitor_central.d") ; 
            hata_monitor_maxabs.save("hata_monitor_maxabs.d") ; 
            test_ham_constr.save("test_ham_constr.d") ; 
            test_mom_constr_r.save("test_mom_constr_r.d") ; 
            test_mom_constr_t.save("test_mom_constr_t.d") ; 
            test_mom_constr_p.save("test_mom_constr_p.d") ; 
	    check_evol.save("evol_equations.d") ;
            
            save("sigma") ;
            
        }


        // Resolution of hyperbolic equations
        // ----------------------------------
	compute_sources(strain_euler) ;

	A_hata_new = A_hata_evol[jtime] 
	  + pdt*( an*source_A_hata_evol[jtime] + anm1*source_A_hata_evol[jtime-1]
		  + anm2*source_A_hata_evol[jtime-2] ) ;
	B_hata_new = B_hata_evol[jtime] 
	  + pdt*( an*source_B_hata_evol[jtime] + anm1*source_B_hata_evol[jtime-1]
		  + anm2*source_B_hata_evol[jtime-2] ) ;

	A_hh_new = A_hh_evol[jtime] 
	  + pdt*( an*source_A_hh_evol[jtime] + anm1*source_A_hh_evol[jtime-1]
		  + anm2*source_A_hh_evol[jtime-2] ) ;

	B_hh_new = B_hh_evol[jtime] 
	  + pdt*( an*source_B_hh_evol[jtime] + anm1*source_B_hh_evol[jtime-1]
		  + anm2*source_B_hh_evol[jtime-2] ) ;
	
	Scalar bc_A = -2.*A_hata_new ;
	bc_A.set_spectral_va().ylm() ;
	evolve_outgoing_BC(pdt, nz_bound, A_hh_evol[jtime], bc_A, xij_a, xijm1_a, 
	 		   Acc, 0) ;
	A_hh_new.match_tau(par_A, &par_mat_A_hh) ;
        
	Scalar bc_B = -2.*B_hata_new ;
	bc_B.set_spectral_va().ylm() ;
	evolve_outgoing_BC(pdt, nz_bound, B_hh_evol[jtime], bc_B, xij_b, xijm1_b, 
	  		   Bcc, -1) ;
	B_hh_new.match_tau(par_B, &par_mat_B_hh) ;
	
        // Boundary conditions for hh and hata
	//------------------------------------
  	Scalar sbcmu = (18*mu_hh_evol[jtime] - 9*mu_hh_evol[jtime-1] 
			+ 2*mu_hh_evol[jtime-2]) / (6*pdt) ;
	if (sbcmu.get_etat() == ETATZERO) {
	    sbcmu.annule_hard() ;
	    sbcmu.set_spectral_base(base_pseudo) ;
	}
   	sbcmu.set_spectral_va().ylm() ;
	tmp = mu_hh_evol[jtime] ;
	if (tmp.get_etat() == ETATZERO) {
	    tmp.annule_hard() ;
	    tmp.set_spectral_base(base_pseudo) ;
	}
   	tmp.set_spectral_va().ylm() ;
  	evolve_outgoing_BC(pdt, nz_bound, tmp, sbcmu, xij_mu_hh, xijm1_mu_hh, 
			   mub_hh, 0) ;
 	mub_hh *= 6*pdt ;

  	Scalar sbckhi = (18*khi_hh_evol[jtime] - 9*khi_hh_evol[jtime-1] 
			 + 2*khi_hh_evol[jtime-2]) / (6*pdt) ;
	if (sbckhi.get_etat() == ETATZERO) {
	    sbckhi.annule_hard() ;
	    sbckhi.set_spectral_base(base_ref) ;
	}
   	sbckhi.set_spectral_va().ylm() ;
	tmp = khi_hh_evol[jtime] ;
	if (tmp.get_etat() == ETATZERO) {
	    tmp.annule_hard() ;
	    tmp.set_spectral_base(base_ref) ;
	}
   	tmp.set_spectral_va().ylm() ;
  	evolve_outgoing_BC(pdt, nz_bound, tmp, sbckhi, xij_ki_hh, xijm1_ki_hh, 
			   khib_hh, 0) ;
 	khib_hh *= 6*pdt ;

  	sbcmu = (18*mu_a_evol[jtime] - 9*mu_a_evol[jtime-1] 
			+ 2*mu_a_evol[jtime-2]) / (6*pdt) ;
	if (sbcmu.get_etat() == ETATZERO) {
	    sbcmu.annule_hard() ;
	    sbcmu.set_spectral_base(base_pseudo) ;
	}
    	sbcmu.set_spectral_va().ylm() ;
	tmp = mu_a_evol[jtime] ;
	if (tmp.get_etat() == ETATZERO) {
	    tmp.annule_hard() ;
	    tmp.set_spectral_base(base_pseudo) ;
	}
   	tmp.set_spectral_va().ylm() ;
  	evolve_outgoing_BC(pdt, nz_bound, tmp, sbcmu, xij_mu_a, xijm1_mu_a, 
 			   mub_hata, 0) ;
 	mub_hata *= 6*pdt ;

  	sbckhi = (18*khi_a_evol[jtime] - 9*khi_a_evol[jtime-1] 
			 + 2*khi_a_evol[jtime-2]) / (6*pdt) ;
	if (sbckhi.get_etat() == ETATZERO) {
	    sbckhi.annule_hard() ;
	    sbckhi.set_spectral_base(base_ref) ;
	}
   	sbckhi.set_spectral_va().ylm() ;
	tmp = khi_a_evol[jtime] ;
	if (tmp.get_etat() == ETATZERO) {
	    tmp.annule_hard() ;
	    tmp.set_spectral_base(base_ref) ;
	}
   	tmp.set_spectral_va().ylm() ;
  	evolve_outgoing_BC(pdt, nz_bound, tmp, sbckhi, xij_ki_a, xijm1_ki_a, 
			   khib_hata, 0) ;
 	khib_hata *= 6*pdt ;

        // Advance in time
        // ---------------
        
        jtime++ ; 
        ttime += pdt ; 
        the_time.update(ttime, jtime, ttime) ; 

	// Setting As and Bs for h^{ij} and \hat{A}^{ij}
        set_AB_hata(A_hata_new, B_hata_new) ;
        set_AB_hh(A_hh_new, B_hh_new) ;

	hij_tt.set_A_tildeB( A_hh_new, B_hh_new, &par_bc_hh, &par_mat_hh ) ;
	for (int i=1; i<=3; i++)
	    for (int j=i; j<=3; j++) 
		for (int l=nz_bound+1; l<nz; l++)
		    hij_tt.set(i,j).set_domain(l) = hijtt_old(i,j).domain(l) ;
	hata_tt.set_A_tildeB( A_hata_new, B_hata_new, &par_bc_hata, &par_mat_hata ) ;
	for (int i=1; i<=3; i++)
	    for (int j=i; j<=3; j++) {
		for (int l=nz_bound+1; l<nz; l++)
		    hata_tt.set(i,j).set_domain(l) = hatatt_old(i,j).domain(l) ;
		hata_tt.set(i,j).set_dzpuis(2) ;
	    }

	// Computation of h^{ij} at new time-step
	hh_det_one(hij_tt, &par_mat_hh) ;

        // Reset of derived quantities
        del_deriv() ;        

	// Update of khi's and mu's
	//-------------------------
	tmp = hij_tt( 1, 1 ) ;
	tmp.mult_r() ; tmp.mult_r() ;
	khi_hh_evol.update( tmp, jtime, the_time[jtime] ) ;
	mu_hh_evol.update( hij_tt.mu(), jtime, the_time[jtime] ) ;
	tmp = hata_tt( 1, 1 ) ;
	tmp.mult_r() ; tmp.mult_r() ;
	khi_a_evol.update( tmp, jtime, the_time[jtime] ) ;
	mu_a_evol.update( hata_tt.mu(), jtime, the_time[jtime] ) ;
	
        // Resolution of elliptic equations
        // --------------------------------
        psi_evol.update(psi_evol[jtime-1], jtime, ttime) ; 
        
	// \hat{A}^{ij} is computed at the new time-step
	compute_X_from_momentum_constraint(hat_S, hata_tt, niter_elliptic) ;
	
	// Iteration on the conformal factor
        for (int k = 0; k < niter_elliptic; k++) {
    
            psi_new = solve_psi(ener_euler) ; 
            psi_new = relax * psi_new + (1.-relax) * psi() ; 
	    set_psi_del_npsi(psi_new) ;   
        }

        set_npsi_del_n(npsi_evol[jtime-1]) ; 

	// Iteration on N*Psi  ## play with the number of iterations...
        npsi_evol.update(psi_evol[jtime-1], jtime, ttime) ; 
	for (int k = 0; k < niter_elliptic; k++) {

            npsi_new = solve_npsi( ener_euler, s_euler ) ; 
            npsi_new = relax * npsi_new + (1.-relax) * npsi() ; 
	    set_npsi_del_n(npsi_new) ; 
	}

	// Iteration on beta ## play with the number of iterations...
        beta_evol.update(beta_evol[jtime-1], jtime, ttime) ;             
	for (int k = 0; k < niter_elliptic; k++) {

            beta_new = solve_beta(method_poisson_vect) ; 
            beta_new = relax * beta_new + (1.-relax) * beta() ;
	    beta_evol.update(beta_new, jtime, ttime) ;             
	}    
                
        des_meridian(vec_X()(1), 0., ray_des, "\\gb\\ur\\d", ngraph0+6,
                     graph_device) ; 
        des_meridian(vec_X()(2), 0., ray_des, "\\gb\\u\\gh\\d", ngraph0+7,
                     graph_device) ; 
        des_meridian(vec_X()(3), 0., ray_des, "\\gb\\u\\gf\\d", ngraph0+8,
                     graph_device) ; 
	tmp = A_hh() ;
	tmp.set_spectral_va().ylm_i() ;
        des_meridian(tmp, 0., ray_des, "A\\dh", ngraph0+9,
                     graph_device) ; 
	tmp = B_hh_new;
	tmp.set_spectral_va().ylm_i() ;
        des_meridian(tmp, 0., ray_des, "B\\dh", ngraph0+10,
                     graph_device) ;
        des_meridian(trh(), 0., ray_des, "tr h", ngraph0+11,
                     graph_device) ; 
        des_meridian(hh()(1,1), 0., ray_des, "h\\urr\\d", ngraph0+12,
                     graph_device) ; 
        des_meridian(hh()(2,3), 0., ray_des, "h\\u\\gh\\gf\\d", ngraph0+13,
                     graph_device) ; 
        des_meridian(hh()(3,3), 0., ray_des, "h\\u\\gf\\gf\\d", ngraph0+14,
                     graph_device) ; 
                
        arrete(nopause) ; 
    }

    par_A.clean_all() ;
    par_B.clean_all() ;
    par_mat_A_hh.clean_all() ;
    par_mat_B_hh.clean_all() ;

    par_bc_hh.clean_all() ;
    par_mat_hh.clean_all() ;

    par_bc_hata.clean_all() ;
    par_mat_hata.clean_all() ;
} 


//***************************************************************************

const Tbl& monitor_scalar(const Scalar& uu, Tbl& resu) {

    assert( resu.get_ndim() == 1) ; 
    assert( resu.get_taille() >= 6) ;
    
    resu.set_etat_qcq() ; 
    
    resu.set(0) = uu.val_grid_point(0,0,0,0) ; 
    resu.set(1) = max(max(uu)) ; 
    resu.set(2) = min(min(uu)) ; 
    
    const Mg3d& mg = *(uu.get_mp().get_mg()) ; 
    
    int nz = mg.get_nzone() ;
    int nzm1 = nz - 1 ;  
    int nr = mg.get_nr(nzm1) ; 
    int nt = mg.get_nt(nzm1) ; 
    int np = mg.get_np(nzm1) ; 
    
    resu.set(3) = uu.val_grid_point(nzm1, 0, 0, nr-1) ; 
    resu.set(4) = uu.val_grid_point(nzm1, 0, nt-1, nr-1) ; 
    resu.set(5) = uu.val_grid_point(nzm1, np/2, nt-1, nr-1) ; 
    
    return resu ;      
}
