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

const Tbl& monitor_scalar(const Scalar& uu, Tbl& resu) ;

void Tslice_dirac_max::evolve(double pdt, int nb_time_steps,
                              int niter_elliptic, 
                              double relax) {

    // Parameters for the d'Alembert equations
    // ----------------------------------------
    int bc = 2 ;    // type of boundary condition : 2 = Bayliss & Turkel outgoing wave
 
    Param par_khi ; 
    par_khi.add_double(pdt) ; 
    par_khi.add_int(bc) ; 
    int *workflag_khi = new int(0) ; // working flag 
    par_khi.add_int_mod(*workflag_khi) ; 
    
    Param par_mu ; 
    par_mu.add_double(pdt) ; 
    par_mu.add_int(bc) ; 
    int *workflag_mu = new int(0) ; // working flag 
    par_mu.add_int_mod(*workflag_mu) ; 

    // Intermediate quantities
    // -----------------------
    const Map& map = nn().get_mp() ; 
    const Base_vect& triad = *(beta().get_triad()) ;
    int nz = map.get_mg()->get_nzone() ; 
    double ray_des = 1.25 * map.val_r(nz-2, 1., 0., 0.) ; // outermost radius
                                                          // for plots
    Scalar n_new(map) ; 
    Scalar q_new(map) ; 
    Vector beta_new(map, CON, triad) ; 
    Scalar khi_new(map) ; 
    Scalar mu_new(map) ; 
    
    // Successive values of various quantities:
    // ---------------------------------------
    Evolution_full<double> m_adm(adm_mass(), jtime, the_time[jtime]) ; 
    Evolution_full<double> test_ham_constr ; 
    Evolution_full<double> test_mom_constr_r ; 
    Evolution_full<double> test_mom_constr_t ; 
    Evolution_full<double> test_mom_constr_p ; 
    Evolution_full<Tbl> nn_monitor ;
    Evolution_full<Tbl> psi_monitor ;
    Evolution_full<Tbl> trh_monitor ;
    Evolution_full<Tbl> beta_monitor_maxabs ;
    Evolution_full<Tbl> hh_monitor_central ;
    Evolution_full<Tbl> hh_monitor_maxabs ;
    Evolution_full<Tbl> aa_monitor_central ;
    Evolution_full<Tbl> aa_monitor_maxabs ;
    Tbl select_scalar(6) ; 
    Tbl select_tens(6) ;
        
    // Evolution loop
    // --------------

    for (int jt = 0; jt < nb_time_steps; jt++) {
    
        double ttime = the_time[jtime] ; 
    
        cout << 
        "==============================================================\n"
        << "  step: " << jtime << "   time = " << the_time[jtime] << endl  
        << "==============================================================\n" ;
    
        cout << *this << endl ; 
        
        // Monitoring
        // ---------- 
        cout << "ADM mass : " << adm_mass() << endl ; 
        m_adm.update(adm_mass(), jtime, the_time[jtime]) ;
        if (jt > 0) des_evol(m_adm, "ADM mass", "Variation of ADM mass", 80) ;          
        
        
        nn_monitor.update(monitor_scalar(nn(), select_scalar), 
                            jtime, the_time[jtime]) ; 
        
        psi_monitor.update(monitor_scalar(psi(), select_scalar), 
                            jtime, the_time[jtime]) ; 
        
        trh_monitor.update(monitor_scalar(trh(), select_scalar), 
                            jtime, the_time[jtime]) ; 
        
        beta_monitor_maxabs.update(maxabs_all_domains(beta()), 
                                    jtime, the_time[jtime]) ; 
        
        hh_monitor_central.update(central_value(hh()), 
                                    jtime, the_time[jtime]) ; 
        
        hh_monitor_maxabs.update(maxabs_all_domains(hh()), 
                                    jtime, the_time[jtime]) ; 
        
        aa_monitor_central.update(central_value(aa()), 
                                    jtime, the_time[jtime]) ; 
        
        aa_monitor_maxabs.update(maxabs_all_domains(aa()), 
                                    jtime, the_time[jtime]) ; 
        
        int check_mod = 2 ; 
        if (jt%check_mod == 0) {

            int jt_graph = jt / check_mod ; 
            
            Tbl tham = check_hamiltonian_constraint() ; 
            double max_error = tham(0,0) ; 
            for (int l=1; l<nz-1; l++) {    // all domains but the last one
                double xx = fabs(tham(0,l)) ;  
                if (xx > max_error) max_error = xx ; 
            }
            test_ham_constr.update(max_error, jt_graph, the_time[jtime]) ; 
            if (jt > 0) des_evol(test_ham_constr, "Absolute error", 
                "Check of Hamiltonian constraint", 81) ; 

            Tbl tmom = check_momentum_constraint() ; 
            max_error = tmom(0,0) ;
            for (int l=1; l<nz-1; l++) {    // all domains but the last one
                double xx = fabs(tmom(0,l)) ;  
                if (xx > max_error) max_error = xx ; 
            }
            test_mom_constr_r.update(max_error, jt_graph, the_time[jtime]) ; 
            if (jt > 0) des_evol(test_mom_constr_r, "Absolute error", 
                "Check of momentum constraint (r comp.)", 82) ; 

            max_error = tmom(1,0) ;
            for (int l=1; l<nz-1; l++) {    // all domains but the last one
                double xx = fabs(tmom(1,l)) ;  
                if (xx > max_error) max_error = xx ; 
            }
            test_mom_constr_t.update(max_error, jt_graph, the_time[jtime]) ; 
            if (jt > 0) des_evol(test_mom_constr_t, "Absolute error", 
                "Check of momentum constraint (\\gh comp.)", 83) ; 

            max_error = tmom(2,0) ;
            for (int l=1; l<nz-1; l++) {    // all domains but the last one
                double xx = fabs(tmom(2,l)) ;  
                if (xx > max_error) max_error = xx ; 
            }
            test_mom_constr_p.update(max_error, jt_graph, the_time[jtime]) ; 
            if (jt > 0) des_evol(test_mom_constr_p, "Absolute error", 
                "Check of momentum constraint (\\gf comp.)", 84) ; 
                
        }

        int save_mod = 10 ; 
        if (jt%save_mod == 0) { 
            m_adm.save("adm_mass.d") ; 
            nn_monitor.save("nn_monitor.d") ;
            psi_monitor.save("psi_monitor.d") ;
            trh_monitor.save("trh_monitor.d") ;
            beta_monitor_maxabs.save("beta_monitor_maxabs.d") ; 
            hh_monitor_central.save("hh_monitor_central.d") ; 
            hh_monitor_maxabs.save("hh_monitor_maxabs.d") ; 
            aa_monitor_central.save("aa_monitor_central.d") ; 
            aa_monitor_maxabs.save("aa_monitor_maxabs.d") ; 
            test_ham_constr.save("test_ham_constr.d") ; 
            test_mom_constr_r.save("test_mom_constr_r.d") ; 
            test_mom_constr_t.save("test_mom_constr_t.d") ; 
            test_mom_constr_p.save("test_mom_constr_p.d") ; 
        }


        // Resolution of hyperbolic equations
        // ----------------------------------
        
        solve_hij(par_khi, par_mu, khi_new, mu_new) ;
        
        // Advance in time
        // ---------------
        
        jtime++ ; 
        ttime += pdt ; 
        the_time.update(ttime, jtime, ttime) ; 
                
        // Setting khi_evol, mu_evol, trh_evol and hh_evol at the new time:
        set_khi_mu(khi_new, mu_new) ;    
        
        // Reset of derived quantities
        del_deriv() ;        

        // Resolution of elliptic equations
        // --------------------------------
        
        // N, Q and beta at the new time are initialized by their
        // values at previous time:
        
        n_evol.update(n_evol[jtime-1], jtime, ttime) ; 
        set_qq_del_psi(qq_evol[jtime-1]) ; 
        beta_evol.update(beta_evol[jtime-1], jtime, ttime) ;             
        
        for (int k = 0; k < niter_elliptic; k++) {
    
            des_meridian(aa()(1,1), 0., ray_des, "A\\urr\\d", 20) ; 
            des_meridian(aa()(2,3), 0., ray_des, "A\\u\\gh\\gf\\d", 21) ; 
            des_meridian(aa()(3,3), 0., ray_des, "A\\u\\gf\\gf\\d", 22) ; 

            n_new = solve_n() ; 
            q_new = solve_q() ; 
            beta_new = solve_beta() ; 
    
            n_new = relax * n_new + (1.-relax) * nn() ; 
            q_new = relax * q_new + (1.-relax) * qq() ; 
            beta_new = relax * beta_new + (1.-relax) * beta() ;
        
            maxabs(n_new - nn(), "Difference between N_new and N") ;     
            maxabs(q_new - qq(), "Difference between Q_new and Q") ;     
            maxabs(beta_new - beta(), "Difference between beta_new and beta") ;     
        
            n_evol.update(n_new, jtime, ttime) ; 
            set_qq_del_psi(q_new) ; 
            beta_evol.update(beta_new, jtime, ttime) ;             

        }
                
        des_meridian(beta()(1), 0., ray_des, "\\gb\\ur\\d", 10) ; 
        des_meridian(beta()(2), 0., ray_des, "\\gb\\u\\gh\\d", 11) ; 
        des_meridian(beta()(3), 0., ray_des, "\\gb\\u\\gf\\d", 12) ; 
        des_meridian(hh()(1,1), 0., ray_des, "h\\urr\\d", 13) ; 
        des_meridian(hh()(2,3), 0., ray_des, "h\\u\\gh\\gf\\d", 14) ; 
        des_meridian(hh()(3,3), 0., ray_des, "h\\u\\gf\\gf\\d", 15) ; 
        
                
        // arrete() ; 

    }

    par_khi.clean_all() ; 
    par_mu.clean_all() ; 
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



