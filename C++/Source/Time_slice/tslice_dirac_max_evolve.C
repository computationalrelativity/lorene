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
#include "metric.h"
#include "evolution.h"
#include "param.h"
#include "graphique.h"
#include "utilitaires.h"

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
    Scalar n_backup(map) ; 
    Scalar q_backup(map) ; 
    Vector beta_backup(map, CON, triad) ; 

    // Evolution loop
    // --------------

    for (int jt = 0; jt < nb_time_steps; jt++) {
    
        double ttime = the_time[jtime] ; 
    
        cout << 
        "==============================================================\n"
        << "  step: " << jtime << "   time = " << the_time[jtime] << endl  
        << "==============================================================\n" ;
    
        cout << *this << endl ; 
        
        // Resolution of elliptic equations
        // --------------------------------
        
        n_backup = nn() ; 
        q_backup = qq() ; 
        beta_backup = beta() ; 
        
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
                
        // Restoring values at current time step (before computing the source 
        //  for h^{ij})
        
        n_evol.update(n_backup, jtime, ttime) ; 
        set_qq_del_psi(q_backup) ; 
        beta_evol.update(beta_backup, jtime, ttime) ;             
        
        // Resolution of hyperbolic equations
        // ----------------------------------
        
        solve_hij(par_khi, par_mu, khi_new, mu_new) ;
        
        // Advance in time
        // ---------------
        
        jtime++ ; 
        ttime += pdt ; 
        the_time.update(ttime, jtime, ttime) ; 
                
        n_evol.update(n_new, jtime, ttime) ; 
        set_qq_del_psi(q_new) ; 
        beta_evol.update(beta_new, jtime, ttime) ;             
        
        // Updates of khi_evol, mu_evol, trh_evol and hh_evol:
        set_khi_mu(khi_new, mu_new) ;           

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
