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

    Scalar n_new(map) ; 
    Scalar q_new(map) ; 
    Vector beta_new(map, CON, triad) ; 
    Scalar khi_new(map) ; 
    Scalar mu_new(map) ; 
    Scalar n_backup(map) ; 
    Scalar q_backup(map) ; 
    Vector beta_backup(map, CON, triad) ; 

    // Recovering the values of khi and mu on previous slices
    // (in order to compute time derivatives)
    // ------------------------------------------------------

    if (!(khi_evol.is_known(jtime-1))) {
    
        khi_new = khi() ; 

        for (int j = jtime-depth+1 ; j <= jtime; j++) {
            khi_evol.downdate(j) ;  // cleaning; to be set up below
        }
        for (int j = jtime-depth+1 ; j <= jtime; j++) {
            // Time derivative = zero : 
            khi_evol.update(khi_new, j, the_time[j]) ;
            
            //## Alternative : time derivative consistant with that of h^{ij} : 
            // khi_evol.update((hh_evol[j]).transverse(ff).tt_part().khi(),
            //                j, the_time[j]) ;
       }
    }

    if (!(mu_evol.is_known(jtime-1))) {
    
        mu_new = mu() ; 

        for (int j = jtime-depth+1 ; j <= jtime; j++) {
            mu_evol.downdate(j) ; // cleaning; to be set up below
        }
        for (int j = jtime-depth+1 ; j <= jtime; j++) {
            // Time derivative = zero : 
            mu_evol.update(mu_new, j, the_time[j]) ;
            
            //## Alternative : time derivative consistant with that of h^{ij} : 
            // mu_evol.update((hh_evol[j]).transverse(ff).tt_part().mu(),
            //                j, the_time[j]) ;
       }
    }



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

        des_meridian(hh()(1,1), 0., 5., "h\\urr\\d", 30) ; 
        des_meridian(hh()(2,3), 0., 5., "h\\u\\gh\\gf\\d", 31) ; 
        des_meridian(hh()(3,3), 0., 5., "h\\u\\gf\\gf\\d", 32) ; 
        
        arrete() ; 

    }

    par_khi.clean_all() ; 
    par_mu.clean_all() ; 
} 
