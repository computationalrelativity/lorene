/*
 *  Virtual methods of class Time_slice and derived classes to
 *  compute the ADM mass
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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

char tslice_adm_mass_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2010/10/20 07:58:09  j_novak
 * Better implementation of the explicit time-integration. Not fully-tested yet.
 *
 * Revision 1.4  2008/12/04 19:36:40  j_novak
 * Removed old tests.
 *
 * Revision 1.3  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.2  2004/05/10 09:11:20  e_gourgoulhon
 * Corrected bug in Time_slice::adm_mass().
 * Reorganized outputs.
 *
 * Revision 1.1  2004/05/09 20:56:29  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */


// C headers
#include <math.h>

// Lorene headers
#include "time_slice.h"

//--------------------
// Time_slice version 
//--------------------

double Time_slice::adm_mass() const {

    if ( !(adm_mass_evol).is_known(jtime) ) {  // a new computation is necessary
    
        const Map& mp = gam_dd().get_mp() ;
        Metric_flat ff(mp, *(gam_dd().get_triad())) ;
        int nz = mp.get_mg()->get_nzone() ; 
        Tbl* tmass = new Tbl(nz) ; 
        tmass->set_etat_qcq() ; 
    
        Vector ww = gam_dd().derive_con(ff).trace(1,2).up(0,ff) 
                    - gam_dd().trace(ff).derive_con(ff) ; 

        for (int l=0; l<nz; l++) {
            double radius = mp.val_r(l, 1., 0., 0.) ;
            tmass->set(l) = ww.flux(radius, ff) / (16.* M_PI) ; 
        }
        
        adm_mass_evol.update(*tmass, jtime, the_time[jtime]) ; 
        
        delete tmass ;  
        
        cout << "Time_slice::adm_mass : " << adm_mass_evol[jtime] << endl ; 
    }
  
    const Tbl& tadm = adm_mass_evol[jtime] ; 
    return tadm(tadm.get_taille()-1) ; 
}


//--------------------------
// Time_slice_conf version 
//--------------------------

double Time_slice_conf::adm_mass() const {

    if ( !(adm_mass_evol).is_known(jtime) ) {  // a new computation is necessary
    
        const Map& mp = psi().get_mp() ;
        int nz = mp.get_mg()->get_nzone() ; 
        Tbl* tmass = new Tbl(nz) ; 
        tmass->set_etat_qcq() ; 
    
        Vector ww = psi().derive_con(ff) 
                    + 0.125* ( hdirac() - (hh().trace(ff)).derive_con(ff) )  ; 

        for (int l=0; l<nz; l++) {
            double radius = mp.val_r(l, 1., 0., 0.) ;
            tmass->set(l) = - ww.flux(radius, ff) / (2.* M_PI) ; 
        }
        
        adm_mass_evol.update(*tmass, jtime, the_time[jtime]) ; 
        
        delete tmass ;  
    
#ifndef NDEBUG
        cout << "Time_slice_conf::adm_mass : " << adm_mass_evol[jtime] << endl ; 
#endif    
    }
  
    const Tbl& tadm = adm_mass_evol[jtime] ; 
    return tadm(tadm.get_taille()-1) ; 
}


//--------------------------
// Tslice_dirac_max version 
//--------------------------

double Tslice_dirac_max::adm_mass() const {

    if ( !(adm_mass_evol).is_known(jtime) ) {  // a new computation is necessary
    
        const Map& mp = psi().get_mp() ;
        int nz = mp.get_mg()->get_nzone() ; 
        Tbl* tmass = new Tbl(nz) ; 
        tmass->set_etat_qcq() ; 
    
        Vector ww = psi().derive_con(ff) 
                    - 0.125* (hh().trace(ff)).derive_con(ff) ;
                    // trh() is not used since it has dzpuis = 4 

        for (int l=0; l<nz; l++) {
            double radius = mp.val_r(l, 1., 0., 0.) ;
            tmass->set(l) = - ww.flux(radius, ff) / (2.* M_PI) ; 
        }
        
        adm_mass_evol.update(*tmass, jtime, the_time[jtime]) ; 
        
        delete tmass ;  
#ifndef NDEBUG
        cout << "Tslice_dirac_max::adm_mass : " << adm_mass_evol[jtime] 
             << endl ; 
#endif

    }
  
    const Tbl& tadm = adm_mass_evol[jtime] ; 
    return tadm(tadm.get_taille()-1) ; 
}




