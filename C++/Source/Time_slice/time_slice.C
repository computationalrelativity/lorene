/*
 *  Methods of class Time_slice
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon, Jose Luis Jaramillo & Jerome Novak
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

char time_slice_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/03/28 21:29:45  e_gourgoulhon
 * Evolution_std's renamed with suffix "_evol"
 * Method gam() modified
 * Added special constructor for derived classes.
 *
 * Revision 1.3  2004/03/26 13:33:02  j_novak
 * New methods for accessing/updating members (nn(), beta(), gam_uu(), k_uu(), ...)
 *
 * Revision 1.2  2004/03/26 08:22:56  e_gourgoulhon
 * Modifications to take into account the new setting of class
 * Evolution.
 *
 * Revision 1.1  2004/03/24 14:57:17  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "tensor.h"
#include "metric.h"
#include "evolution.h"



			    //--------------//
			    // Constructors //
			    //--------------//


// Standard constructor (Hamiltonian-like)
// ---------------------------------------
Time_slice::Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor kk_in,
               int depth_in) 
               : depth(depth_in),
		 scheme_order(depth_in-1),
                 jtime(0),
                 the_time(0., depth_in),
                 gam_dd_evol(depth_in),
                 gam_uu_evol(depth_in),
                 k_dd_evol(depth_in),
                 k_uu_evol(depth_in),
                 n_evol(lapse_in, depth_in),
                 beta_evol(shift_in, depth_in) {
                                  
    double time_init = the_time[jtime] ; 

    if (gamma_in.get_index_type(0) == COV) {
        gam_dd_evol.update(gamma_in, jtime, time_init) ; 
    }
    else {
        gam_uu_evol.update(gamma_in, jtime, time_init) ; 
    }
                 
    if (kk_in.get_index_type(0) == COV) {
        k_dd_evol.update(kk_in, jtime, time_init) ; 
    }
    else {
        k_uu_evol.update(kk_in, jtime, time_init) ; 
    }
                 
    set_der_0x0() ; 

}
                 
                 
// Standard constructor (Lagrangian-like)
// ---------------------------------------
Time_slice::Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Evolution_std<Sym_tensor>& gamma_in) 
               : depth(gamma_in.get_size()), 
 		 scheme_order(gamma_in.get_size()-1),
		 jtime(0),
                 the_time(0., gamma_in.get_size()),
                 gam_dd_evol( gamma_in.get_size() ),
                 gam_uu_evol( gamma_in.get_size() ),
                 k_dd_evol( gamma_in.get_size() ),
                 k_uu_evol( gamma_in.get_size() ),
                 n_evol(lapse_in, gamma_in.get_size() ),
                 beta_evol(shift_in, gamma_in.get_size() ) {

    cerr << 
    "Time_slice constuctor from evolution of gamma not implemented yet !\n" ;
    abort() ; 
                 
    set_der_0x0() ; 

}                 
                 
// Constructor as standard time slice of flat spacetime (Minkowski)
//-----------------------------------------------------------------               
Time_slice::Time_slice(const Map& mp, const Base_vect& triad, int depth_in)  
               : depth(depth_in),
		 scheme_order(depth_in-1),
                 jtime(0),
                 the_time(0., depth_in),
                 gam_dd_evol(depth_in),
                 gam_uu_evol(depth_in),
                 k_dd_evol(depth_in),
                 k_uu_evol(depth_in),
                 n_evol(depth_in),
                 beta_evol(depth_in) {
                 
    double time_init = the_time[jtime] ; 
    
    const Base_vect_spher* ptriad_s = 
        dynamic_cast<const Base_vect_spher*>(&triad) ;                  
    bool spher = (ptriad_s != 0x0) ; 
    
    if (spher) {
        gam_dd_evol.update( mp.flat_met_spher().cov(), jtime, time_init) ;  
    }         
    else {
        assert( dynamic_cast<const Base_vect_cart*>(&triad) != 0x0) ; 
        gam_dd_evol.update( mp.flat_met_cart().cov(), jtime, time_init) ;                           
    }


    // K_ij identically zero:
    Sym_tensor ktmp(mp, COV, triad) ;
    ktmp.set_etat_zero() ; 
    k_dd_evol.update(ktmp, jtime, time_init) ;  

    // Lapse identically one:
    Scalar tmp(mp) ; 
    tmp.set_etat_one() ; 
    n_evol.update(tmp, jtime, time_init) ; 
    
    // shift identically zero:
    Vector btmp(mp, CON, triad) ;
    btmp.set_etat_zero() ; 
    beta_evol.update(btmp, jtime, time_init) ;  
    
    set_der_0x0() ; 
}
                 

// Copy constructor
// ----------------
Time_slice::Time_slice(const Time_slice& tin) 
               : depth(tin.depth),
		 scheme_order(tin.scheme_order),
                 jtime(tin.jtime),
                 the_time(tin.the_time),
                 gam_dd_evol(tin.gam_dd_evol),
                 gam_uu_evol(tin.gam_uu_evol),
                 k_dd_evol(tin.k_dd_evol),
                 k_uu_evol(tin.k_uu_evol),
                 n_evol(tin.n_evol),
                 beta_evol(tin.beta_evol) {
                 
    set_der_0x0() ; 
}

// Special constructor for derived classes
//----------------------------------------               
Time_slice::Time_slice(int depth_in)  
               : depth(depth_in),
		 scheme_order(depth_in-1),
                 jtime(0),
                 the_time(0., depth_in),
                 gam_dd_evol(depth_in),
                 gam_uu_evol(depth_in),
                 k_dd_evol(depth_in),
                 k_uu_evol(depth_in),
                 n_evol(depth_in),
                 beta_evol(depth_in) {
                 
    set_der_0x0() ; 
}
                 



			    //--------------//
			    //  Destructor  //
			    //--------------//

Time_slice::~Time_slice(){

    Time_slice::del_deriv() ; 

}

                //---------------------//
                //  Memory management  //
                //---------------------//

void Time_slice::del_deriv() const {

    if (p_gamma != 0x0) delete p_gamma ; 
    
    set_der_0x0() ;
}


void Time_slice::set_der_0x0() const {

    p_gamma = 0x0 ; 
    
}


                //-----------------------//
                // Mutators / assignment //
                //-----------------------//

void Time_slice::operator=(const Time_slice& tin) {

    depth = tin.depth;
    scheme_order = tin.scheme_order ;
    jtime = tin.jtime;
    the_time = tin.the_time;
    gam_dd_evol = tin.gam_dd_evol;
    gam_uu_evol = tin.gam_uu_evol;
    k_dd_evol = tin.k_dd_evol;
    k_uu_evol = tin.k_uu_evol;
    n_evol = tin.n_evol;
    beta_evol = tin.beta_evol; 
    
    del_deriv() ; 
    
}


                //------------------//
                //      output      //
                //------------------//

ostream& operator<<(ostream& flux, const Time_slice& sigma) {

    flux << '\n' ;
    flux << "Lorene class : " << typeid(sigma).name() << '\n' ; 
    int jlast = sigma.jtime ; 
    flux << "Time label t = " << sigma.the_time[jlast] << '\n' ; 
    flux << "Index of time step j = " << jlast << '\n' ; 
    flux << "------------------------------------------------------------\n" 
        <<  "Max. of absolute values of the various fields in each domain: \n" ;
    if (sigma.gam_dd_evol.is_known(jlast)) {
        maxabs( sigma.gam_dd_evol[jlast], "gam_dd", flux) ;
    }
    if (sigma.gam_uu_evol.is_known(jlast)) {
        maxabs( sigma.gam_uu_evol[jlast], "gam_uu", flux) ;
    }
    if (sigma.k_dd_evol.is_known(jlast)) {
        maxabs( sigma.k_dd_evol[jlast], "k_dd", flux) ;
    }
    if (sigma.k_uu_evol.is_known(jlast)) {
        maxabs( sigma.k_uu_evol[jlast], "k_uu", flux) ;
    }
    if (sigma.n_evol.is_known(jlast)) {
        maxabs( sigma.n_evol[jlast], "N", flux) ;
    }
    if (sigma.beta_evol.is_known(jlast)) {
        maxabs( sigma.beta_evol[jlast], "beta", flux) ;
    }

    if (sigma.p_gamma != 0x0) flux << *sigma.p_gamma << endl ; 
    
    return flux ; 

}




















                
                
                
                

