/*
 *  Methods of class Tslice_dirac_max
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

char tslice_dirac_max_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/03/30 14:00:31  j_novak
 * New class Tslide_dirac_max (first version).
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


// Constructor from conformal decomposition
// ----------------------------------------

Tslice_dirac_max::Tslice_dirac_max(const Scalar& lapse_in, const Vector& shift_in,
            const Metric_flat& ff_in, const Scalar& psi_in, 
            const Sym_tensor_trans& hh_in, const Sym_tensor aa_in, 
            int depth_in) 
  : Time_slice_conf( lapse_in, shift_in, ff_in, psi_in, hh_in, aa_in, 
		     0*lapse_in, depth_in), 
    khi_evol(hh_in.tt_part().khi(), depth_in), 
    mu_evol(hh_in.tt_part().mu(), depth_in),
    trh_evol(hh_in.the_trace(), depth_in) {

}
                 



// Copy constructor
// ----------------

Tslice_dirac_max::Tslice_dirac_max(const Tslice_dirac_max& tin) 
                    : Time_slice_conf(tin), 
                      khi_evol(tin.khi_evol), 
                      mu_evol(tin.mu_evol),
                      trh_evol(tin.trh_evol) {

}
			    //--------------//
			    //  Destructor  //
			    //--------------//

Tslice_dirac_max::~Tslice_dirac_max(){

}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Tslice_dirac_max::operator=(const Tslice_dirac_max& tin) {

    Time_slice_conf::operator=(tin) ; 

    khi_evol = tin.khi_evol ; 
    mu_evol = tin.mu_evol ; 
    trh_evol = tin.trh_evol ; 
       
}


                //----------------------------------------------------//
                //  Update of fields from base class Time_slice_conf  //
                //----------------------------------------------------//


const Sym_tensor& Tslice_dirac_max::hh() const {

  if (!( hh_evol.is_known(jtime) ) ) {
    assert (khi_evol.is_known(jtime)) ;
    assert (mu_evol.is_known(jtime)) ;
    assert (trh_evol.is_known(jtime)) ;

    const Map& mapping = khi_evol[jtime].get_mp() ;

    Sym_tensor_tt hh_tt_tmp(mapping, *ff.get_triad(), ff) ;

    hh_tt_tmp.set_khi_mu( khi(), mu() ) ;
    
    Sym_tensor_trans hh_t_tmp(mapping, *ff.get_triad(), ff) ;

    hh_t_tmp.set_tt_trace( hh_tt_tmp, trh() ) ;

    Vector vec_tmp(mapping, CON,  *ff.get_triad()) ; //The longitudinal part 
    vec_tmp.set_etat_zero() ; // is zero (Dirac gauge).

    Sym_tensor hh_tmp(mapping, CON, *ff.get_triad()) ;
    
    hh_tmp.set_longit_trans(vec_tmp, hh_t_tmp) ;
    
    hh_evol.update(hh_tmp, jtime, the_time[jtime] ) ;
    
  }
  
  return hh_evol[jtime] ; 

}


const Scalar& Tslice_dirac_max::trk() const {

    if( !(trk_evol.is_known(jtime)) ) {

      Scalar resu(ff.get_mp()) ;
      resu.set_etat_zero() ;
      
      trk_evol.update(resu, jtime, the_time[jtime]) ;

    } 
    
    return trk_evol[jtime] ; 

}


const Vector& Tslice_dirac_max::hdirac() const {

    if (p_hdirac == 0x0) {
        p_hdirac = new Vector(ff.get_mp(), CON, ff.get_triad() ) ;
	p_hdirac->set_etat_zero() ;
    }
    
    return *p_hdirac ; 

}




                //-----------------------------------//
                //  Update of fields from this class //
                //-----------------------------------//


const Scalar& Tslice_dirac_max::khi() const {

    if (!( khi_evol.is_known(jtime) ) ) {

        assert( hh_evol.is_known(jtime) ) ; 

	
        
        khi_evol.update( hh().transverse(ff).tt_part().khi(), 
			 jtime, the_time[jtime] ) ; 
    }

    return khi_evol[jtime] ;

} 

const Scalar& Tslice_dirac_max::mu() const {

    if (!( mu_evol.is_known(jtime) ) ) {

      assert( hh_evol.is_known(jtime) ) ; 
      
        
        mu_evol.update( hh().transverse(ff).tt_part().mu(), 
			jtime, the_time[jtime] ) ; 
    }

    return mu_evol[jtime] ;


}


const Scalar& Tslice_dirac_max::trh() const {

  int it_max = 100 ;
  int it ;

  double precis = 1.e-12 ;

  if( !(trh_evol.is_known(jtime)) ) {
    
    assert ( khi_evol.is_known(jtime) ) ;
    assert ( mu_evol.is_known(jtime) );
    
    const Map& mapping = khi_evol[jtime].get_mp() ;
    
    Sym_tensor_tt hijtt(mapping, *ff.get_triad(), ff) ;
    hijtt.set_khi_mu( khi_evol[jtime], mu_evol[jtime] ) ;
    
    Sym_tensor_trans hij = hijtt ; //Start of iteration
    Scalar htrace_jm1(mapping) ;
    htrace_jm1.set_etat_zero() ;
    
    for (it=0; it<it_max; it++) {
      
      Scalar htrace = hij(1,1) * hij(2,3) * hij(2,3) 
	+ hij(2,2) * hij(1,3) * hij(1,3) + hij(3,3) * hij(1,2) * hij(1,2)
	- 2 * hij(1,2) * hij(1,3) * hij(2,3) - hij(1,1) * hij(2,2) * hij(3,3)
	+ hij(1,2) * hij(1,2) + hij(1,3) * hij(1,3) + hij(2,3) * hij(2,3)
	- hij(1,1) * hij(2,2) - hij(1,1) * hij(3,3) - hij(2,2) * hij(3,3) ;

      hij.set_tt_trace( hijtt, htrace) ; //Solves the Poisson equation for Phi

      double diff = max(max(abs(htrace - htrace_jm1))) ;
      cout << "Tslide_dirac_max::trh() : " 
	   << "Iteration : " << it << " difference : " << diff << endl ;
      if (diff < precis) break ;
      else htrace_jm1 = htrace ;
    }
    
    if (it == it_max) {
      cout << " No convergence reached in Tslide_dirac_max::trh() " << '\n' ;
      cout << "Required accuracy : " << precis << endl ;
      abort() ;
    }

    trh_evol.update(hij.the_trace(), jtime, the_time[jtime] ) ;

  }
    
  return trh_evol[jtime] ; 

}



                //------------------//
                //      output      //
                //------------------//

ostream& Tslice_dirac_max::operator>>(ostream& flux) const {

    Time_slice_conf::operator>>(flux) ; 

    flux << "Dirac gauge and maximal slicing" << '\n' ;

    if (khi_evol.is_known(jtime)) {
        maxabs( khi_evol[jtime], "Khi", flux) ;
    }
    if (mu_evol.is_known(jtime)) {
        maxabs( mu_evol[jtime], "Mu", flux) ;
    }
    if (trh_evol.is_known(jtime)) {
        maxabs( trh_evol[jtime], "tr h", flux) ;
    }
    
    return flux ; 

}



                
                
                
                

