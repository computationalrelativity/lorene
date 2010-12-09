/*
 * Methods of class Binary_xcts to compute global quantities
 * (see file binary_xcts.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

char binary_global_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2010/12/09 10:39:41  m_bejger
 * Further corrections to integral quantities
 *
 * Revision 1.6  2010/10/26 19:16:26  m_bejger
 * Cleanup of some diagnostic messages
 *
 * Revision 1.5  2010/10/25 15:02:08  m_bejger
 * mass_kom_vol() corrected
 *
 * Revision 1.4  2010/10/24 21:45:24  m_bejger
 * mass_adm() corrected
 *
 * Revision 1.3  2010/06/17 14:48:14  m_bejger
 * Minor corrections
 *
 * Revision 1.2  2010/06/04 19:54:19  m_bejger
 * Minor corrections, mass volume integrals need to be checked out
 *
 * Revision 1.1  2010/05/04 07:35:54  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "nbr_spx.h"
#include "binary_xcts.h"
#include "unites.h"
#include "metric.h"

		    //---------------------------------//
		    //		ADM mass                   //
		    //---------------------------------//

double Binary_xcts::mass_adm() const {
    
  using namespace Unites ;
  
  if (p_mass_adm == 0x0) {	    // a new computation is required
    
    p_mass_adm = new double ; 
	*p_mass_adm = 0 ; 
	
    const Map_af map0 (et[0]->get_mp()) ;
    const Metric& flat = (et[0]->get_flat()) ;

    Vector dpsi((et[0]->get_Psi()).derive_cov(flat)) ;

    dpsi.change_triad(map0.get_bvect_spher()) ;

    Scalar integrand ( dpsi(1) ) ;

    *p_mass_adm = - 2.* map0.integrale_surface_infini(integrand)/ qpig ;
    
	}

    return *p_mass_adm ; 
    
}

			//---------------------------------//
			//          ADM mass               //
			//          (volume integral)      // 
			//---------------------------------//

double Binary_xcts::mass_adm_vol() const {

  using namespace Unites ;

  double massadm = 0. ;

  for (int i=0; i<=1; i++) {	    // loop on the stars

    // Declaration of all fields

      const Scalar& psi(et[i]->get_Psi()) ;
      Scalar psi5 = pow(psi, 5.) ;
      psi5.std_spectral_base() ;
	  
	  Scalar spsi7 = pow(psi, -7.) ; 
      spsi7.std_spectral_base() ;

      const Scalar& ener_euler = et[i]->get_ener_euler() ;
      const Scalar& hacar_auto = et[i]->get_hacar_auto() ;
      const Scalar& hacar_comp = et[i]->get_hacar_comp() ;
 
	  Scalar source = psi5 % ener_euler 
	  			    + spsi7 % (hacar_auto + hacar_comp)/(4.*qpig) ;  
	  
      source.std_spectral_base() ;

      massadm += source.integrale() ;
  }
  
  return massadm ;
}

		    //---------------------------------//
		    //		Komar mass                 //
		    //---------------------------------//

double Binary_xcts::mass_kom() const {
    
  using namespace Unites ;

  if (p_mass_kom == 0x0) {	    // a new computation is requireed

    p_mass_kom = new double ;    
    *p_mass_kom = 0 ; 

    const Scalar& logn = et[0]->get_logn() ;
    const Metric& flat = et[0]->get_flat() ; 

	Map_af map0 (et[0]->get_mp()) ; 
	  
    Vector vect = logn.derive_con(flat) ; 
    vect.change_triad(map0.get_bvect_spher()) ;

    Scalar integrant (vect(1)) ;
    
    *p_mass_kom = map0.integrale_surface_infini (integrant) / qpig ;
       
  }	// End of the case where a new computation was necessary
        
  return *p_mass_kom ; 
    
}

double Binary_xcts::mass_kom_vol() const {
    
  using namespace Unites ;

  double masskom = 0.;

  for (int i=0; i<=1; i++) {	    // loop on the stars

      const Scalar& Psi  = et[i]->get_Psi() ;	  
      const Scalar& psi4 = et[i]->get_psi4() ;
      const Scalar& ener_euler = et[i]->get_ener_euler() ;
      const Scalar& s_euler = et[i]->get_s_euler() ;
	  const Scalar& chi  = et[i]->get_chi() ; 
	  const Scalar& hacar_auto = et[i]->get_hacar_auto() ;
      const Scalar& hacar_comp = et[i]->get_hacar_comp() ;  

	  Scalar psi4chi = psi4 % chi ; 
	  psi4chi.std_spectral_base() ; 

	  Scalar source = 0.5* ener_euler * (psi4chi + pow(Psi, 5.) ) 
		  + psi4chi * s_euler + pow(Psi, -7.) * (7.*chi/Psi + 1.) 
		  			* (hacar_auto + hacar_comp)/ (8.*qpig) ; 
	  source.std_spectral_base() ;
	  

      masskom += source.integrale() ;
	  
  }

  return masskom ;

}

			//---------------------------------//
		    //	 Total angular momentum        //
		    //---------------------------------//
	
//## to be checked 
const Tbl& Binary_xcts::angu_mom() const {

  using namespace Unites ;
	
	if (p_angu_mom == 0x0) {	    // a new computation is requireed
    
	p_angu_mom = new Tbl(3) ; 
	p_angu_mom->annule_hard() ;	// fills the double array with zeros

	// Reference Cartesian vector basis of the Absolute frame
	Base_vect_cart bvect_ref(0.) ; 	// 0. = parallel to the Absolute frame
	
	for (int i=0; i<=1; i++) {	    // loop on the stars

		const Map& mp = et[i]->get_mp() ; 
		int nzm1 = mp.get_mg()->get_nzone() - 1 ; 
		
		// Function exp(-(r-r_0)^2/sigma^2)
		// --------------------------------
		
		double r0 = mp.val_r(nzm1-1, 1, 0, 0) ;
		double sigma = r0 ;
		
		Scalar rr (mp) ;
		rr = mp.r ;
		
		Scalar ff (mp) ;
		ff = exp( -(rr - r0)*(rr - r0)/sigma/sigma ) ;
		for (int ii=0; ii<nzm1; ii++)
		  ff.set_domain(ii) = 1. ;
		ff.set_outer_boundary(nzm1, 0) ;
		ff.std_spectral_base() ;

		// Azimuthal vector d/dphi 
		Vector vphi(mp, CON, bvect_ref) ; 		
		Scalar yya (mp) ;
		yya = mp.ya ;
		Scalar xxa (mp) ;
		xxa = mp.xa ;
		vphi.set(1) = - yya * ff ; 	// phi^X
		vphi.set(2) = xxa * ff ; 
		vphi.set(3) = 0 ;  

		vphi.set(1).set_outer_boundary(nzm1, 0) ;
		vphi.set(2).set_outer_boundary(nzm1, 0) ;
	
		vphi.std_spectral_base() ; 
		vphi.change_triad(mp.get_bvect_cart()) ; 
		
		// Matter part
		// -----------
		const Scalar& ee = et[i]->get_ener_euler() ;  
		const Scalar& pp = et[i]->get_press() ;
		Scalar rho = pow(et[i]->get_psi4(), 2.5) * (ee + pp) ; 
		rho.std_spectral_base() ;

		Vector jmom = rho * (et[i]->get_u_euler()) ; 
	    const Metric& flat = et[i]->get_flat() ;

		Vector vphi_cov = vphi.up_down(flat) ;
		
		Scalar integrand = contract(jmom, 0, vphi_cov, 0) ; 
		      
		p_angu_mom->set(2) += integrand.integrale() ;
		
	}  // End of the loop on the stars

    }	// End of the case where a new computation was necessary
  
  	return *p_angu_mom ; 
  
}

		    //---------------------------------//
		    //		Total energy	           //
		    //---------------------------------//

double Binary_xcts::total_ener() const {
   
    if (p_total_ener == 0x0) {	    // a new computation is requireed
	
	p_total_ener = new double ; 
	    
	    *p_total_ener = mass_adm() - star1.mass_b() - star2.mass_b() ; 
	    
    }	// End of the case where a new computation was necessary
    
    return *p_total_ener ; 
    
}


		    //---------------------------------//
		    //	 Error on the virial theorem   //
		    //---------------------------------//

double Binary_xcts::virial() const {
    
    if (p_virial == 0x0) {	    // a new computation is requireed
	
	p_virial = new double ; 
	    
	    *p_virial = 1. - mass_kom() / mass_adm() ; 
	    
	}
    
    return *p_virial ; 
    
}

		    //---------------------------------//
		    //	 Error on the virial theorem   // 
			//   (volume version)              //
		    //---------------------------------//

double Binary_xcts::virial_vol() const {
    
    if (p_virial_vol == 0x0) {	    // a new computation is requireed
	
	p_virial_vol = new double ; 
	    
	    *p_virial_vol = 1. - mass_kom_vol() / mass_adm_vol() ; 
	    
	}
    
    return *p_virial_vol ; 
    
}
