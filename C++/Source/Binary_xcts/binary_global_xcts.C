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

    Vector dpsi(et[0]->get_dcov_Psi()) ;
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

  double massadm ;
  massadm = 1. ;

  //## : to be done 
  
  for (int i=0; i<=1; i++) {	    // loop on the stars

    // Declaration of all fields
      const Scalar& psi(et[i]->get_Psi()) ;

      Scalar psi4 = pow(psi, 4.) ; 
      psi4.std_spectral_base() ;
	  
	  Scalar spsi20 = pow(psi, -20.) ; 
      spsi20.std_spectral_base() ;

      const Scalar& ener_euler = et[i]->get_ener_euler() ;
      const Scalar& hacar_auto = et[i]->get_hacar_auto() ;
      const Scalar& hacar_comp = et[i]->get_hacar_comp() ;
 
      // Source in IWM approximation 
	  Scalar source =  psi4 % psi % ener_euler ; 
	                + spsi20 % (hacar_auto + hacar_comp)/(4.*qpig) ;  
	  
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
   
    Map_af map0 (et[0]->get_mp()) ; 
    
    Vector vect = (et[0]->get_dcov_chi())/(et[0]->get_chi()) 
                - (et[0]->get_dcov_Psi())/(et[0]->get_Psi()) ;
    
    vect.change_triad(map0.get_bvect_spher()) ;
    Scalar integrant (vect(1)) ;
    
    *p_mass_kom = map0.integrale_surface_infini (integrant) / qpig ;
   
   
    
  }	// End of the case where a new computation was necessary
        
  return *p_mass_kom ; 
    
}

double Binary_xcts::mass_kom_vol() const {
    
  using namespace Unites ;

  double masskom ;
  masskom = 1. ;

  //## : to be done 
  /*
  for (int i=0; i<=1; i++) {	    // loop on the stars

     // Declaration of all fields
      const Scalar& psi4 = et[i]->get_psi4() ;
      const Scalar& ener_euler = et[i]->get_ener_euler() ;
      const Scalar& s_euler = et[i]->get_s_euler() ;
      const Scalar& kcar_auto = et[i]->get_kcar_auto() ;
      const Scalar& kcar_comp = et[i]->get_kcar_comp() ;
      const Metric& flat = et[i]->get_flat() ;
      const Sym_tensor& hij = et[i]->get_hij() ;
      const Scalar& logn = et[i]->get_logn_auto() + et[i]->get_logn_comp() ;
      const Scalar& logn_auto = et[i]->get_logn_auto() ;
      Scalar nn = exp(logn) ;
      nn.std_spectral_base() ;
      
      const Tensor& dcov_logn_auto = logn_auto.derive_cov(flat) ;
      const Vector& dcov_logn = et[i]->get_dcov_logn() ;
      const Vector& dcon_logn = et[i]->get_dcon_logn() ;
      const Vector& dcov_phi = et[i]->get_dcov_phi() ;
      Tensor dcovdcov_logn_auto = (logn_auto.derive_cov(flat))
	.derive_cov(flat) ;
      dcovdcov_logn_auto.inc_dzpuis() ;

      Scalar source = qpig * psi4 % (ener_euler + s_euler) ;
      source += psi4 % (kcar_auto + kcar_comp) ;

      source = source / qpig * nn  ;
  
      source.std_spectral_base() ;

      masskom += source.integrale() ;
	  
  }

  */
  return masskom ;

}


			//---------------------------------//
		    //	 Total angular momentum        //
		    //---------------------------------//

const Tbl& Binary_xcts::angu_mom() const {

  using namespace Unites ;
	
	//## : to be done 	
	
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
		double sigma = 1.*r0 ;
		
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
		const Scalar& psi = et[i]->get_Psi() ; 
		Scalar rho = pow(psi, 0.625) * (ee + pp) ; 
		rho.std_spectral_base() ;

		Vector jmom = rho * (et[i]->get_u_euler()) ; 
				
		const Metric_flat flat (mp.flat_met_cart()) ; 
		
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
