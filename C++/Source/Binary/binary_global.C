/*
 * Methods of class Binary to compute global quantities
 *
 * (see file binary.h for documentation)
 */

/*
 *   Copyright (c) 2004 Francois Limousin
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char binary_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.12  2005/09/15 14:41:04  e_gourgoulhon
 * The total angular momentum is now computed via a volume integral.
 *
 * Revision 1.11  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.10  2005/04/08 12:36:45  f_limousin
 * Just to avoid warnings...
 *
 * Revision 1.9  2005/02/17 17:35:00  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.8  2004/07/21 11:46:24  f_limousin
 * Add function mass_adm_vol() to compute the ADM mass of the system
 * with a volume integral instead of a surface one.
 *
 * Revision 1.7  2004/05/25 14:25:53  f_limousin
 * Add the virial theorem for conformally flat configurations.
 *
 * Revision 1.6  2004/03/31 12:44:54  f_limousin
 * Minor modifs.
 *
 * Revision 1.5  2004/03/25 10:29:01  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2004/02/27 10:25:30  f_limousin
 * Modif. to avoid an error in compilation.
 *
 * Revision 1.3  2004/02/27 10:03:04  f_limousin
 * The computation of mass_adm() and mass_komar() is now OK !
 *
 * Revision 1.2  2004/01/20 15:21:36  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */


// Headers C
#include "math.h"

// Headers Lorene
#include "nbr_spx.h"
#include "binary.h"
#include "unites.h"

		    //---------------------------------//
		    //		ADM mass	       //
		    //---------------------------------//

double Binary::mass_adm() const {
    
  using namespace Unites ;
  if (p_mass_adm == 0x0) {	    // a new computation is requireed
    
    p_mass_adm = new double ; 
	    
    *p_mass_adm = 0 ; 
    
    // First Method
    /*
    const Metric& gamij = (et[0]->get_gamma()) ;
    const Metric& flat = (et[0]->get_flat()) ;
    const Map_af map0 (et[0]->get_mp()) ;
    
    Sym_tensor gamij_cov = gamij.cov() ;
    const Tensor& dcov_gamij = gamij_cov.derive_cov(flat) ;
    
    Vector dgamma_1 =  contract( flat.con(), 1, contract(contract(flat.con()
       			   , 0, dcov_gamij, 2), 0, 2), 0) ;
    dgamma_1.change_triad(map0.get_bvect_spher()) ;
    Scalar integrant_1 (dgamma_1(1)) ;
     
    Vector dgamma_2 =  contract( flat.con(), 1, contract(contract(flat.con()
			   , 0, dcov_gamij, 0), 0, 1), 0) ;
    dgamma_2.change_triad(map0.get_bvect_spher()) ;
    Scalar integrant_2 (dgamma_2(1)) ;

    *p_mass_adm = map0.integrale_surface_infini (integrant_1 - integrant_2)/(4*qpig) ;
    */
    // Second Method
    
    const Map_af map0 (et[0]->get_mp()) ;
    const Metric& flat = (et[0]->get_flat()) ;

    Vector dpsi(0.5*(et[0]->get_lnq() - 
    		     et[0]->get_logn()).derive_cov(flat)) ;
  
    Vector ww (0.125*(contract(et[0]->get_hij().derive_cov(flat), 1, 2) 
    		      - (et[0]->get_hij().trace(flat)).derive_con(flat))) ;
    
    dpsi.change_triad(map0.get_bvect_spher()) ;
    ww.change_triad(map0.get_bvect_spher()) ;

    // ww = 0 in Dirac gauge (Eq 174 of BGGN)
    Scalar integrand (dpsi(1) + 0*ww(1)) ;

    *p_mass_adm = map0.integrale_surface_infini (integrand) / (-qpig/2.) ;
    
    }	// End of the case where a new computation was necessary
    
    return *p_mass_adm ; 
    
}

double Binary::mass_adm_vol() const {

  using namespace Unites ;

  double massadm ;
  massadm = 0. ;

  for (int i=0; i<=1; i++) {	    // loop on the stars

      const Scalar& psi4 = et[i]->get_psi4() ;
      Scalar psi (pow(psi4, 0.25)) ;
      psi.std_spectral_base() ;
      const Scalar& ener_euler = et[i]->get_ener_euler() ;
      const Scalar& kcar_auto = et[i]->get_kcar_auto() ;
      const Scalar& kcar_comp = et[i]->get_kcar_comp() ;
      const Metric& gtilde = et[i]->get_gtilde() ;
      const Metric& flat = et[i]->get_flat() ;
      const Sym_tensor& hij = et[i]->get_hij() ;
      const Sym_tensor& hij_auto = et[i]->get_hij_auto() ;
      const Vector& dcov_logn = et[i]->get_dcov_logn() ;
      const Vector& dcov_phi = et[i]->get_dcov_phi() ;
      const Vector& dcov_lnq = 2*dcov_phi + dcov_logn ;
      const Scalar& lnq_auto = et[i]->get_lnq_auto() ;
      const Scalar& logn_auto = et[i]->get_logn_auto() ;
      const Scalar& phi_auto = 0.5 * (lnq_auto - logn_auto) ;

      const Tensor& dcov_hij_auto = hij_auto.derive_cov(flat) ;
      const Tensor& dcov_gtilde = gtilde.cov().derive_cov(flat) ;
      const Tensor& dcov_phi_auto = phi_auto.derive_cov(flat) ;
      const Tensor& dcov_logn_auto = logn_auto.derive_cov(flat) ;
      const Tensor& dcov_lnq_auto = lnq_auto.derive_cov(flat) ;
      Tensor dcovdcov_lnq_auto = lnq_auto.derive_cov(flat).derive_cov(flat) ;
      dcovdcov_lnq_auto.inc_dzpuis() ;
      Tensor dcovdcov_logn_auto = logn_auto.derive_cov(flat).derive_cov(flat) ;
      dcovdcov_logn_auto.inc_dzpuis() ;
      
      Scalar source = - psi4 % (qpig*ener_euler + (kcar_auto + kcar_comp)/4.) 
	- 0*2*contract(contract(gtilde.con(), 0, dcov_phi, 0), 0, dcov_phi_auto, 0, true) ;
      
      source += 4*contract(hij, 0, 1, dcov_logn * dcov_phi_auto, 0, 1) +
	2*contract(hij, 0, 1, dcov_phi * dcov_phi_auto, 0, 1) +
	0.0625 * contract(gtilde.con(), 0, 1, contract(
			   dcov_hij_auto, 0, 1, dcov_gtilde, 0, 1), 0, 1) - 
           0.125 * contract(gtilde.con(), 0, 1, contract(dcov_hij_auto, 
			      0, 1, dcov_gtilde, 0, 2), 0, 1) -
	 contract(hij,0,1,dcovdcov_lnq_auto + dcov_lnq_auto*dcov_lnq,0,1) +
	 contract(hij,0,1,dcovdcov_logn_auto + dcov_logn_auto*dcov_logn,0,1) ;

      source = source * psi ;

      source.std_spectral_base() ;

      massadm += - source.integrale()/qpig ;

/*
      Scalar temp1 (contract(phi_auto.derive_con(flat), 0, dcov_phi, 0)) ;
      Scalar temp3 (contract(dcon_phi, 0, dcov_phi, 0)) ;
      Scalar temp4 (contract(dcon_logn, 0, dcov_phi_auto, 0)) ;
      Scalar temp5 (contract(logn_auto.derive_con(flat), 0, dcov_phi, 0)) ;
      Scalar temp6 (contract(hij, 0, 1, dcov_logn * dcov_phi_auto, 0, 1)) ;
      Scalar temp7 (contract(hij, 0, 1, dcov_logn_auto * dcov_phi, 0, 1)) ;
      Scalar temp8 (contract(hij_auto, 0, 1, dcov_logn * dcov_phi, 0, 1)) ;
	  
      cout << "temp1 = " << -temp1.integrale() << endl ;
      cout << "temp3 = " << -temp3.integrale() << endl ;
      cout << "temp4 = " << -temp4.integrale() << endl ;
      cout << "temp5 = " << -temp5.integrale() << endl ;
      cout << "temp6 = " << -temp6.integrale() << endl ;
      cout << "temp7 = " << -temp7.integrale() << endl ;
      cout << "temp8 = " << -temp8.integrale() << endl ;
*/
   
  }

  return massadm ;

}

		    //---------------------------------//
		    //		Komar mass	       //
		    //---------------------------------//

double Binary::mass_kom() const {
    
  using namespace Unites ;

  if (p_mass_kom == 0x0) {	    // a new computation is requireed
    
    p_mass_kom = new double ; 
      
    *p_mass_kom = 0 ; 
    
    const Tensor& logn = et[0]->get_logn() ;
    const Metric& flat = (et[0]->get_flat()) ;
    const Sym_tensor&  hij = (et[0]->get_hij()) ;
    Map_af map0 (et[0]->get_mp()) ; 
    
    Vector vect = logn.derive_con(flat) + 
                       contract(hij, 1, logn.derive_cov(flat), 0) ;
    vect.change_triad(map0.get_bvect_spher()) ;
    Scalar integrant (vect(1)) ;
    
    *p_mass_kom = map0.integrale_surface_infini (integrant) / qpig ;
    
  }	// End of the case where a new computation was necessary
    
  return *p_mass_kom ; 
    
}

double Binary::mass_kom_vol() const {
    
  using namespace Unites ;

  double masskom ;
  masskom = 0. ;

  for (int i=0; i<=1; i++) {	    // loop on the stars
      const Scalar& psi4 = et[i]->get_psi4() ;
      const Scalar& ener_euler = et[i]->get_ener_euler() ;
      const Scalar& s_euler = et[i]->get_s_euler() ;
      const Scalar& kcar_auto = et[i]->get_kcar_auto() ;
      const Scalar& kcar_comp = et[i]->get_kcar_comp() ;
      const Metric& gtilde = et[i]->get_gtilde() ;
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
      source += - 0*contract(dcov_logn_auto, 0, dcon_logn, 0, true) 
	  - 2. * contract(contract(gtilde.con(), 0, dcov_phi, 0), 0, dcov_logn_auto, 0, true) ;
      source += - contract(hij, 0, 1, dcovdcov_logn_auto + 
			   dcov_logn_auto*dcov_logn, 0, 1) ;

      source = source / qpig * nn  ;
  
      source.std_spectral_base() ;

      masskom += source.integrale() ;
	  
  }

  return masskom ;

}


		    //---------------------------------//
		    //	 Total angular momentum        //
		    //---------------------------------//

const Tbl& Binary::angu_mom() const {

  using namespace Unites ;

  /*
    if (p_angu_mom == 0x0) {	    // a new computation is requireed
	
      p_angu_mom = new Tbl(3) ; 
      
      p_angu_mom->annule_hard() ;	// fills the double array with zeros
  
      const Sym_tensor& kij_auto = et[0]->get_tkij_auto() ;
      const Sym_tensor& kij_comp = et[0]->get_tkij_comp() ;
      const Tensor& psi4 = et[0]->get_psi4() ;
      const Map_af map0 (kij_auto.get_mp()) ;

      Sym_tensor kij = (kij_auto + kij_comp) / psi4 ;
      kij.change_triad(map0.get_bvect_cart()) ;
 
      // X component
      // -----------

      Vector vect_x(et[0]->get_mp(), CON, map0.get_bvect_cart()) ;      
       
      for (int i=1; i<=3; i++) {

	  Scalar kij_1 = kij(3, i) ;
	  Scalar kij_2 = kij(2, i) ;
	  	  
	  kij_1.mult_rsint() ;
	  Valeur vtmp = kij_1.get_spectral_va().mult_sp() ;
	  kij_1.set_spectral_va() = vtmp ; 
	  
	  kij_2.mult_r() ;
	  vtmp = kij_2.get_spectral_va().mult_ct() ;
	  kij_2.set_spectral_va() = vtmp ; 

	  vect_x.set(i) = kij_1 - kij_2 ;
      }
 
      vect_x.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_x (vect_x(1)) ;
      
      p_angu_mom->set(0) = map0.integrale_surface_infini (integrant_x) 
	                  / (8*M_PI) ;
      
      // Y component
      // -----------
      
      Vector vect_y(et[0]->get_mp(), CON, map0.get_bvect_cart()) ;      
   
      for (int i=1; i<=3; i++) {

	  Scalar kij_1 = kij(1, i) ;
	  Scalar kij_2 = kij(3, i) ;	  
	  
	  kij_1.mult_r() ;
	  Valeur vtmp = kij_1.get_spectral_va().mult_ct() ;
	  kij_1.set_spectral_va() = vtmp ; 
	  
	  kij_2.mult_rsint() ;
	  vtmp = kij_2.get_spectral_va().mult_cp() ;
	  kij_2.set_spectral_va() = vtmp ; 
	 
	  vect_y.set(i) = kij_1 - kij_2 ;
      }

      vect_y.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_y (vect_y(1)) ;
      
      p_angu_mom->set(1) = map0.integrale_surface_infini (integrant_y) 
	                  / (8*M_PI) ;
      
      // Z component
      // -----------

      Vector vect_z(et[0]->get_mp(), CON, map0.get_bvect_cart()) ;      

      for (int i=1; i<=3; i++) {

	  Scalar kij_1 = kij(2, i) ;
	  Scalar kij_2 = kij(1, i) ;	  	  

	  kij_1.mult_rsint() ;
	  Valeur vtmp = kij_1.get_spectral_va().mult_cp() ;
	  kij_1.set_spectral_va() = vtmp ; 

	  kij_2.mult_rsint() ;
	  vtmp =  kij_2.get_spectral_va().mult_sp() ;
	  kij_2.set_spectral_va() = vtmp ;

	  vect_z.set(i) = kij_1 - kij_2 ;
      }
       
      vect_z.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_z (vect_z(1)) ;
      
      p_angu_mom->set(2) = map0.integrale_surface_infini (integrant_z) 
	                 ;// (8*M_PI) ;
      
      
    }	// End of the case where a new computation was necessary
  */
  

	/*  
  if (p_angu_mom == 0x0) {	    // a new computation is requireed
    p_angu_mom = new Tbl(3) ; 
    p_angu_mom->annule_hard() ;	// fills the double array with zeros
    p_angu_mom->set(0) = 0. ;
    p_angu_mom->set(1) = 0. ;

    // Alignement 
    double orientation_un = et[0]->get_mp().get_rot_phi() ;
    assert ((orientation_un==0) || (orientation_un==M_PI)) ;
    double orientation_deux = et[1]->get_mp().get_rot_phi() ;
    assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
    int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;
    
    // Construction of an auxiliar grid and mapping
    int nzones = et[0]->get_mp().get_mg()->get_nzone() ;
    double* bornes = new double [nzones+1] ;
    double courant = (et[0]->get_mp().get_ori_x()-et[0]->get_mp().get_ori_x())+1 ;
    for (int i=nzones-1 ; i>0 ; i--) {
      bornes[i] = courant ;
      courant /= 2. ;
    }
    bornes[0] = 0 ;
    bornes[nzones] = __infinity ;
    
    Map_af mapping (*(et[0]->get_mp().get_mg()), bornes) ;
    
    delete [] bornes ; 
    
    // Construction of k_total
    Sym_tensor k_total (mapping, CON, mapping.get_bvect_cart()) ;
    
    Vector shift_un (mapping, CON, mapping.get_bvect_cart()) ;
    Vector shift_deux (mapping, CON, mapping.get_bvect_cart()) ;
    
    Vector beta_un (et[0]->get_beta_auto()) ;
    Vector beta_deux (et[1]->get_beta_auto()) ;
    beta_un.change_triad(et[0]->get_mp().get_bvect_cart()) ;
    beta_deux.change_triad(et[1]->get_mp().get_bvect_cart()) ;
    beta_un.std_spectral_base() ;
    beta_deux.std_spectral_base() ;
    
    shift_un.set(1).import(beta_un(1)) ;
    shift_un.set(2).import(beta_un(2)) ;
    shift_un.set(3).import(beta_un(3)) ;
 
    shift_deux.set(1).import(same_orient*beta_deux(1)) ;
    shift_deux.set(2).import(same_orient*beta_deux(2)) ;
    shift_deux.set(3).import(beta_deux(3)) ;
    
    Vector shift_tot (shift_un+shift_deux) ;
    shift_tot.std_spectral_base() ;
    shift_tot.annule(0, nzones-2) ;
    
    
    // Substract the residuals
    shift_tot.inc_dzpuis(2) ;
    shift_tot.dec_dzpuis(2) ;
    
    
    Sym_tensor temp_gamt (et[0]->get_gtilde().cov()) ;
    temp_gamt.change_triad(mapping.get_bvect_cart()) ;
    Metric gamt_cart (temp_gamt) ;
    
    k_total = shift_tot.ope_killing_conf(gamt_cart) / 2. ;
    
    for (int lig=1 ; lig<=3 ; lig++)
    for (int col=lig ; col<=3 ; col++)
      k_total.set(lig, col).mult_r_ced() ;
    
    Vector vecteur_un (mapping, CON, mapping.get_bvect_cart()) ;
    for (int i=1 ; i<=3 ; i++)
      vecteur_un.set(i) = k_total(1, i) ;
    vecteur_un.change_triad (mapping.get_bvect_spher()) ;
    Scalar integrant_un (vecteur_un(1)) ;
    
    Vector vecteur_deux (mapping, CON, mapping.get_bvect_cart()) ;
    for (int i=1 ; i<=3 ; i++)
      vecteur_deux.set(i) = k_total(2, i) ;
    vecteur_deux.change_triad (mapping.get_bvect_spher()) ;
    Scalar integrant_deux (vecteur_deux(1)) ;
    
    // Multiplication by y and x :
    integrant_un.set_spectral_va() = integrant_un.get_spectral_va()
      .mult_st() ;
    integrant_un.set_spectral_va() = integrant_un.get_spectral_va()
      .mult_sp() ;
    
    integrant_deux.set_spectral_va() = integrant_deux.get_spectral_va()
      .mult_st() ;
    integrant_deux.set_spectral_va() = integrant_deux.get_spectral_va()
      .mult_cp() ;
    
    p_angu_mom->set(2) = mapping.integrale_surface_infini (-integrant_un
					 +integrant_deux) / (2*qpig) ;

  }

	*/
	
	if (p_angu_mom == 0x0) {	    // a new computation is requireed
    
	p_angu_mom = new Tbl(3) ; 
    p_angu_mom->annule_hard() ;	// fills the double array with zeros

	// Reference Cartesian vector basis of the Absolute frame
	Base_vect_cart bvect_ref(0.) ; 	// 0. = parallel to the Absolute frame
	
	for (int i=0; i<=1; i++) {	    // loop on the stars

		const Map& mp = et[i]->get_mp() ; 
		int nzm1 = mp.get_mg()->get_nzone() - 1 ; 
		
		// Azimuthal vector d/dphi 
		Vector vphi(mp, CON, bvect_ref) ; 		
		vphi.set(1) = - mp.ya ; 	// phi^X
		vphi.set(2) = mp.xa ; 
		vphi.set(3) = 0 ;  
		vphi.annule_domain(nzm1) ; 	//##  phi set to zero in the CED	
		vphi.std_spectral_base() ; 
	    vphi.change_triad(mp.get_bvect_cart()) ; 
		
		// Matter part
		// -----------
		const Scalar& ee = et[i]->get_ener_euler() ;  // E
		const Scalar& pp = et[i]->get_press() ;	// p
		const Scalar& psi4 = et[i]->get_psi4() ; // Psi^4
		Scalar rho = pow(psi4, double(2.5)) * (ee+pp) ; 
		rho.std_spectral_base() ;

		Vector jmom = rho * (et[i]->get_u_euler()) ; 
				
		const Metric& gtilde = et[i]->get_gtilde() ; 
		
		Vector vphi_cov = vphi.up_down(gtilde) ;
		
		Scalar integrand = contract(jmom, 0, vphi_cov, 0) ; 
		      
		p_angu_mom->set(2) += integrand.integrale() ;

		// Extrinsic curvature part (0 if IWM)
		// -----------------------------------
		
		const Sym_tensor& aij = et[i]->get_tkij_auto() ;
		
		Sym_tensor kill_phi = vphi_cov.ope_killing(gtilde) ; 
		
		rho = pow(psi4, double(1.5)) ;   //## check the power !
		rho.std_spectral_base() ;
		
		integrand = rho * contract(aij, 0, 1, kill_phi, 0, 1) ; 
		
		p_angu_mom->set(2) += integrand.integrale() / (4*qpig) ;
		
		
	}  // End of the loop on the stars

    }	// End of the case where a new computation was necessary
  
  	return *p_angu_mom ; 
  
}



		    //---------------------------------//
		    //		Total energy	       //
		    //---------------------------------//

double Binary::total_ener() const {
    /*
    if (p_total_ener == 0x0) {	    // a new computation is requireed
	
	p_total_ener = new double ; 
	    
	    *p_total_ener = mass_adm() - star1.mass_b() - star2.mass_b() ; 
	    
    }	// End of the case where a new computation was necessary
    
    */
    return *p_total_ener ; 
    
}


		    //---------------------------------//
		    //	 Error on the virial theorem   //
		    //---------------------------------//

double Binary::virial() const {
    
    if (p_virial == 0x0) {	    // a new computation is requireed
	
	p_virial = new double ; 
	    
	    *p_virial = 1. - mass_kom_vol() / mass_adm_vol() ; 
	    
	}
    
    return *p_virial ; 
    
}
